from ase import Atoms
from ase.optimize import *
from ase.constraints import *
from ase.io import write, read
from ase.io.trajectory import Trajectory
from ase.calculators.socketio import SocketIOCalculator
from ase.vibrations import Vibrations
from ase.calculators.mixing import SumCalculator
from ase.calculators.singlepoint import SinglePointCalculator
from molecule.molecule import Molecule, Group, ATOMTYPES
from sella import Sella, Constraints, IRC
from fireworks import *
from fireworks.core.rocket_launcher import rapidfire
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.queue.queue_launcher import rapidfire as rapidfirequeue
from fireworks.utilities.fw_serializers import load_object_from_file
from fireworks.core.fworker import FWorker
import fireworks.fw_config
from pysidt.sidt import *
import pynta.models
from pynta.transitionstate import get_unique_optimized_adsorbates,determine_TS_construction,get_unique_TS_structs,generate_constraints_harmonic_parameters,get_unique_TS_templates_site_pairings
from pynta.utils import *
from pynta.calculator import run_harmonically_forced, map_harmonically_forced, add_sella_constraint
from pynta.mol import *
from pynta.coveragedependence import get_unstable_pairs, mol_to_atoms, get_configs_for_calculation, get_cov_energies_configs_concern_tree, get_configurations, get_unique_adsorbate_admols, train_sidt_cov_dep_regressor, process_calculation
from pynta.geometricanalysis import *
from pynta.adsorbate import construct_initial_guess_files
from pynta.postprocessing import postprocess, write_rmg_libraries
import numpy as np
import json
import copy
from copy import deepcopy
import sys
import shutil
import time
import logging
import signal
from contextlib import contextmanager
from copy import deepcopy
from joblib import Parallel, delayed
import os 

def finalize_calc(atoms):
    """Cleanly stop a persistent calculator process (e.g. VaspInteractive keeps a live VASP
    process alive across ionic steps and must be finalized once the run is done). This is a
    no-op for calculators without a finalize() method, so it is safe to call for any software."""
    calc = getattr(atoms, "calc", None)
    if calc is not None and hasattr(calc, "finalize"):
        try:
            calc.finalize()
        except Exception:
            pass

def check_socket_software(software_name,socket):
    """Guard against incompatible calculator/socket combinations. VaspInteractive keeps the VASP
    process alive across steps by itself and must not be wrapped in pynta's SocketIOCalculator
    (which keeps the process alive by a different mechanism). main.py raises the same error at
    workflow-construction time; this is a defensive check for fireworks built directly."""
    if software_name == "VaspInteractive" and socket:
        raise ValueError(
            "software='VaspInteractive' is incompatible with socket=True: both keep the "
            "calculator process alive across steps, but by different mechanisms. "
            "Use VaspInteractive with socket=False (set nsw in software_kwargs >= the "
            "optimizer's max steps), or use software='Vasp' with socket=True."
        )

def enforce_vasp_pbc(atoms, software_name):
    """VASP is a plane-wave code and requires full 3D periodicity. Force
    pbc=(True,True,True) for Vasp/VaspInteractive regardless of the pbc stored
    in the geometry or configured on the Pynta object (e.g. a slab run set up
    with pbc=(T,T,F) that is later refined with VASP). No-op for other software."""
    if isinstance(software_name, str) and software_name.lower() in ("vasp", "vaspinteractive"):
        atoms.pbc = (True, True, True)
    return atoms

class OptimizationTask(FiretaskBase):
    def run_task(self, fw_spec):
        raise NotImplementedError

class EnergyTask(FiretaskBase):
    def run_task(self, fw_spec):
        raise NotImplementedError

class VibrationTask(FiretaskBase):
    def run_task(self, fw_spec):
        raise NotImplementedError

class CollectTask(FiretaskBase):
    def run_task(self, fw_spec):
        raise NotImplementedError

@explicit_serialize
class DoNothingTask(FiretaskBase):
    def run_task(self, fw_spec):
        return FWAction()

def optimize_firework(xyz,software,label,opt_method=None,sella=None,socket=False,order=0,software_kwargs={},opt_kwargs={},
                      run_kwargs={},constraints=[],parents=[],out_path=None,time_limit_hrs=np.inf,fmaxhard=0.0,ignore_errors=False,
                      target_site_num=None,metal=None,facet=None,priority=1,allow_fizzled_parents=False):
    d = {"xyz" : xyz, "software" : software,"label" : label}
    if opt_method: d["opt_method"] = opt_method
    if software_kwargs: d["software_kwargs"] = software_kwargs
    if opt_kwargs: d["opt_kwargs"] = opt_kwargs
    if run_kwargs: d["run_kwargs"] = run_kwargs
    if constraints: d["constraints"] = constraints
    sella = True if sella or (sella is None and order != 0) else False
    if not sella: assert order == 0 #without Sella only minization is possible
    d["order"] = order
    d["sella"] = sella
    d["socket"] = socket
    d["time_limit_hrs"] = time_limit_hrs
    d["fmaxhard"] = fmaxhard
    d["ignore_errors"] = ignore_errors
    d["target_site_num"] = target_site_num
    d["metal"] = metal
    d["facet"] = facet
    t1 = MolecularOptimizationTask(d)
    directory = os.path.dirname(xyz)
    if out_path is None: out_path = os.path.join(directory,label+".xyz")
    t2 = FileTransferTask({'files': [
        {'src': label+'.xyz', 'dest': out_path}, 
        {'src': label+'.traj', 'dest': os.path.join(directory, label+".traj")},
    ],
    'mode': 'copy', 'ignore_errors': ignore_errors})


    return Firework([t1,t2],parents=parents,name=label+"opt",spec={"_allow_fizzled_parents": allow_fizzled_parents,"_priority": priority})

@explicit_serialize
class MolecularOptimizationTask(OptimizationTask):
    required_params = ["software","label"]
    optional_params = ["software_kwargs","opt_method",
        "opt_kwargs","run_kwargs", "constraints","sella","order","socket","time_limit_hrs","fmaxhard","ignore_errors","target_site_num","metal","facet"]
    def run_task(self, fw_spec):
        errors = []
        software_kwargs = deepcopy(self["software_kwargs"]) if "software_kwargs" in self.keys() else dict()
        socket = self["socket"] if "socket" in self.keys() else False
        check_socket_software(self["software"],socket)
        if socket:
            unixsocket = "ase_"+self["software"].lower()+"_"+self["label"]+"_"+self["xyz"].replace("/","_").replace(".","_")
            socket_address = os.path.join("/tmp","ipi_"+unixsocket)
            if "command" in software_kwargs.keys() and "{unixsocket}" in software_kwargs["command"]:
                software_kwargs["command"] = software_kwargs["command"].format(unixsocket=unixsocket)

        software = to_ase_software(self["software"],software_kwargs)
       
        opt_kwargs = deepcopy(self["opt_kwargs"]) if "opt_kwargs" in self.keys() else dict()
        opt_method = name_to_ase_opt(self["opt_method"]) if "opt_method" in self.keys() else BFGS
        run_kwargs = deepcopy(self["run_kwargs"]) if "run_kwargs" in self.keys() else dict()
        sella = self["sella"] if "sella" in self.keys() else False
        order = self["order"] if "order" in self.keys() else 0
        time_limit_hrs = self["time_limit_hrs"] if "time_limit_hrs" in self.keys() else np.inf
        fmaxhard = self["fmaxhard"] if "fmaxhard" in self.keys() else 0.0
        target_site_num = self["target_site_num"] if "target_site_num" in self.keys() else None
        metal = self["metal"] if "metal" in self.keys() else None
        facet = self["facet"] if "facet" in self.keys() else None
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else False

        label = self["label"]
        xyz = self['xyz']
        suffix = os.path.split(xyz)[-1].split(".")[-1]

        try:
            if suffix == "xyz":
                sp = read(xyz)
            elif suffix == "traj": #take last point on trajectory
                sp = Trajectory(xyz)[-1]
            else: #assume xyz
                sp = read(xyz)
        except Exception as e:
            if not ignore_errors:
                raise e
            else:
                errors.append(e)
                return FWAction(stored_data={"error": errors,"converged": False})

        if socket and os.path.exists(socket_address):
            os.unlink(socket_address)

        enforce_vasp_pbc(sp,self["software"])
        sp.calc = SocketIOCalculator(software,log=sys.stdout,unixsocket=unixsocket) if socket else software

        constraints = deepcopy(self["constraints"]) if "constraints" in self.keys() else []

        if not sella:
            assert order == 0
            out_constraints = []
            for c in constraints:
                if isinstance(c,dict):
                    constraint = construct_constraint(c)
                    out_constraints.append(constraint)
                elif c == "freeze half slab":
                    out_constraints.append(FixAtoms([
                        atom.index for atom in sp if atom.position[2] < sp.cell[2, 2] / 2.
                    ]))
                elif c.split()[0] == "freeze" and c.split()[1] == "all": #ex: "freeze all Cu"
                    sym = c.split()[2]
                    out_constraints.append(FixAtoms(
                        indices=[atom.index for atom in sp if atom.symbol == sym]
                        ))
                elif c.split()[0] == "freeze" and c.split()[1] == "up" and c.split()[2] == "to":
                    n = int(c.split()[3])
                    out_constraints.append(FixAtoms(
                        indices=list(range(n))
                        ))
                
            sp.set_constraint(out_constraints)
            opt_kwargs["trajectory"] = label+".traj"

            opt = opt_method(sp,**opt_kwargs)

            try:
                if np.isinf(time_limit_hrs):
                    opt.run(**run_kwargs)
                else:
                    with limit_time(time_limit_hrs*60.0*60.0):
                        opt.run(**run_kwargs)
            except Exception as e:
                if not ignore_errors:
                    raise e
                else:
                    errors.append(e)
                # if socket: #auto debugging example
                #     import logging
                #     logging.error(e)
                #     errors.append(e)
                #     sp.calc.close()
                #     fw = restart_opt_firework(self,fw_spec["_tasks"])
                #     return FWAction(stored_data={"debugged_error": errors},detours=[fw])
                # else:
                #     if not ignore_errors:
                #         raise e
                #     else:
                #         errors.append(e)

        else:
            out_constraints = []
            for c in constraints:
                if isinstance(c,dict):
                    constraint = construct_constraint(c)
                    out_constraints.append(constraint)
                elif c == "freeze half slab":
                    out_constraints.append(FixAtoms([
                        atom.index for atom in sp if atom.position[2] < sp.cell[2, 2] / 2.
                    ]))
                elif c.split()[0] == "freeze" and c.split()[1] == "all": #ex: "freeze all Cu"
                    sym = c.split()[2]
                    out_constraints.append(FixAtoms(
                        indices=[atom.index for atom in sp if atom.symbol == sym]
                        ))
                elif c.split()[0] == "freeze" and c.split()[1] == "up" and c.split()[2] == "to":
                    n = int(c.split()[3])
                    out_constraints.append(FixAtoms(
                        indices=list(range(n))
                        ))

            sp.set_constraint(out_constraints)
            
            opt = Sella(sp,trajectory=label+".traj",order=order)
            try:
                if np.isinf(time_limit_hrs):
                    opt.run(**run_kwargs)
                else:
                    with limit_time(time_limit_hrs*60.0*60.0):
                        opt.run(**run_kwargs)
            except Exception as e:
                if not ignore_errors:
                    raise e
                else:
                    errors.append(e)

        if socket:
            try:
                sp.calc.close()
            except Exception as e:
                if self["software"] == "Espresso":
                    pass #Espresso tends to error even after socket calculations finish correctly
                else:
                    if not ignore_errors:
                        raise e
                    else:
                        errors.append(e)

        finalize_calc(sp)  #stop a persistent calc (e.g. VaspInteractive) now the optimization is done; results are cached

        try:
            converged = opt.converged()
        except TypeError:
            converged = opt.converged(opt.atoms.get_forces().ravel())

        if not converged: #optimization has converged
            fmax = np.inf
            try:
                fmax = get_fmax(sp)
            except:
                pass
            try:
                tr = Trajectory(label+".traj")
                fmaxmin = np.inf
                minind = len(tr)
                for i,spt in enumerate(tr):
                    if "freeze half slab" in constraints:
                        spt.set_constraint(FixAtoms([
                                        atom.index for atom in spt if atom.position[2] < spt.cell[2, 2] / 2.
                                    ]))
                    fmaxt = get_fmax(spt)
                    if fmaxt < fmaxmin:
                        minind = i
                        fmaxmin = fmaxt

                if fmaxmin < fmax:
                    fmax = fmaxmin
                    sp = tr[minind]
            except:
                pass
            if fmax < fmaxhard:
                converged = True
            else:
                e = ValueError("Did not converge fmax below fmaxhard: {0} > {1}".format(fmax,fmaxhard))
                if not ignore_errors:
                    raise e
                else:
                    errors.append(e)

        #In principle this was a good idea, but occupied site detection at defaults fails more often than the problem this solves occurs
        # if converged and target_site_num: #optimization converged to correct structure, for now just check has correct number of occupied sites
        #     cas = SlabAdsorptionSites(sp,facet,allow_6fold=False,composition_effect=False,
        #                     label_sites=True,
        #                     surrogate_metal=metal)
        #     adcov = SlabAdsorbateCoverage(sp,adsorption_sites=cas)
        #     sites = adcov.get_sites()
        #     occ = [site for site in sites if site["occupied"]]
        #     if target_site_num != len(occ):
        #         converged = False
        #         e = StructureError
        #         if not ignore_errors:
        #             raise e
        #         else:
        #             errors.append(e)

        if converged:
            if self["software"] == "XTB" and "initial_charges" in sp.arrays.keys():
                del sp.arrays["initial_charges"]
            
            # Cache the energy before detaching calculator
            energy = sp.get_potential_energy()
            forces = sp.get_forces() if hasattr(sp.calc, 'get_forces') else None
            
            sp_to_write = sp.copy()
            sp_to_write.calc = SinglePointCalculator(sp_to_write, energy=energy, forces=forces)
            
            for key in list(sp_to_write.arrays.keys()):
                if key not in ('positions', 'numbers') and len(sp_to_write.arrays[key]) != len(sp_to_write):
                    del sp_to_write.arrays[key]
            write(label+".xyz", sp_to_write)

        return FWAction(stored_data={"error": errors,"converged": converged})

@explicit_serialize
class MolecularOptimizationFailTask(OptimizationTask):
    required_params = ["software","label"]
    optional_params = ["software_kwargs","opt_method",
        "opt_kwargs","run_kwargs"]
    def run_task(self, fw_spec):
        print(fw_spec)
        software_kwargs = deepcopy(self["software_kwargs"]) if "software_kwargs" in self.keys() else dict()
        
        software = to_ase_software(self["software"],software_kwargs)

        opt_kwargs = deepcopy(self["opt_kwargs"]) if "opt_kwargs" in self.keys() else dict()
        opt_method = name_to_ase_opt(self["opt_method"]) if "opt_method" in self.keys() else BFGS
        run_kwargs = deepcopy(self["run_kwargs"]) if "run_kwargs" in self.keys() else dict()

        label = self["label"]
        xyz = self['xyz']
        sp = read(xyz)

        enforce_vasp_pbc(sp,self["software"])
        sp.calc = software
        opt = opt_method(sp,trajectory=label+".traj")
        opt.run(fmax=0.02,steps=2)
        finalize_calc(sp)

        if not opt.converged():
            fw = restart_opt_firework(self,fw_spec["_tasks"])
            return FWAction(detours=[fw])

        write(label+".xyz", sp)

        return FWAction()

def energy_firework(xyz,software,label,software_kwargs={},parents=[],out_path=None,ignore_errors=False):
    d = {"xyz" : xyz, "software" : software, "label" : label}
    if software_kwargs: d["software_kwargs"] = software_kwargs
    d["ignore_errors"] = ignore_errors
    t1 = MolecularEnergyTask(d)
    directory = os.path.dirname(xyz)
    if out_path is None: out_path = os.path.join(directory,label+"_energy.json")
    t2 = FileTransferTask({'files': [{'src': label+'_energy.json', 'dest': out_path}], 'mode': 'copy', 'ignore_errors': ignore_errors})
    return Firework([t1,t2],parents=parents,name=label+"energy")

@explicit_serialize
class MolecularEnergyTask(EnergyTask):
    required_params = ["software","label"]
    optional_params = ["software_kwargs","energy_kwargs","ignore_errors"]
    def run_task(self, fw_spec):
        xyz = self['xyz']
        software = to_ase_software(self["software"],software_kwargs)
        label = self["label"]
        software_kwargs = deepcopy(self["software_kwargs"]) if "software_kwargs" in self.keys() else dict()
        energy_kwargs = deepcopy(self["energy_kwargs"]) if "energy_kwargs" in self.keys() else dict()
        ignore_errors = deepcopy(self["ignore_errors"]) if "ignore_errors" in self.keys() else False

        sp = None
        try:
            sp = read(xyz)
            enforce_vasp_pbc(sp,self["software"])
            sp.calc = software(**software_kwargs)
            en = sp.get_potential_energy(**energy_kwargs)
            with open(label+'_energy.json', 'a') as file:
                wr = json.dump(en,file)
        except Exception as e:
            if not ignore_errors:
                raise e
            else:
                return FWAction(stored_data={"error": e}, exit=True)
        finally:
            finalize_calc(sp)

        return FWAction()

def vibrations_firework(xyz,software,label,software_kwargs={},parents=[],out_path=None,constraints=[],socket=False,ignore_errors=False):
    d = {"xyz" : xyz, "software" : software, "label" : label, "socket": socket}
    if software_kwargs: d["software_kwargs"] = software_kwargs
    if constraints: d["constraints"] = constraints
    d["ignore_errors"] = ignore_errors
    t1 = MolecularVibrationsTask(d)
    directory = os.path.dirname(xyz)
    if out_path is None:
        out_path = directory
    elif "." in os.path.split(out_path)[1]:
        out_path = os.path.split(out_path)[0]
    t2 = FileTransferTask({'files': [{'src': label+'_vib.json', 'dest': os.path.join(out_path,label+'_vib.json')},
        {'src':'vib.0.traj','dest': os.path.join(out_path,"vib.0.traj")},],
        'mode': 'copy', 'ignore_errors': ignore_errors})
    t3 = FileTransferTask({'files': [{'src':'vib','dest':os.path.join(out_path,"vib")}],
        'mode': 'copytree', 'ignore_errors': ignore_errors})
    return Firework([t1,t2,t3],parents=parents,name=label+"vib")

@explicit_serialize
class MolecularVibrationsTask(VibrationTask):
    required_params = ["software","label"]
    optional_params = ["software_kwargs","constraints","ignore_errors","socket"]
    def run_task(self, fw_spec):
        indices = None
        xyz = self['xyz']
        label = self["label"]
        software_kwargs = deepcopy(self["software_kwargs"]) if "software_kwargs" in self.keys() else dict()
        software = to_ase_software(self["software"],software_kwargs)
        ignore_errors = deepcopy(self["ignore_errors"]) if "ignore_errors" in self.keys() else False
        socket = self["socket"] if "socket" in self.keys() else False
        check_socket_software(self["software"],socket)
        if socket:
            unixsocket = "ase_"+self["software"].lower()+"_"+self["label"]+"_"+self["xyz"].replace("/","_").replace(".","_")
            socket_address = os.path.join("/tmp","ipi_"+unixsocket)
            if "command" in software_kwargs.keys() and "{unixsocket}" in software_kwargs["command"]:
                software_kwargs["command"] = software_kwargs["command"].format(unixsocket=unixsocket)

        if socket and os.path.exists(socket_address):
            os.unlink(socket_address)

        sp = None
        try:
            sp = read(xyz)
            enforce_vasp_pbc(sp,self["software"])
            sp.calc = SocketIOCalculator(software,log=sys.stdout,unixsocket=unixsocket) if socket else software

            constraints = deepcopy(self["constraints"]) if "constraints" in self.keys() else []
            out_constraints = []
            for c in constraints:
                if isinstance(c,dict):
                    constraint = construct_constraint(c)
                    out_constraints.append(constraint)
                elif c == "freeze half slab":
                    out_constraints.append(FixAtoms([
                        atom.index for atom in sp if atom.position[2] < sp.cell[2, 2] / 2.
                    ]))
                    indices = [atom.index for atom in sp if atom.position[2] > sp.cell[2, 2] / 2.]
                elif c.split()[0] == "freeze" and c.split()[1] == "all": #ex: "freeze all Cu"
                    sym = c.split()[2]
                    out_constraints.append(FixAtoms(
                        indices=[atom.index for atom in sp if atom.symbol == sym]
                        ))
                    indices = [atom.index for atom in sp if atom.symbol != sym]
                elif c.split()[0] == "freeze" and c.split()[1] == "up" and c.split()[2] == "to":
                    n = int(c.split()[3])
                    out_constraints.append(FixAtoms(
                        indices=list(range(n))
                        ))
                    indices = [ i for i in range(len(sp)) if i >= n]
                else:
                    raise ValueError
                
            sp.set_constraint(out_constraints)
            
            vib = Vibrations(sp,indices=indices)
            vib.run()

            if socket:
                try:
                    sp.calc.close()
                except Exception as e:
                    if self["software"] == "Espresso":
                        pass #Espresso tends to error even after socket calculations finish correctly
                    else:
                        if not ignore_errors:
                            raise e
                        else:
                            errors.append(e)

            vib.write_mode(n=0)

            vibdata = vib.get_vibrations()
            vibdict = {"hessian": vibdata.get_hessian_2d().tolist(),
                    "frequencies": [str(x) for x in vibdata.get_frequencies().tolist()]}
            with open(label+"_vib.json","w") as fout:
                json.dump(vibdict,fout)

        except Exception as e:
            if not ignore_errors:
                raise e
            else:
                return FWAction(stored_data={"error": e}, exit=True)
        finally:
            finalize_calc(sp)

        return FWAction()

@explicit_serialize
class MolecularAdsorbateEstimate(FiretaskBase):
    required_params = ["mol","mol_name","slab_path","path","metal","facet","sites","site_adjacency","repeats",
                        "single_site_bond_params_lists", "single_sites_lists",
                        "double_site_bond_params_lists","double_sites_lists",
                        "vib_obj_dict","nslab","Eharmtol","Eharmfiltertol","Nharmmin","pbc",
                        "harm_f_software", "harm_f_software_kwargs","opt_software","opt_software_kwargs","opt_constraints","vib_constraints",
                        "fmaxopt","socket"]
    optional_params = ["out_path","spawn_jobs","nprocs","postprocess"]
    def run_task(self, fw_spec):
        spawn_jobs = self["spawn_jobs"] if "spawn_jobs" in self.keys() else False
        nprocs = self["nprocs"] if "nprocs" in self.keys() else 1
        path = self["path"]
        repeats = self["repeats"]
        
        mol = Molecule().from_adjacency_list(self["mol"])
        mol_name = self["mol_name"]
        
        opt_software = self["opt_software"]
        opt_software_kwargs = self["opt_software_kwargs"]
        opt_constraints = self["opt_constraints"]
        vib_constraints = self["vib_constraints"]
        
        single_site_bond_params_lists = self["single_site_bond_params_lists"]
        single_sites_lists = self["single_sites_lists"]
        for q in single_sites_lists:
            for s in q:
                s["position"] = np.array(s["position"])
                s["normal"] = np.array(s["normal"])
        double_site_bond_params_lists = self["double_site_bond_params_lists"]
        double_sites_lists = self["double_sites_lists"]
        for q in double_sites_lists:
            for s in q:
                s["position"] = np.array(s["position"])
                s["normal"] = np.array(s["normal"])
            
        metal = self["metal"]
        facet = self["facet"]
        nslab = self["nslab"]
        pbc = self["pbc"]
        sites = self["sites"]
        for s in sites:
            s["position"] = np.array(s["position"])
            s["normal"] = np.array(s["normal"])
        
        site_adjacency = {int(k):[int(x) for x in v] for k,v in self["site_adjacency"].items()}
        Nharmmin = self["Nharmmin"]
        Eharmtol = self["Eharmtol"]
        Eharmfiltertol = self["Eharmfiltertol"]
        slab_path = self["slab_path"]
        slab = read(slab_path)
        harm_f_software = self["harm_f_software"]
        harm_f_software_kwargs = self["harm_f_software_kwargs"]
        fmaxopt = self["fmaxopt"]
        socket = self["socket"]
        postprocess = self["postprocess"] if "postprocess" in self.keys() else False
        
        xyzs = construct_initial_guess_files(mol,mol_name,path,slab,metal,
                               single_site_bond_params_lists,single_sites_lists,double_site_bond_params_lists,double_sites_lists,
                               Eharmtol,Eharmfiltertol,Nharmmin,sites,site_adjacency,pbc,nslab,harm_f_software,harm_f_software_kwargs,nprocs)
        
        if spawn_jobs:
            optfws = []
            optfws2 = []
            previbxyzs = []
            for xyz in xyzs:
                prefix = os.path.split(os.path.split(xyz)[0])[1]
                previbxyz = os.path.join(os.path.split(xyz)[0],prefix+".xyz")
                previbxyzs.append(previbxyz)
                fwopt = optimize_firework(os.path.join(path,"Adsorbates",mol_name,prefix,prefix+"_init.xyz"),
                    opt_software,"weakopt_"+prefix,socket=socket,
                    opt_method="MDMin",opt_kwargs={'dt': 0.05},software_kwargs=opt_software_kwargs,
                    run_kwargs={"fmax" : 0.5, "steps" : 70},parents=[],constraints=opt_constraints,
                    ignore_errors=True, metal=metal, facet=facet, priority=3.5)
                fwopt2 = optimize_firework(os.path.join(path,"Adsorbates",mol_name,str(prefix),"weakopt_"+prefix+".xyz"),
                    opt_software,prefix,socket=socket,
                    opt_method="QuasiNewton",software_kwargs=opt_software_kwargs,
                    run_kwargs={"fmax" : fmaxopt, "steps" : 70},parents=[fwopt],constraints=opt_constraints,
                    ignore_errors=True, metal=metal, facet=facet, priority=3,
                    allow_fizzled_parents=True)
                optfws.append(fwopt)
                optfws2.append(fwopt2)

                vib_obj_dict = {"software": opt_software, "label": mol_name, "software_kwargs": opt_software_kwargs,
                    "constraints": vib_constraints}

            cfw = collect_firework(previbxyzs,True,["vibrations_firework"],[vib_obj_dict],["vib.json"],[],parents=optfws2,label=mol_name,detour=False,
                                   postprocess=postprocess,metal=metal,facet=facet,sites=self["sites"],site_adjacency=self["site_adjacency"],repeats=self["repeats"],slab_path=slab_path)
            newwf = Workflow(optfws+optfws2+[cfw],name=mol_name+"_optvib")
            return FWAction(detours=newwf)
        else:
            return FWAction()


@explicit_serialize
class MolecularTSEstimate(FiretaskBase):
    required_params = ["rxn","ts_path","slab_path","adsorbates_path","rxns_file","path","metal","facet","sites","site_adjacency","repeats",
                        "name_to_adjlist_dict", "gratom_to_molecule_atom_maps",
                        "gratom_to_molecule_surface_atom_maps","irc_mode",
                        "vib_obj_dict","opt_obj_dict","nslab","Eharmtol","Eharmfiltertol","Nharmmin","max_num_hfsp_opts","max_dist_hfsp","surrogate_metal",
                        "harm_f_software", "harm_f_software_kwargs"]
    optional_params = ["out_path","spawn_jobs","nprocs","IRC_obj_dict","postprocess"]
    def run_task(self, fw_spec):
        out_path = self["out_path"] if "out_path" in self.keys() else ts_path
        spawn_jobs = self["spawn_jobs"] if "spawn_jobs" in self.keys() else False
        nprocs = self["nprocs"] if "nprocs" in self.keys() else 1
        IRC_obj_dict = self["IRC_obj_dict"] if "IRC_obj_dict" in self.keys() else None

        ts_path = self["ts_path"]
        rxn = self["rxn"]
        index = rxn["index"]
        metal = self["metal"]
        facet = self["facet"]
        nslab = self["nslab"]
        sites = self["sites"]
        repeats = self["repeats"]
        for s in sites:
            s["position"] = np.array(s["position"])
            s["normal"] = np.array(s["normal"])
        
        site_adjacency = {int(k):[int(x) for x in v] for k,v in self["site_adjacency"].items()}
        Eharmtol = self["Eharmtol"]
        Eharmfiltertol = self["Eharmfiltertol"]
        Nharmmin = self["Nharmmin"]
        max_num_hfsp_opts = self["max_num_hfsp_opts"]
        max_dist_hfsp = self["max_dist_hfsp"]
        slab_path = self["slab_path"]
        surrogate_metal = self["surrogate_metal"]
        slab = read(slab_path)
        irc_mode = self["irc_mode"]
        harm_f_software = self["harm_f_software"]
        harm_f_software_kwargs = self["harm_f_software_kwargs"]

        adsorbates_path = self["adsorbates_path"]
        postprocess = self["postprocess"] if "postprocess" in self.keys() else False 

        gratom_to_molecule_atom_maps = dict()
        gratom_to_molecule_surface_atom_maps = dict()
        for ad in os.listdir(adsorbates_path):
            if os.path.isdir(os.path.join(adsorbates_path,ad)):
                with open(os.path.join(adsorbates_path,ad,"info.json"),'r') as f:
                    info = json.load(f)
                    gratom_to_molecule_atom_maps[info["name"]] = {int(k):v for k,v in info["atom_to_molecule_atom_map"].items()}
                    gratom_to_molecule_surface_atom_maps[info["name"]] = {int(k):v for k,v in info["gratom_to_molecule_surface_atom_map"].items()}

        reactants = Molecule().from_adjacency_list(rxn["reactant"])
        reactants.multiplicity = reactants.get_radical_count() + 1
        products = Molecule().from_adjacency_list(rxn["product"])
        products.multiplicity = products.get_radical_count() + 1

        reactant_names = rxn["reactant_names"]
        product_names = rxn["product_names"]
        rxn_name = rxn["reaction"]
        family_name = rxn.get("reaction_family", None)

        mol_dict = {name: Molecule().from_adjacency_list(adj.replace("multiplicity -187","")) for name,adj in self["name_to_adjlist_dict"].items()}

        for sm,mol in mol_dict.items():
            mol.multiplicity = mol.get_radical_count() + 1

        reactant_mols = [mol_dict[name] for name in reactant_names]
        product_mols = [mol_dict[name] for name in product_names]

        adsorbates = get_unique_optimized_adsorbates(rxn,adsorbates_path,mol_dict,gratom_to_molecule_surface_atom_maps,sites,nslab)
        
        if any(len(v) == 0 for v in adsorbates.values()):
            raise ValueError("Missing reactant or product for reaction: {}".format({k:len(v) for k,v in adsorbates.items()}))

        forward,species_names = determine_TS_construction(reactant_names,
                    reactant_mols,product_names,product_mols)

        ordered_adsorbates = [adsorbates[name] for name in species_names]

        ts_dict = {"forward": forward, "name": rxn_name, "reactants": reactants.to_adjacency_list(), "products": products.to_adjacency_list(),
            "species_names": species_names, "nslab": nslab}

        num_surf_sites = [len(mol_dict[name].get_surface_sites()) for name in species_names]
        
        if forward:
            reverse_names = product_names
        else:
            temp = products
            products = reactants
            reactants = temp
            reverse_names = reactant_names

        mols = [mol_dict[name] for name in species_names]
        ts_dict["mols"] = [mol.to_adjacency_list() for mol in mols]
        ts_dict["ads_sizes"] = [ads_size(mol) for mol in mols]
        template_mol_map = get_template_mol_map(reactants,mols)
        ts_dict["template_mol_map"] = template_mol_map
        ts_dict["reverse_names"] = reverse_names
        ts_dict["molecule_to_atom_maps"] = [{value:key for key,value in gratom_to_molecule_atom_maps[name].items()} for name in species_names]
        ts_dict["family_name"] = family_name 
        
        with open(os.path.join(ts_path,"info.json"),'w') as f:
            json.dump(ts_dict,f)

        molecule_to_atom_maps = [{value:key for key,value in gratom_to_molecule_atom_maps[name].items()} for name in species_names]
        template_to_ase = {i:get_ase_index(i,template_mol_map,molecule_to_atom_maps,
                    nslab,[ads_size(mol) for mol in mols]) for i in range(len(reactants.atoms))}
        ase_to_mol_num = {}
        for tind,aind in template_to_ase.items():
            if aind:
                for i,mol_map in enumerate(template_mol_map):
                    if tind in mol_map.keys():
                        ase_to_mol_num[aind] = i
                        break

        tsstructs,tsmols,neighbor_sites,ninds = get_unique_TS_structs(adsorbates,species_names,slab,sites,site_adjacency,nslab,num_surf_sites,mol_dict,
                                 gratom_to_molecule_atom_maps,gratom_to_molecule_surface_atom_maps,
                                 facet,metal)

        print("number of TS guesses pre-empty-sites and multiple mappings:")
        print(len(tsstructs))

        nsites = len([a for a in reactants.atoms if a.is_surface_site() and len(a.bonds) == 0])
        
        unique_tsstructs,unique_tsmols,target_sites,label_site_mappings = get_unique_TS_templates_site_pairings(tsstructs,
                                        tsmols,reactants,products,nsites,slab,neighbor_sites,ninds,sites,nslab)
        
        print("number of unique TS guesses:")
        print(len(unique_tsstructs))
        
        tsstructs_out,constraint_lists,atom_bond_potential_lists,site_bond_potential_lists = generate_constraints_harmonic_parameters(
                                            unique_tsstructs,unique_tsmols,label_site_mappings,adsorbates,slab,reactants,
                                             products,rxn.get("reaction_family", None),template_reversed=(not forward),
                                            ordered_names=species_names,reverse_names=reverse_names,
                                            mol_dict=mol_dict,gratom_to_molecule_atom_maps=gratom_to_molecule_atom_maps,
                                            gratom_to_molecule_surface_atom_maps=gratom_to_molecule_surface_atom_maps,
                                            nslab=nslab,facet=facet,metal=metal,slab_sites=sites,site_adjacency=site_adjacency)

        print("number of TS guesses with empty sites and multiple mappings:")
        print(len(tsstructs_out))

        if max_num_hfsp_opts or max_dist_hfsp:
            sorted_inds = index_site_bond_potential_lists_by_site_distances(site_bond_potential_lists).tolist()
            if max_dist_hfsp:
                max_dists = [get_max_site_dist(site_bond_potential_lists[i]) for i in sorted_inds]
                dist_inds = [sorted_inds[k] for k,d in enumerate(max_dists) if d <= max_dist_hfsp]
            else:
                dist_inds = sorted_inds
            if max_num_hfsp_opts:
                count_inds = sorted_inds[:max_num_hfsp_opts]
            else:
                count_inds = sorted_inds
            inds = dist_inds if len(dist_inds) <= len(count_inds) else count_inds
            tsstructs_out = [tsstructs_out[ind] for ind in inds]
            atom_bond_potential_lists = [atom_bond_potential_lists[ind] for ind in inds]
            site_bond_potential_lists = [site_bond_potential_lists[ind] for ind in inds]
            constraint_lists = [constraint_lists[ind] for ind in inds]
            print("number of TS guesses after filtering by max distance between sites")
            print(len(tsstructs_out))

        inputs = [ (tsstructs_out[j],atom_bond_potential_lists[j],site_bond_potential_lists[j],nslab,constraint_lists[j],ts_path,j,molecule_to_atom_maps,ase_to_mol_num,harm_f_software,harm_f_software_kwargs) for j in range(len(tsstructs_out))]
        outputs = Parallel(n_jobs=nprocs)(delayed(map_harmonically_forced)(inp) for inp in inputs)

        xyzs = [output[2] for output in outputs if output[0]]
        Es = [output[1] for output in outputs if output[0]]

        xyzs,Es = filter_nonunique_TS_guess_indices(xyzs,Es) #remove identical guesses (that will just get filtered out later in the collect resulting in less guesses)

        Einds = np.argsort(np.array(Es))
        Emin = np.min(np.array(Es))
        xyzsout = []
        for Eind in Einds:
            if Es[Eind]/Emin < Eharmtol: #include all TSs with energy close to Emin
                xyzsout.append(xyzs[Eind])
            elif Es[Eind]/Emin > Eharmfiltertol: #if the energy is much larger than Emin skip it
                continue
            elif len(xyzsout) < Nharmmin: #if the energy isn't similar, but isn't much larger include the smallest until Nharmmin is reached
                xyzsout.append(xyzs[Eind])

        if spawn_jobs:
            if irc_mode == "relaxed" or irc_mode == "fixed":
                irc_obj_dict_forward = deepcopy(self["IRC_obj_dict"])
                irc_obj_dict_forward["forward"] = True
                irc_obj_dict_reverse = deepcopy(self["IRC_obj_dict"])
                irc_obj_dict_reverse["forward"] = False

                ctask = MolecularCollect({"xyzs":xyzsout,"check_symm":True,"fw_generators": ["optimize_firework",["vibrations_firework","IRC_firework","IRC_firework"]],
                    "fw_generator_dicts": [self["opt_obj_dict"],[self["vib_obj_dict"],irc_obj_dict_forward,irc_obj_dict_reverse]],
                    "out_names": ["opt.xyz",["vib.json","irc_forward.traj","irc_reverse.traj"]],"future_check_symms": [True,False], "label": "TS"+str(index)+"_"+rxn_name,
                    "postprocess": postprocess, "metal": metal, "facet": facet, "sites": self["sites"], "site_adjacency": self["site_adjacency"],
                    "slab_path": slab_path, "repeats": repeats})
                cfw = Firework([ctask],name="TS"+str(index)+"_"+rxn_name+"_collect",spec={"_allow_fizzled_parents": True, "_priority": 5})
                newwf = Workflow([cfw],name='rxn_'+str(index)+str(rxn_name))
                return FWAction(detours=newwf) #using detour allows us to inherit children from the original collect to the subsequent collects
                
            #if irc_mode == "skip":
            else:
                ctask = MolecularCollect({"xyzs":xyzsout,"check_symm":True,"fw_generators": ["optimize_firework",["vibrations_firework"]],
                        "fw_generator_dicts": [self["opt_obj_dict"],[self["vib_obj_dict"]]],
                        "out_names": ["opt.xyz",["vib.json"]],"future_check_symms": [True,False], "label": "TS"+str(index)+"_"+rxn_name, "repeats": repeats})
                cfw = Firework([ctask],name="TS"+str(index)+"_"+rxn_name+"_collect",spec={"_allow_fizzled_parents": True, "_priority": 5})
                newwf = Workflow([cfw],name='rxn_'+str(index)+str(rxn_name))
                return FWAction(detours=newwf) #using detour allows us to inherit children from the original collect to the subsequent collects
        else:
            return FWAction()

def collect_firework(xyzs,check_symm,fw_generators,fw_generator_dicts,out_names,future_check_symms,repeats,parents=[],label="",allow_fizzled_parents=False,
                     detour=True,postprocess=False,metal=None,facet=None,sites=None,site_adjacency=None,slab_path=None):
    task = MolecularCollect({"xyzs": xyzs, "check_symm": check_symm, "fw_generators": fw_generators,
        "fw_generator_dicts": fw_generator_dicts, "out_names": out_names, "future_check_symms": future_check_symms, "label": label, 
        "detour": detour, "postprocess": postprocess, "metal": metal, "facet": facet, "sites": sites, "site_adjacency": site_adjacency,
        "slab_path": slab_path, "repeats": repeats})
    return Firework([task],parents=parents,name=label+"collect",spec={"_allow_fizzled_parents": allow_fizzled_parents,"_priority": 5})

@explicit_serialize
class MolecularCollect(CollectTask):
    required_params = ["xyzs","check_symm","fw_generators","fw_generator_dicts","out_names","future_check_symms","label","repeats"]
    optional_params = ["detour","postprocess","metal","facet","sites","site_adjacency","slab_path"]
    def run_task(self, fw_spec):
        xyzs = [xyz for xyz in self["xyzs"] if os.path.exists(xyz)] #if the associated task errored a file may not be present
        if self["check_symm"]:
            xyzs = get_unique_sym(xyzs) #only unique structures
        if len(xyzs) == 0:
            raise ValueError("No xyzs to collect")

        fw_generators = self["fw_generators"]
        fw_generator_dicts = self["fw_generator_dicts"]
        out_names = self["out_names"]
        future_check_symms = self["future_check_symms"]
        detour = self["detour"] if "detour" in self.keys() else True
        postprocess = self["postprocess"] if "postprocess" in self.keys() else False 
        metal = self["metal"] if "metal" in self.keys() else None 
        facet = self["facet"] if "facet" in self.keys() else None 
        sites = self["sites"] if "sites" in self.keys() else None 
        site_adjacency = self["site_adjacency"] if "site_adjacency" in self.keys() else None 
        slab_path = self["slab_path"] if "slab_path" in self.keys() else None
        repeats = self["repeats"] if "repeats" in self.keys() else None
        
        for i in range(len(fw_generators)):
            if not isinstance(fw_generators[i],list):
                fw_generators[i] = [fw_generators[i]]
            if not isinstance(fw_generator_dicts[i],list):
                fw_generator_dicts[i] = [fw_generator_dicts[i]]
            if not isinstance(out_names[i],list):
                out_names[i] = [out_names[i]]


        fws = []
        out_xyzs = []
        for i,fw_generator in enumerate(fw_generators[0]):
            fw_generator = globals()[fw_generator]
            for xyz in xyzs:
                d = deepcopy(fw_generator_dicts[0][i])
                d["xyz"] = xyz
                d["out_path"] = os.path.join(os.path.split(xyz)[0],out_names[0][i])
                d["label"] = out_names[0][i]
                d["ignore_errors"] = True
                out_xyzs.append(d["out_path"])
                fw = fw_generator(**d)
                if not isinstance(fw,list):
                    fw = [fw]
                fws.extend(fw)

        if len(fw_generators) > 1:
            task = MolecularCollect({"xyzs": out_xyzs,"check_symm": future_check_symms[0],
                    "fw_generators": fw_generators[1:],"fw_generator_dicts": fw_generator_dicts[1:],
                    "out_names": out_names[1:],"future_check_symms": future_check_symms[1:],"label": self["label"], 
                    "postprocess": postprocess, "metal": metal, "facet": facet, "sites": sites, "site_adjacency": site_adjacency,
                    "slab_path": slab_path, "repeats": repeats})
            cfw = Firework([task],parents=fws,name=self["label"]+"collect",spec={"_allow_fizzled_parents":True,"_priority": 4})
            newwf = Workflow(fws+[cfw],name=self["label"]+"collect"+str(-len(self["fw_generators"])))
            if detour:
                return FWAction(detours=newwf) #using detour allows us to inherit children from the original collect to the subsequent collects
            else:
                return FWAction(additions=newwf)
        else: #last round 
            if postprocess:
                pfw = postprocessing_firework(os.path.split(os.path.split(xyzs[0])[0])[0],metal,facet,sites,site_adjacency,slab_path,self["label"],repeats,parents=fws,
                    priority=10, allow_fizzled_parents=True)
                newwf = Workflow(fws+[pfw],name=self["label"]+"collectfinal")
                if detour:
                    return FWAction(detours=newwf)
                else:
                    return FWAction(additions=newwf)
            else:
                if detour:
                    return FWAction(detours=fws)
                else:
                    return FWAction(additions=fws)

def postprocessing_firework(path,metal,facet,sites,site_adjacency,slab_path,label,repeats,parents=[],priority=10,allow_fizzled_parents=False):
    """
    generates firework that marks the species/TS corresponding to path as complete and then postprocesses everything marked as complete 
    in the Pynta run
    """
    d = {"path" : path, "metal": metal, "facet": facet, "sites": sites, "site_adjacency": site_adjacency, "slab_path": slab_path, "repeats": repeats}
    
    t1 = PostprocessingTask(d)
    return Firework([t1],parents=parents,name=label+"_postprocessing",spec={"_allow_fizzled_parents": allow_fizzled_parents,"_priority": priority})

@explicit_serialize
class PostprocessingTask(FiretaskBase):
    """
    firework that marks the species/TS corresponding to path as complete and then postprocesses everything marked as complete 
    in the Pynta run
    """
    required_params = ["path","metal","facet","sites","site_adjacency","slab_path","repeats"]
    optional_params = []
    def run_task(self, fw_spec):
        path = self["path"]
        metal = self["metal"]
        facet = self["facet"]
        sites = self["sites"]
        repeats = self["repeats"]
        for s in sites:
            s["position"] = np.array(s["position"])
            s["normal"] = np.array(s["normal"])
        
        site_adjacency = {int(k):[int(x) for x in v] for k,v in self["site_adjacency"].items()}
        slab_path = self["slab_path"]
        
        is_TS = os.path.split(path)[1][:2] == "TS"
        
        open(os.path.join(path,"complete.sgnl"),'a').close()
        
        time.sleep(1) #wait one second, in case of any synchronization issues
        
        if is_TS:
            pynta_path = os.path.split(path)[0]
        else:
            pynta_path = os.path.split(os.path.split(path)[0])[0]
            
        spc_dict,ts_dict,spc_dict_thermo = postprocess(pynta_path,metal,facet,sites,site_adjacency,slab_path=slab_path,repeats=repeats,check_finished=True)
        
        write_rmg_libraries(pynta_path,spc_dict,spc_dict_thermo,ts_dict,metal,facet)
    
def TSnudge_firework(xyz,label,forward_path=None,reverse_path=None,spawn_jobs=False,software=None,opt_method=None,sella=False,
        socket=False,software_kwargs={},opt_kwargs={},run_kwargs={},constraints=[],parents=[],out_path=None,ignore_errors=False):
        """
        xyz is really the vibrational output data file
        out_path is a dud variable here
        """
        task = MolecularTSNudge(vib_traj=xyz,label=label,forward_path=forward_path,reverse_path=reverse_path,spawn_jobs=spawn_jobs,software=software,
            opt_method=opt_method,sella=sella,socket=socket,software_kwargs=software_kwargs,opt_kwargs=opt_kwargs,run_kwargs=run_kwargs,
            constraints=constraints,ignore_errors=ignore_errors)
        fw = Firework([task],parents=[],name=label+"TSnudge")
        return fw

@explicit_serialize
class MolecularTSNudge(FiretaskBase):
    required_params = ["vib_traj","label"]
    optional_params = ["forward_path","reverse_path","spawn_jobs","software","opt_method","sella","socket",
            "software_kwargs", "opt_kwargs", "run_kwargs", "constraints", "ignore_errors"]
    def run_task(self, fw_spec):
        forward_path = self["forward_path"] if "forward_path" in self.keys() and self["forward_path"] else os.path.join(os.path.split(self["vib_traj"])[0],self["label"]+"_forward")
        reverse_path = self["reverse_path"] if "reverse_path" in self.keys() and self["reverse_path"] else os.path.join(os.path.split(self["vib_traj"])[0],self["label"]+"_reverse")
        spawn_jobs = self["spawn_jobs"] if "spawn_jobs" in self.keys() else False
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else False

        AfterTS.get_forward_and_reverse(
            self["vib_traj"],
            forward_path,
            reverse_path)

        if spawn_jobs:
            fwf = optimize_firework(forward_path+".xyz",self["software"],self["label"],opt_method=self["opt_method"],sella=self["sella"],
                socket=self["socket"],order=0,software_kwargs=self["software_kwargs"],opt_kwargs=self["opt_kwargs"],
                                  run_kwargs=self["run_kwargs"],constraints=self["constraints"],parents=[],
                                  out_path=os.path.join(os.path.split(forward_path)[0],self["label"]+"_forward_optimized.xyz"),
                                  ignore_errors=ignore_errors)
            fwr = optimize_firework(reverse_path+".xyz",self["software"],self["label"],opt_method=self["opt_method"],sella=self["sella"],
                socket=self["socket"],order=0,software_kwargs=self["software_kwargs"],opt_kwargs=self["opt_kwargs"],
                                  run_kwargs=self["run_kwargs"],constraints=self["constraints"],parents=[],
                                  out_path=os.path.join(os.path.split(reverse_path)[0], self["label"]+"_reverse_optimized.xyz"),
                                  ignore_errors=ignore_errors)
            return FWAction(additions=[fwf,fwr])
        else:
            return FWAction()

def IRC_firework(xyz,label,out_path=None,spawn_jobs=False,software=None,
        socket=False,software_kwargs={},opt_kwargs={},run_kwargs={},constraints=[],parents=[],ignore_errors=False,forward=True):
        if out_path is None: out_path = os.path.join(directory,label+"_irc.traj")
        t1 = MolecularIRC(xyz=xyz,label=label,software=software,
            socket=socket,software_kwargs=software_kwargs,opt_kwargs=opt_kwargs,run_kwargs=run_kwargs,
            constraints=constraints,ignore_errors=ignore_errors,forward=forward)
        t2 = FileTransferTask({'files': [{'src': label+'_irc.traj', 'dest': out_path}], 'mode': 'copy', 'ignore_errors' : ignore_errors})
        fw = Firework([t1,t2],parents=[],name=label+"_IRC",spec={"_priority": 0.5})
        return fw

@explicit_serialize
class MolecularIRC(FiretaskBase):
    required_params = ["xyz","label"]
    optional_params = ["software","socket",
            "software_kwargs", "opt_kwargs", "run_kwargs", "constraints", "ignore_errors", "forward"]
    def run_task(self, fw_spec):
        errors = []
        software_kwargs = deepcopy(self["software_kwargs"]) if "software_kwargs" in self.keys() else dict()
        socket = self["socket"] if "socket" in self.keys() else False
        check_socket_software(self["software"],socket)
        if socket:
            unixsocket = "ase_"+self["software"].lower()+"_"+self["label"]+"_"+self["xyz"].replace("/","_").replace(".","_")
            socket_address = os.path.join("/tmp","ipi_"+unixsocket)
            if "command" in software_kwargs.keys() and "{unixsocket}" in software_kwargs["command"]:
                software_kwargs["command"] = software_kwargs["command"].format(unixsocket=unixsocket)

        software = to_ase_software(self["software"],software_kwargs)

        opt_kwargs = deepcopy(self["opt_kwargs"]) if "opt_kwargs" in self.keys() else dict()
        run_kwargs = deepcopy(self["run_kwargs"]) if "run_kwargs" in self.keys() else dict()
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else False
        forward = self["forward"] if "forward" in self.keys() else False

        label = self["label"]
        xyz = self['xyz']
        suffix = os.path.split(xyz)[-1].split(".")[-1]

        try:
            if suffix == "xyz":
                sp = read(xyz)
            elif suffix == "traj": #take last point on trajectory
                sp = Trajectory(xyz)[-1]
            else: #assume xyz
                sp = read(xyz)
        except Exception as e:
            if not ignore_errors:
                raise e
            else:
                errors.append(e)

        if socket and os.path.exists(socket_address):
            os.unlink(socket_address)

        enforce_vasp_pbc(sp,self["software"])
        sp.calc = SocketIOCalculator(software,log=sys.stdout,unixsocket=unixsocket) if socket else software

        constraints = deepcopy(self["constraints"]) if "constraints" in self.keys() else []
        
        out_constraints = []
        for c in constraints:
            if c == "freeze half slab":
                out_constraints.append(FixAtoms([
                    atom.index for atom in sp if atom.position[2] < sp.cell[2, 2] / 2.
                ]))
            elif c.split()[0] == "freeze" and c.split()[1] == "all": #ex: "freeze all Cu"
                sym = c.split()[2]
                out_constraints.append(FixAtoms(
                    indices=[atom.index for atom in sp if atom.symbol == sym]
                    ))
            elif c.split()[0] == "freeze" and c.split()[1] == "up" and c.split()[2] == "to":
                n = int(c.split()[3])
                out_constraints.append(FixAtoms(
                    indices=list(range(n))
                    ))
            else:
                raise ValueError("Could not interpret constraint: {}".format(c))
            
        sp.set_constraint(out_constraints)
        
        opt = IRC(sp,trajectory=label+"_irc.traj",dx=0.1,eta=1e-4,gamma=0.4)
        try:
            if forward:
                run_kwargs["direction"] = "forward"
                opt.run(**run_kwargs)
            else:
                run_kwargs["direction"] = "reverse"
                opt.run(**run_kwargs)
        except Exception as e:
            if not ignore_errors:
                raise e
            else:
                errors.append(e)

        if socket:
            try:
                sp.calc.close()
            except Exception as e:
                if self["software"] == "Espresso":
                    pass #Espresso tends to error even after socket calculations finish correctly
                else:
                    if not ignore_errors:
                        raise e
                    else:
                        errors.append(e)

        finalize_calc(sp)

        if not opt.converged():
            e = ValueError
            if not ignore_errors:
                raise e
            else:
                errors.append(e)

        if len(errors) == 0:
            pass
        else:
            return FWAction(stored_data={"error": errors})

        return FWAction()

def HFSP_firework(xyz,atom_bond_potentials,site_bond_potentials,nslab,constraints,molecule_to_atom_maps,ase_to_mol_num,
                      out_path=None,label="",parents=[],ignore_errors=False):
    d = {"xyz": xyz, "atom_bond_potentials": atom_bond_potentials, "site_bond_potentials": site_bond_potentials, 
         "nslab": nslab, "constraints": constraints, "molecule_to_atom_maps": molecule_to_atom_maps, "ase_to_mol_num": ase_to_mol_num,
        "label": label, "ignore_errors": ignore_errors}
    t1 = MolecularHFSP(d)
    directory = os.path.dirname(xyz)
    if out_path is None: out_path = os.path.join(directory,label+".xyz")
    t2 = FileTransferTask({'files': [{'src': label+'.xyz', 'dest': out_path}, {'src': "xtbharm.traj", 'dest': os.path.join(directory,label+".traj")}],
            'mode': 'copy', 'ignore_errors' : ignore_errors})
    return Firework([t1,t2],parents=parents,name=label+"HFSP")

@explicit_serialize
class MolecularHFSP(OptimizationTask):
    required_params = ["xyz","atom_bond_potentials","site_bond_potentials","nslab","constraints","molecule_to_atom_maps","ase_to_mol_num"]
    optional_params = ["label","ignore_errors","method"]

    def run_task(self, fw_spec):
        label = self["label"] if "label" in self.keys() else "xtb"
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else False
        method = self["method"] if "method" in self.keys() else "GFN1-xTB"
        
        atom_bond_potentials = self["atom_bond_potentials"]
        site_bond_potentials = self["site_bond_potentials"]
        nslab = self["nslab"]
        molecule_to_atom_maps = [{int(k):v for k,v in x.items()} for x in self["molecule_to_atom_maps"]]
        ase_to_mol_num = {int(k):v for k,v in self["ase_to_mol_num"].items()}
        constraints = self["constraints"]
        xyz = self['xyz']
        
        errors = []
        
        suffix = os.path.split(xyz)[-1].split(".")[-1]
        
        try:
            if suffix == "xyz":
                sp = read(xyz)
            elif suffix == "traj": #take last point on trajectory
                sp = Trajectory(xyz)[-1]
            else: #assume xyz
                sp = read(xyz)
        except Exception as e:
            if not ignore_errors:
                raise e
            else:
                errors.append(e)

        spout,Eharm,Fharm = run_harmonically_forced(sp,atom_bond_potentials,site_bond_potentials,nslab,
                    molecule_to_atom_maps=molecule_to_atom_maps,ase_to_mol_num=ase_to_mol_num,
                    method="GFN1-xTB",constraints=constraints)
        if spout:
            if "initial_charges" in sp.arrays.keys(): #avoid bug in ase
                del sp.arrays["initial_charges"]
            write(label+".xyz", spout)
            converged = True 
        else:
            converged = False
        
        return FWAction(stored_data={"error": errors,"converged": converged})

def calculate_configruation_energies_firework(admol_name,tree_file,path,coadname,coad_stable_sites,Nocc_isolated,
                                              coadmol_E_dict,concern_energy_tol=None,out_path=None,parents=[],iter=0,ignore_errors=False):
    d = {"admol_name": admol_name,"tree_file": tree_file,"path": path,"coadname": coadname,
         "coad_stable_sites": coad_stable_sites,"Nocc_isolated": Nocc_isolated,"coadmol_E_dict": coadmol_E_dict, 
         "concern_energy_tol": concern_energy_tol}
    t1 = CalculateConfigurationEnergiesTask(d)
    if out_path is None: 
        out_path_energy = os.path.join(os.path.split(tree_file)[0],"Ncoad_energy_"+admol_name+"_"+coadname+".json")
        out_path_config = os.path.join(os.path.split(tree_file)[0],"Ncoad_config_"+admol_name+"_"+coadname+".json")
        out_path_concern = os.path.join(os.path.split(tree_file)[0],"configs_of_concern_"+admol_name+"_"+coadname+".json")
    else:
        out_path_energy = os.path.join(out_path,"Ncoad_energy_"+admol_name+"_"+coadname+".json")
        out_path_config = os.path.join(out_path,"Ncoad_config_"+admol_name+"_"+coadname+".json")
        out_path_concern = os.path.join(out_path,"configs_of_concern_"+admol_name+"_"+coadname+".json")
    t2 = FileTransferTask({'files': [{'src': "Ncoad_energy_"+admol_name+"_"+coadname+".json", 'dest': out_path_energy},
                                     {'src': "Ncoad_config_"+admol_name+"_"+coadname+".json", 'dest': out_path_config},
                                     {'src': "configs_of_concern_"+admol_name+"_"+coadname+".json", 'dest': out_path_concern}], 'mode': 'copy', 'ignore_errors': ignore_errors})
    return Firework([t1,t2],parents=parents,name=admol_name+"_"+coadname+"_energies"+str(iter))

@explicit_serialize
class CalculateConfigurationEnergiesTask(FiretaskBase):
    required_params = ["admol_name","tree_file","path","coad_stable_sites","Nocc_isolated","coadmol_E_dict","coadname"]
    optional_params = ["concern_energy_tol","ignore_errors"]
    def run_task(self, fw_spec):
        admol_name = self['admol_name']
        tree_file = self["tree_file"]
        path= self["path"]
        coad_stable_sites = self["coad_stable_sites"]
        Nocc_isolated = self["Nocc_isolated"]
        coadname = self["coadname"]
        coadmol_E_dict = {Molecule().from_adjacency_list(k):v for k,v in self["coadmol_E_dict"].items()}
        concern_energy_tol = self["concern_energy_tol"] if "concern_energy_tol" in self.keys() else None
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else False
        
        try:
            nodes = read_nodes(tree_file)
            tree = MultiEvalSubgraphIsomorphicDecisionTreeRegressor([adsorbate_interaction_decomposition],
                                                        nodes=nodes)
            
            with open(os.path.join(path,"Configurations",admol_name+"_"+coadname+".json"),'r') as f:
                configs = [Molecule().from_adjacency_list(m,check_consistency=False) for m in json.load(f)]
            
            Ncoad_energy_dict,Ncoad_config_dict,configs_of_concern_admol = get_cov_energies_configs_concern_tree(tree, configs, coad_stable_sites, Nocc_isolated, concern_energy_tol, 
                                                coadmol_E_dict=coadmol_E_dict)
            with open("Ncoad_energy_"+admol_name+"_"+coadname+".json",'w') as f:
                json.dump(Ncoad_energy_dict,f)
            with open("Ncoad_config_"+admol_name+"_"+coadname+".json",'w') as f:
                json.dump(Ncoad_config_dict,f)
            with open("configs_of_concern_"+admol_name+"_"+coadname+".json",'w') as f:
                json.dump([tuple([v[0].to_adjacency_list(),v[1],v[2],v[3]]) for v in configs_of_concern_admol.values()],f)
                
        except Exception as e:
            if not ignore_errors:
                raise e
            else:
                return FWAction(stored_data={"error": e}, exit=True)

        return FWAction()

def extract_datum_firework(d, pynta_dir, slab_path, metal, facet, sites, site_adjacency,
                           ad_energy_dict_path, coadmol_E_dict_path, is_ad, parents=[], ignore_errors=True):
    spec = {
        "d": d,
        "pynta_dir": pynta_dir,
        "slab_path": slab_path,
        "metal": metal,
        "facet": facet,
        "sites": sites,
        "site_adjacency": {str(k): v for k,v in site_adjacency.items()} if isinstance(site_adjacency, dict) else site_adjacency,
        "ad_energy_dict_path": ad_energy_dict_path,
        "coadmol_E_dict_path": coadmol_E_dict_path,
        "is_ad": is_ad,
        "ignore_errors": ignore_errors,
    }
    t = ExtractDatumTask(spec)
    return Firework([t], parents=parents,
                    name="Extract Datum "+os.path.basename(d.rstrip("/")),
                    spec={"_priority": 1, "_allow_fizzled_parents": True})

@explicit_serialize
class ExtractDatumTask(FiretaskBase):
    required_params = ["d", "pynta_dir", "slab_path", "metal", "facet", "sites", "site_adjacency",
                       "ad_energy_dict_path", "coadmol_E_dict_path", "is_ad"]
    optional_params = ["ignore_errors"]

    def run_task(self, fw_spec):
        d = self["d"]
        pynta_dir = self["pynta_dir"]
        slab_path = self["slab_path"]
        metal = self["metal"]
        facet = self["facet"]
        is_ad = self["is_ad"]
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else True

        sites = []
        for site in self["sites"]:
            site = dict(site)
            site["normal"] = np.array(site["normal"])
            site["position"] = np.array(site["position"])
            site["indices"] = tuple(site["indices"])
            sites.append(site)
        site_adjacency = {int(k):[int(x) for x in v] for k,v in self["site_adjacency"].items()}

        slab = read(slab_path)

        with open(self["ad_energy_dict_path"]) as f:
            ad_energy_dict = {Molecule().from_adjacency_list(k,check_consistency=False): v
                              for k,v in json.load(f).items()}
        with open(self["coadmol_E_dict_path"]) as f:
            coadmol_E_dict = {Molecule().from_adjacency_list(k,check_consistency=False): v
                              for k,v in json.load(f).items()}

        out = {"valid": False}
        try:
            datum_E, datums_stability = process_calculation(d, ad_energy_dict, slab, metal, facet, sites,
                                                            site_adjacency, pynta_dir, coadmol_E_dict,
                                                            max_dist=3.0, rxn_alignment_min=0.7,
                                                            coad_disruption_tol=1.1,
                                                            out_file_name="out", init_file_name="init",
                                                            vib_file_name="vib_vib", is_ad=is_ad)
            out["valid"] = True
            if datum_E is not None:
                out["datum_E"] = {"mol": datum_E.mol.to_adjacency_list(), "value": float(datum_E.value)}
            else:
                out["datum_E"] = None
            out["datums_stability"] = [{"mol": ds.mol.to_adjacency_list(), "value": bool(ds.value)}
                                       for ds in datums_stability]
        except Exception as e:
            if not ignore_errors:
                raise
            out["valid"] = False
            out["reason"] = type(e).__name__

        with open(os.path.join(d, "datum.json"), "w") as f:
            json.dump(out, f)

        return FWAction()

def train_covdep_model_firework(path,admol_name_path_dict,admol_name_structure_dict,sites,site_adjacency,
                                pynta_dir, metal, facet, slab_path, calculation_directories, coadnames,
                                coad_stable_sites, software, software_kwargs, software_kwargs_TS, freeze_ind,
                                fmaxopt, parents=[],Ncalc_per_iter=6,iter=0,max_iters=6,concern_energy_tol=None,
                                ignore_errors=False, max_coadsorbates=None,sidt_isolated_delta_model=None,
                                sidt_covdep_delta_model=None,ts_frac=None,adsorbate_site_energy_cutoff=0.0):
    d = {"path": path, "admol_name_path_dict": admol_name_path_dict, "admol_name_structure_dict": {k : v.to_adjacency_list() for k,v in admol_name_structure_dict.items()},
         "sites": sites, "site_adjacency": {str(k):v for k,v in site_adjacency.items()}, "pynta_dir": pynta_dir, "metal": metal, "facet": facet, "slab_path": slab_path,
         "calculation_directories": calculation_directories, "coadnames": coadnames, "coad_stable_sites": coad_stable_sites,
        "Ncalc_per_iter": Ncalc_per_iter, "iter": iter, "max_iters": max_iters, "software": software, "software_kwargs": software_kwargs, "software_kwargs_TS": software_kwargs_TS, "freeze_ind": freeze_ind,
        "fmaxopt": fmaxopt, "concern_energy_tol": concern_energy_tol, "ignore_errors": ignore_errors, "max_coadsorbates": max_coadsorbates,
        "sidt_isolated_delta_model": sidt_isolated_delta_model, "sidt_covdep_delta_model": sidt_covdep_delta_model, "ts_frac": ts_frac,
        "adsorbate_site_energy_cutoff": adsorbate_site_energy_cutoff}
    t1 = TrainCovdepModelTask(d)
    return Firework([t1],parents=parents,name="Training Model "+str(iter),spec={"_allow_fizzled_parents":True, "_priority": 4})

@explicit_serialize
class TrainCovdepModelTask(FiretaskBase):
    required_params = ["path","admol_name_path_dict","admol_name_structure_dict","sites","site_adjacency", "pynta_dir", "metal", "facet",
                       "slab_path", "calculation_directories", "coadnames", "coad_stable_sites", "Ncalc_per_iter", "iter", "max_iters", "software", 
                       "software_kwargs", "software_kwargs_TS", "freeze_ind", "fmaxopt"]
    optional_params = ["concern_energy_tol","ignore_errors", "max_coadsorbates","sidt_isolated_delta_model","sidt_covdep_delta_model","ts_frac","adsorbate_site_energy_cutoff"]
    def run_task(self, fw_spec):
        path = self["path"]
        admol_name_path_dict = self["admol_name_path_dict"]
        admol_name_structure_dict = {k: Molecule().from_adjacency_list(v,check_consistency=False) for k,v in self["admol_name_structure_dict"].items()}
        max_coadsorbates = self["max_coadsorbates"] if "max_coadsorbates" in self.keys() else None
        adsorbate_site_energy_cutoff = self["adsorbate_site_energy_cutoff"] if "adsorbate_site_energy_cutoff" in self.keys() else 0.0
        sites = []
        for site in self["sites"]:
            site["normal"] = np.array(site["normal"])
            site["position"] = np.array(site["position"])
            site["indices"] = tuple(site["indices"])
            sites.append(site)
        site_adjacency = {int(k):[int(x) for x in v] for k,v in self["site_adjacency"].items()}
        pynta_dir = self["pynta_dir"]
        metal = self["metal"]
        facet = self["facet"]
        slab_path = self["slab_path"]
        calculation_directories = self["calculation_directories"]
        coadnames  = self["coadnames"]
        Ncalc_per_iter = self["Ncalc_per_iter"]
        iter = self["iter"]
        coad_stable_sites = self["coad_stable_sites"]
        concern_energy_tol = self["concern_energy_tol"] if "concern_energy_tol" in self.keys() else None
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else False
        software_kwargs = self["software_kwargs"]
        software = self["software"]
        freeze_ind = self["freeze_ind"]
        software_kwargs_TS = self["software_kwargs_TS"]
        fmaxopt = self["fmaxopt"]
        max_iters = self["max_iters"]
        
        if "sidt_isolated_delta_model" in self.keys() and self["sidt_isolated_delta_model"] is not None:
            nodes_isolated = read_nodes(self["sidt_isolated_delta_model"])
            def bond_decomposition_adsorbate(mol):
                pairs = []
                bonds = mol.get_all_edges()
                for bond in bonds:
                    if not (bond.atom1.is_surface_site() and bond.atom2.is_surface_site()):
                        pairs.append((mol.atoms.index(bond.atom1), mol.atoms.index(bond.atom2)))

                structs = []
                for pair in pairs:
                    m = mol.copy(deep=True)
                    for ind in pair:
                        m.atoms[ind].label = "*"
                    structs.append(m)

                return structs
            
            sidt_isolated_delta = MultiEvalSubgraphIsomorphicDecisionTreeRegressor(
                    bond_decomposition_adsorbate,
                    nodes=nodes_isolated,
                    root_group = Group().from_adjacency_list("""1 * R u0 px cx {2,[vdW,R,S,D,T,Q]}
                    2 * Rx u0 px cx {1,[vdW,R,S,D,T,Q]}"""),
                    r=[ATOMTYPES[x] for x in r_atoms],
                    r_bonds=[0, 0.05, 1, 2, 3, 4],
                    r_un=r_un,
                    r_site=r_site,
                    r_morph=r_morph,
                    fract_nodes_expand_per_iter=0.1,
            )
        else:
            sidt_isolated_delta = None 
        
        if "sidt_covdep_delta_model" in self.keys() and self["sidt_covdep_delta_model"] is not None:
            from pynta.coveragedependence import adsorbate_interaction_decomposition
            nodes_covdep = read_nodes(self["sidt_covdep_delta_model"])
            sidt_covdep_delta = MultiEvalSubgraphIsomorphicDecisionTreeRegressor([adsorbate_interaction_decomposition],
                                                        nodes=nodes_covdep,
                                                        r=[ATOMTYPES[x] for x in r_atoms],
                                                        r_bonds=[1,2,3,4,0.05],
                                                        r_un=r_un,
                                                        r_site=r_site,
                                                        r_morph=r_morph,
                                                        max_structures_to_generate_extensions=100,
                                                        fract_nodes_expand_per_iter=0.025,
                                                        iter_max=2,
                                                        iter_item_cap=100,
                                                        weigh_node_selection_by_occurrence=True,
                                                        )
        else:
            sidt_covdep_delta = None
        
        coads = [admol_name_structure_dict[coadname] for coadname in coadnames]
        coad_simples = {coadname: remove_slab(coads[i]) for i,coadname in enumerate(coadnames)}
        
        coad_paths = {coadname: os.path.join(pynta_dir,"Adsorbates",coadname) for coadname in coadnames}
        slab = read(slab_path)
        nslab = len(slab)

        r_site = list(set([a.site for admol in admol_name_structure_dict.values() for a in admol.atoms]))
        r_morph = list(set([a.morphology for admol in admol_name_structure_dict.values() for a in admol.atoms]))
        r_atoms = list(set([a.element.symbol for admol in admol_name_structure_dict.values() for a in admol.atoms]))
        r_un = list(set([a.radical_electrons for admol in admol_name_structure_dict.values() for a in admol.atoms]))
        r_lone_pairs = list(set([a.lone_pairs for admol in admol_name_structure_dict.values() for a in admol.atoms]))

        ad_energy_dict = get_lowest_adsorbate_energies(os.path.join(pynta_dir,"Adsorbates"),sidt_isolated_delta=sidt_isolated_delta)
        coad_Es = {coadname: get_adsorbate_energies(coad_paths[coadname],sidt_isolated_delta=sidt_isolated_delta)[0] for coadname in coadnames}

        # If coadmol_E_dict JSON caches exist AND Configurations/ exists, we can skip the expensive
        # generate_allowed_structure_site_structures + coadmol_E/stability_dict rebuild. The stability
        # dict is only consulted when Configurations/ is being built; the E dict is loaded from cache
        # (lookups in it use is_isomorphic, so JSON-rehydrated Molecule keys are equivalent).
        coadmol_E_dict_paths = {coadname: os.path.join(path, "coadmol_E_dict_" + coadname + ".json") for coadname in coadnames}
        configs_dir_exists = os.path.exists(os.path.join(path, "Configurations"))
        all_E_caches_exist = all(os.path.exists(p) for p in coadmol_E_dict_paths.values())

        coadmol_E_dicts = dict()
        coadmol_stability_dicts = dict()
        if configs_dir_exists and all_E_caches_exist:
            for coadname in coadnames:
                with open(coadmol_E_dict_paths[coadname]) as f:
                    coadmol_E_dicts[coadname] = {Molecule().from_adjacency_list(k, check_consistency=False): v
                                                 for k,v in json.load(f).items()}
            allowed_structure_site_structures = None
        else:
            allowed_structure_site_structures = generate_allowed_structure_site_structures(os.path.join(pynta_dir,"Adsorbates"),sites,site_adjacency,nslab,max_dist=np.inf)
            for coadname in coadnames:
                coad_path = coad_paths[coadname]
                coadmol_E_dict = dict()
                coadmol_stability_dict = dict()
                for p in os.listdir(coad_path):
                    if p == "info.json" or (p not in coad_Es[coadname].keys()):
                        continue
                    admol_init,neighbor_sites_init,ninds_init = generate_adsorbate_2D(read(os.path.join(coad_path,p,p+"_init.xyz")),sites,site_adjacency,nslab,max_dist=np.inf,allowed_structure_site_structures=allowed_structure_site_structures)
                    admol,neighbor_sites,ninds = generate_adsorbate_2D(read(os.path.join(coad_path,p,p+".xyz")),sites,site_adjacency,nslab,max_dist=np.inf,allowed_structure_site_structures=allowed_structure_site_structures)
                    out_struct = split_adsorbed_structures(admol,clear_site_info=False)[0]
                    out_struct_init = split_adsorbed_structures(admol_init,clear_site_info=False)[0]
                    coadmol_E_dict[out_struct] = coad_Es[coadname][p]
                    if admol_init.is_isomorphic(admol,save_order=True):
                        coadmol_stability_dict[out_struct_init] = True
                    else:
                        coadmol_stability_dict[out_struct_init] = False
                coadmol_E_dicts[coadname] = coadmol_E_dict
                coadmol_stability_dicts[coadname] = coadmol_stability_dict

        if not os.path.exists(os.path.join(path,"Configurations")):
            info_paths = {adname: os.path.join(os.path.split(os.path.split(p)[0])[0],"info.json") for adname,p in admol_name_path_dict.items()}
            imag_freq_paths = {adname: os.path.join(os.path.split(p)[0],"vib.json_vib.json") for adname,p in admol_name_path_dict.items()}
            unstable_pairs = get_unstable_pairs(os.path.join(path,'pairs'),
                                os.path.join(pynta_dir,"Adsorbates"),
                                sites,site_adjacency,nslab,max_dist=np.inf,show=False, infopath_dict=info_paths, imag_freq_path_dict=imag_freq_paths)
            
            os.makedirs(os.path.join(path,"Configurations"))
            for coadname in coadnames:
                for admol_name,admol in admol_name_structure_dict.items():
                    # Seed get_configurations with every unique stable adsorbate geometry
                    # (the same multi-site set the pairs explore), not just the single
                    # lowest-energy base. Without this the candidate set fixes the adsorbate
                    # at one site and only varies coadsorbates, so adsorbate-site
                    # arrangements that become favorable under coverage are never enumerated.
                    # TSs have a single geometry, so they keep their single base structure.
                    ad_dir = os.path.join(pynta_dir, "Adsorbates", admol_name)
                    if os.path.isdir(ad_dir):
                        base_admols = get_unique_adsorbate_admols(ad_dir, sites, site_adjacency, nslab,
                            allowed_structure_site_structures=allowed_structure_site_structures,
                            energy_cutoff=adsorbate_site_energy_cutoff)
                        if not base_admols:
                            base_admols = [admol]
                    else:
                        base_admols = [admol]

                    # concatenate across base structures; duplicate configs are harmless
                    # (the per-coverage minimum is unaffected) and avoids an O(n^2) dedup
                    configs = []
                    for base in base_admols:
                        configs.extend(get_configurations(base, coad_simples[coadname], coad_stable_sites[coadname],  coadmol_stability_dict=coadmol_stability_dicts[coadname], unstable_groups=unstable_pairs,
                            coadmol_E_dict=coadmol_E_dicts[coadname], max_coadsorbates=max_coadsorbates))

                    with open(os.path.join(path,"Configurations",admol_name+"_"+coadname+".json"),'w') as f:
                        json.dump([x.to_adjacency_list() for x in configs],f)
        
        if iter > 1:
            with open(os.path.join(path,"pairs_datums.json"),'r') as f:
                pairs_datums = [Datum(mol=Molecule().from_adjacency_list(d["mol"],check_consistency=False), value=d["value"]) for d in json.load(f)]
            with open(os.path.join(path,"Iterations",str(iter-1),"cumulative_sample_datums.json"),'r') as f:
                old_sample_datums = [Datum(mol=Molecule().from_adjacency_list(d["mol"],check_consistency=False), value=d["value"]) for d in json.load(f)]
        elif iter == 1:
            with open(os.path.join(path,"pairs_datums.json"),'r') as f:
                pairs_datums = [Datum(mol=Molecule().from_adjacency_list(d["mol"],check_consistency=False), value=d["value"]) for d in json.load(f)]
            old_sample_datums = [] 
        
        if iter == 0:
            computed_configs = []
        else:
            with open(os.path.join(path,"Iterations",str(iter-1),"computed_configurations.json"),'r') as f:
                computed_configs = [Molecule().from_adjacency_list(x,check_consistency=False) for x in json.load(f)]
        
        new_datums_E = []
        new_computed_configs = []
        for d in calculation_directories:
            with open(os.path.join(d,"info.json"),'r') as f:
                info = json.load(f)
            adjlist = info['adjlist']
            coadname = info["coadname"]
            init_config = Molecule().from_adjacency_list(adjlist,check_consistency=False)
            new_computed_configs.append(init_config)
            datum_path = os.path.join(d, "datum.json")
            if os.path.exists(datum_path):
                with open(datum_path) as f:
                    cached = json.load(f)
                if cached.get("valid") and cached.get("datum_E") is not None:
                    datum_E = Datum(Molecule().from_adjacency_list(cached["datum_E"]["mol"], check_consistency=False),
                                    cached["datum_E"]["value"])
                    new_datums_E.append(datum_E)
                    # original behavior appends datum_E.mol when init != out (the isomorphism check);
                    # process_calculation records that as a value=false entry in datums_stability,
                    # so checking for one is equivalent and avoids the isomorphism call entirely.
                    if any(not ds["value"] for ds in cached.get("datums_stability", [])):
                        new_computed_configs.append(datum_E.mol)
            else:
                datum_E,datums_stability = process_calculation(d,ad_energy_dict,slab,metal,facet,sites,site_adjacency,pynta_dir,coadmol_E_dicts[coadname],max_dist=3.0,rxn_alignment_min=0.7,
                    coad_disruption_tol=1.1,out_file_name="out",init_file_name="init",vib_file_name="vib_vib",is_ad=None,sidt_isolated_delta=sidt_isolated_delta,sidt_covdep_delta=sidt_covdep_delta)
                if datum_E:
                    new_datums_E.append(datum_E)
                if datum_E and not datum_E.mol.is_isomorphic(init_config,save_order=True):
                    new_computed_configs.append(datum_E.mol)
        
        if not os.path.exists(os.path.join(path,"Iterations",str(iter))):
            os.makedirs(os.path.join(path,"Iterations",str(iter)))
        
        with open(os.path.join(path,"Iterations",str(iter),"computed_configurations.json"),'w') as f:
            json.dump([x.to_adjacency_list() for x in computed_configs+new_computed_configs],f)
        
        if iter == 0:
            pairs_datums = new_datums_E
            with open(os.path.join(path,"pairs_datums.json"),'w') as f:
                json.dump([{"mol": d.mol.to_adjacency_list(),"value": d.value} for d in pairs_datums],f)
            sampling_datums = []
        else:
            sampling_datums = old_sample_datums + new_datums_E
            with open(os.path.join(path,"Iterations",str(iter),"cumulative_sample_datums.json"),'w') as f:
                json.dump([{"mol": d.mol.to_adjacency_list(),"value": d.value}for d in sampling_datums],f)
            
        Nconfigs = len(admol_name_structure_dict)
        Ncoads = 1
        tree = train_sidt_cov_dep_regressor(pairs_datums,sampling_datums,r_site=r_site,r_morph=r_morph,
                            r_atoms=r_atoms,r_un=r_un,r_lone_pairs=r_lone_pairs,node_fract_training=0.7)
        
        tree_file = os.path.join(path,"Iterations",str(iter),"regressor.json")
        write_nodes(tree,tree_file)
        
        config_E_fws = []
        for coadname in coadnames:
            for admol_name in admol_name_path_dict.keys():
                st = admol_name_structure_dict[admol_name]
                Nocc_isolated = len([a for a in st.atoms if a.is_surface_site() and any(not a2.is_surface_site() for a2 in a.bonds.keys())])
                fw = calculate_configruation_energies_firework(admol_name,tree_file,path,coadname,coad_stable_sites[coadname],Nocc_isolated,
                                                {k.to_adjacency_list(): v for k,v in coadmol_E_dicts[coadname].items()},concern_energy_tol=concern_energy_tol,parents=[],iter=iter,ignore_errors=ignore_errors)
                config_E_fws.append(fw)
        
        ts_frac = self["ts_frac"] if "ts_frac" in self.keys() else None
        sidt_isolated_delta_model_param = self["sidt_isolated_delta_model"] if "sidt_isolated_delta_model" in self.keys() else None
        sidt_covdep_delta_model_param = self["sidt_covdep_delta_model"] if "sidt_covdep_delta_model" in self.keys() else None
        scfw = select_calculations_firework(path,admol_name_path_dict,admol_name_structure_dict,sites,site_adjacency,
                            pynta_dir, metal, facet, slab_path, calculation_directories, coadnames,
                            coad_stable_sites, software, software_kwargs, software_kwargs_TS, freeze_ind, fmaxopt, parents=config_E_fws,Ncalc_per_iter=Ncalc_per_iter,iter=iter,
                            max_iters=max_iters,concern_energy_tol=concern_energy_tol,ignore_errors=ignore_errors,ts_frac=ts_frac,
                            max_coadsorbates=max_coadsorbates,sidt_isolated_delta_model=sidt_isolated_delta_model_param,sidt_covdep_delta_model=sidt_covdep_delta_model_param)
        
        newwf = Workflow(config_E_fws+[scfw],name="Select Calculations "+str(iter))
        
        return FWAction(detours=newwf)

def select_calculations_firework(path,admol_name_path_dict,admol_name_structure_dict,sites,site_adjacency,
                                pynta_dir, metal, facet, slab_path, calculation_directories, coadnames,
                                coad_stable_sites, software, software_kwargs, software_kwargs_TS, freeze_ind, fmaxopt, parents=[],Ncalc_per_iter=6,iter=0,max_iters=6,concern_energy_tol=None,ignore_errors=False,ts_frac=None,max_coadsorbates=None,sidt_isolated_delta_model=None,sidt_covdep_delta_model=None):
    d = {"path": path,"admol_name_path_dict": admol_name_path_dict,"admol_name_structure_dict": {k:v.to_adjacency_list() for k,v in admol_name_structure_dict.items()},
         "sites": sites, "site_adjacency": {str(k): v for k,v in site_adjacency.items()}, "pynta_dir": pynta_dir, "metal": metal, "facet": facet, "slab_path": slab_path,
         "calculation_directories": calculation_directories, "coadnames": coadnames, "coad_stable_sites": coad_stable_sites,
         "Ncalc_per_iter": Ncalc_per_iter, "software": software, "iter": iter, "max_iters": max_iters,
                       "software_kwargs": software_kwargs, "software_kwargs_TS": software_kwargs_TS, "freeze_ind": freeze_ind,
                       "fmaxopt": fmaxopt, "concern_energy_tol": concern_energy_tol, "ignore_errors": ignore_errors, "ts_frac": ts_frac,
                       "max_coadsorbates": max_coadsorbates, "sidt_isolated_delta_model": sidt_isolated_delta_model, "sidt_covdep_delta_model": sidt_covdep_delta_model}
    t1 = SelectCalculationsTask(d)
    return Firework([t1],parents=parents,name="Selecting Calculations "+str(iter),spec={"_priority": 4})

@explicit_serialize
class SelectCalculationsTask(FiretaskBase):
    required_params = ["path","admol_name_path_dict","admol_name_structure_dict","sites","site_adjacency", "pynta_dir", "metal", "facet",
                       "slab_path", "calculation_directories", "coadnames", "coad_stable_sites", "iter", "software", "max_iters",
                       "software_kwargs", "software_kwargs_TS", "freeze_ind", "fmaxopt"]
    optional_params = ["concern_energy_tol","ignore_errors","ts_frac","max_coadsorbates","sidt_isolated_delta_model","sidt_covdep_delta_model"]
    def run_task(self, fw_spec):
        path = self["path"]
        admol_name_path_dict = self["admol_name_path_dict"]
        admol_name_structure_dict = {k: Molecule().from_adjacency_list(v,check_consistency=False) for k,v in self["admol_name_structure_dict"].items()}
        sites = []
        sites = []
        for site in self["sites"]:
            site["normal"] = np.array(site["normal"])
            site["position"] = np.array(site["position"])
            site["indices"] = tuple(site["indices"])
            sites.append(site)
        site_adjacency = {int(k):[int(x) for x in v] for k,v in self["site_adjacency"].items()}
        pynta_dir = self["pynta_dir"]
        metal = self["metal"]
        facet = self["facet"]
        slab_path = self["slab_path"]
        calculation_directories = self["calculation_directories"]
        coadnames  = self["coadnames"]
        Ncalc_per_iter = self["Ncalc_per_iter"]
        iter = self["iter"]
        coad_stable_sites = self["coad_stable_sites"]
        concern_energy_tol = self["concern_energy_tol"] if "concern_energy_tol" in self.keys() else None
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else False
        software_kwargs = self["software_kwargs"]
        software = self["software"]
        freeze_ind = self["freeze_ind"]
        software_kwargs_TS = self["software_kwargs_TS"]
        fmaxopt = self["fmaxopt"]
        max_iters = self["max_iters"]
        
        if iter == max_iters: #terminate
            return FWAction()
        
        coads = {coadname: admol_name_structure_dict[coadname] for coadname in coadnames}
        
        coad_paths = {coadname: os.path.join(pynta_dir,"Adsorbates",coadname) for coadname in coadnames}
        
        slab = read(slab_path)
        nslab = len(slab)
        
        allowed_structure_site_structures = generate_allowed_structure_site_structures(os.path.join(pynta_dir,"Adsorbates"),sites,site_adjacency,nslab,max_dist=np.inf)
        
        ad_energy_dict = get_lowest_adsorbate_energies(os.path.join(pynta_dir,"Adsorbates"))
        coad_Es = {coadname: get_adsorbate_energies(coad_paths[coadname])[0] for coadname in coadnames}
        coadmol_E_dicts = dict()
        coadmol_stability_dicts = dict()
        for coadname in coadnames:
            coadmol_E_dict = dict()
            coadmol_stability_dict = dict()
            for p in os.listdir(coad_paths[coadname]):
                if p == "info.json" or (p not in coad_Es[coadname].keys()):
                    continue
                admol_init,neighbor_sites_init,ninds_init = generate_adsorbate_2D(read(os.path.join(coad_paths[coadname],p,p+"_init.xyz")),sites,site_adjacency,nslab,max_dist=np.inf)
                admol,neighbor_sites,ninds = generate_adsorbate_2D(read(os.path.join(coad_paths[coadname],p,p+".xyz")),sites,site_adjacency,nslab,max_dist=np.inf)
                out_struct = split_adsorbed_structures(admol,clear_site_info=False)[0]
                out_struct_init = split_adsorbed_structures(admol_init,clear_site_info=False)[0]
                coadmol_E_dict[out_struct] = coad_Es[coadname][p] 
                if admol_init.is_isomorphic(admol,save_order=True):
                    coadmol_stability_dict[out_struct_init] = True
                else:
                    coadmol_stability_dict[out_struct_init] = False
            
            coadmol_E_dicts[coadname] = coadmol_E_dict
            coadmol_stability_dicts[coadname] = coadmol_stability_dict

        #cache ad_energy_dict and per-coadname coadmol_E_dict to disk so per-sample
        #ExtractDatumTask workers don't rebuild them (rebuild requires generate_adsorbate_2D
        #with max_dist=np.inf on every adsorbate -- slow)
        ad_energy_dict_path = os.path.join(path, "ad_energy_dict.json")
        with open(ad_energy_dict_path, "w") as f:
            json.dump({m.to_adjacency_list(): v for m,v in ad_energy_dict.items()}, f)
        coadmol_E_dict_paths = {}
        for coadname in coadnames:
            p = os.path.join(path, "coadmol_E_dict_" + coadname + ".json")
            with open(p, "w") as f:
                json.dump({m.to_adjacency_list(): v for m,v in coadmol_E_dicts[coadname].items()}, f)
            coadmol_E_dict_paths[coadname] = p

        #load configurations and Ncoad_energies
        configs_of_concern_by_coad_admol = dict()
        Ncoad_energy_by_coad_admol = dict()
        for coadname in coadnames:
            configs_of_concern_by_coad_admol[coadname] = dict()
            Ncoad_energy_by_coad_admol[coadname] = dict()
            for admol_name,st in admol_name_structure_dict.items():
                config_path = os.path.join(path,"Iterations",str(iter),"configs_of_concern_"+admol_name+"_"+coadname+".json")
                with open(config_path,'r') as f:
                    configs_of_concern_by_coad_admol[coadname][admol_name] = [(Molecule().from_adjacency_list(k[0],check_consistency=False),k[1],k[2],k[3]) for k in json.load(f)]
                Ncoad_energy_path = os.path.join(path,"Iterations",str(iter),"Ncoad_energy_"+admol_name+"_"+coadname+".json")
                with open(Ncoad_energy_path,'r') as f:
                    Ncoad_energy_by_coad_admol[coadname][admol_name] = {int(k):v for k,v in json.load(f).items()}
                
        #load tree
        nodes = read_nodes(os.path.join(path,"Iterations",str(iter),"regressor.json"))
        tree = MultiEvalSubgraphIsomorphicDecisionTreeRegressor([adsorbate_interaction_decomposition],
                                            nodes=nodes)
        
        #load computed configs
        with open(os.path.join(path,"Iterations",str(iter),"computed_configurations.json"),'r') as f:
            computed_configs = [Molecule().from_adjacency_list(x,check_consistency=False) for x in json.load(f)]
        
        
        ts_frac = self["ts_frac"] if "ts_frac" in self.keys() else None
        configs_for_calculation,coad_admol_to_config_for_calculation = get_configs_for_calculation(configs_of_concern_by_coad_admol,Ncoad_energy_by_coad_admol,admol_name_structure_dict,coadnames,computed_configs,tree,Ncalc_per_iter,ts_frac=ts_frac)

        os.makedirs(os.path.join(path,"Iterations",str(iter),"Samples"))
        assert len(configs_for_calculation) > 0, configs_for_calculation
        sample_fws = []
        calculation_directories = []
        for i,config in enumerate(configs_for_calculation):
            adname = None
            coadname = None
            breaking = False
            for cname in coad_admol_to_config_for_calculation.keys():
                for admol_name,config_list in coad_admol_to_config_for_calculation[cname].items():
                    if any(x is config for x in config_list):
                        adname = admol_name
                        coadname = cname
                        breaking = True
                        break
                if breaking:
                    break
            else:
                logging.error(config.to_adjacency_list())
                logging.error(coad_admol_to_config_for_calculation.keys())
                logging.error([[admol_name for admol_name,config_list in coad_admol_to_config_for_calculation[cname].items()] for cname in coad_admol_to_config_for_calculation.keys()])
                raise ValueError
            
            partial_admol = admol_name_structure_dict[adname]
            admol_path = admol_name_path_dict[adname]
            partial_atoms = read(admol_path)
            init_atoms = mol_to_atoms(config,slab,sites,metal,partial_atoms=partial_atoms,partial_admol=partial_admol)
            os.makedirs(os.path.join(path,"Iterations",str(iter),"Samples",str(i)))
            init_path = os.path.join(path,"Iterations",str(iter),"Samples",str(i),"init.xyz")
            write(init_path,init_atoms)
            calculation_directories.append(os.path.split(init_path)[0])
            json_out = {"adjlist": config.to_adjacency_list(), "xyz": admol_path, "coadname": coadname}
            with open(os.path.join(os.path.split(init_path)[0],'info.json'),'w') as f:
                json.dump(json_out,f)
            
            sample_dir = os.path.split(init_path)[0]
            if not any(bd.get_order_str() == 'R' for bd in config.get_all_edges()):
                fwopt = optimize_firework(init_path,
                        software,"weakopt",
                        opt_method="MDMin",opt_kwargs={'dt': 0.05,"trajectory": "weakopt.traj"},software_kwargs=software_kwargs,order=0,
                        run_kwargs={"fmax" : 0.5, "steps" : 30},parents=[],
                            constraints=["freeze up to {}".format(freeze_ind)],
                        ignore_errors=True, metal=metal, facet=facet, priority=3)
                fwopt2 = optimize_firework(os.path.join(sample_dir,"weakopt.xyz"),
                                software,"out",
                                opt_method="QuasiNewton",opt_kwargs={"trajectory": "out.traj"},software_kwargs=software_kwargs,order=0,
                                run_kwargs={"fmax" : fmaxopt, "steps" : 70},parents=[fwopt],
                                constraints=["freeze up to {}".format(freeze_ind)],
                                ignore_errors=True, metal=metal, facet=facet, priority=2)

                fwvib = vibrations_firework(os.path.join(sample_dir,"out.xyz"),
                                            software,"vib",software_kwargs=software_kwargs,parents=[fwopt2],
                                            constraints=["freeze up to "+str(nslab)])
                fwextract = extract_datum_firework(sample_dir, pynta_dir, slab_path, metal, facet, sites, site_adjacency,
                                                   ad_energy_dict_path, coadmol_E_dict_paths[coadname], is_ad=True,
                                                   parents=[fwvib], ignore_errors=True)
                sample_fws.extend([fwopt,fwopt2,fwvib,fwextract])
            else:
                fwopt = optimize_firework(init_path,
                        software,"out", sella=True,
                        opt_kwargs={"trajectory": "out.traj"},software_kwargs=software_kwargs_TS,
                        order=1,
                        run_kwargs={"fmax" : fmaxopt, "steps" : 70},parents=[],
                            constraints=["freeze up to {}".format(freeze_ind)],
                        ignore_errors=True, metal=metal, facet=facet, priority=3)

                fwvib = vibrations_firework(os.path.join(sample_dir,"out.xyz"),
                                            software,"vib",software_kwargs=software_kwargs,parents=[fwopt],
                                            constraints=["freeze up to "+str(nslab)])
                fwextract = extract_datum_firework(sample_dir, pynta_dir, slab_path, metal, facet, sites, site_adjacency,
                                                   ad_energy_dict_path, coadmol_E_dict_paths[coadname], is_ad=False,
                                                   parents=[fwvib], ignore_errors=True)
                sample_fws.extend([fwopt,fwvib,fwextract])
            
        max_coadsorbates_next = self["max_coadsorbates"] if "max_coadsorbates" in self.keys() else None
        sidt_isolated_delta_model_next = self["sidt_isolated_delta_model"] if "sidt_isolated_delta_model" in self.keys() else None
        sidt_covdep_delta_model_next = self["sidt_covdep_delta_model"] if "sidt_covdep_delta_model" in self.keys() else None
        tfw = train_covdep_model_firework(path,admol_name_path_dict,admol_name_structure_dict,sites,site_adjacency,
                            pynta_dir, metal, facet, slab_path, calculation_directories, coadnames,
                            coad_stable_sites, software, software_kwargs, software_kwargs_TS, freeze_ind, fmaxopt, parents=sample_fws,
                            Ncalc_per_iter=Ncalc_per_iter,iter=iter+1,max_iters=max_iters,concern_energy_tol=concern_energy_tol,ignore_errors=ignore_errors,ts_frac=ts_frac,
                            max_coadsorbates=max_coadsorbates_next,sidt_isolated_delta_model=sidt_isolated_delta_model_next,sidt_covdep_delta_model=sidt_covdep_delta_model_next)
        newwf = Workflow(sample_fws+[tfw],name="Train Iteration "+str(iter+1))
        
        return FWAction(detours=newwf)
        

class StructureError(Exception): pass
class TimeLimitError(Exception): pass

@contextmanager
def limit_time(t): #seconds
    def signal_handler(signum, frame):
        raise TimeLimitError
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(int(t))
    try:
        yield
    finally:
        signal.alarm(0)

def index_site_bond_potential_lists_by_site_distances(site_bond_potential_lists):
    if len(site_bond_potential_lists[0]) == 1: #if only one site can't do this analysis
        return np.array(list(range(len(site_bond_potential_lists))))

    max_dists = np.array([get_max_site_dist(L) for L in site_bond_potential_lists])
    return np.argsort(max_dists)

def get_max_site_dist(site_bond_potential):
    vecs = [np.array(d["site_pos"]) for d in site_bond_potential]
    max_dist = 0
    for i in range(len(vecs)):
        for j in range(i):
            if i == j:
                continue
            else:
                v = np.linalg.norm(vecs[i]-vecs[j])
                if v > max_dist:
                    max_dist = v
    return max_dist

def get_task_index(task_dict,task_list):
    """
    get the index of a task in the task list
    """
    for i,d in enumerate(task_list):
        if d['_fw_name'] == task_dict['_fw_name']:
            return i
    else:
        raise IndexError

def restart_opt_firework(task,task_list):
    """
    generate a firework to restart an optimization firework from the trajectory
    file based on the optimization task and the full task list
    """
    traj_file = task["label"]+".traj"
    shutil.copy(traj_file,os.path.join(os.path.split(task["xyz"])[0],traj_file))
    d = deepcopy(task.as_dict())
    d["xyz"] = os.path.join(os.path.split(task["xyz"])[0],traj_file)
    new_task = reconstruct_task(d,orig_task=task)
    return reconstruct_firework(new_task,task,task_list,full=True)

def reconstruct_task(task_dict,orig_task=None):
    """
    regenerate a task based on the task_dict
    """
    name = task_dict["_fw_name"]
    if orig_task:
        fcn = orig_task.__class__
    else:
        fcn = globals()[name]
    d = copy.deepcopy(task_dict)
    return fcn(d)

def reconstruct_firework(new_task,old_task,task_list,full=True):
    """
    reconstruct a firework replacing an old_task with new_task based on the task list
    full indicates whether all of the tasks should be included or just those starting
    from the replaced task
    """
    task_index = get_task_index(old_task.as_dict(),task_list)
    tasks = []
    for i,d in enumerate(task_list):
        if i == task_index:
            tasks.append(new_task)
        elif full or i > task_index:
            tasks.append(reconstruct_task(d))
    return Firework(tasks)

def get_fizzled_fws(lpad):
    return [lpad.get_fw_by_id(i) for i in lpad.get_fw_ids_in_wfs() if lpad.get_fw_by_id(i).state == "FIZZLED"]

def get_completed_fws(lpad):
    return [lpad.get_fw_by_id(i) for i in lpad.get_fw_ids_in_wfs() if lpad.get_fw_by_id(i).state == "COMPLETED"]

def get_fw_traceback_task(fizzfw):
    """
    get the traceback and task errored on for the fizzled firework
    """
    trace = fizzfw.to_dict()["launches"][0]["action"]["stored_data"]["_exception"]["_stacktrace"]
    task_dict = fizzfw.to_dict()["launches"][0]["action"]["stored_data"]["_task"]
    task_list = fizzfw.to_dict()["spec"]["_tasks"]
    task_index = get_task_index(task_dict,task_list)
    task = fizzfw.tasks[task_index]
    return trace,task

def debug_fizzled(fw,trace,task):
    """
    generate a new firework to replace the fizzled firework
    """
    launch_dir = fw.as_dict()["launches"][0]["launch_dir"]
    if issubclass(type(task),OptimizationTask):
        if "ValueError" in trace:
            path = os.path.split(task["xyz"])[0]
            shutil.copy(os.path.join(launch_dir, task["label"]+".traj"),path)
            task_dict = task.as_dict()
            task_dict["xyz"] = os.path.join(path, task["label"]+".traj")
            new_task = reconstruct_task(task_dict,orig_task=task)
            return reconstruct_firework(new_task,task,fw.as_dict()["spec"]["_tasks"])
    else:
        raise ValueError

def restart_wf(lpad,queue):
    """
    lpad is a LaunchPad object
    queue is a boolean indicating if running in a queue or not
    """
    if queue:
        fworker = FWorker.from_file(fireworks.fw_config.FWORKER_LOC)
        qadapter = load_object_from_file(fireworks.fw_config.QUEUEADAPTER_LOC)

    fw_ids = lpad.get_fw_ids()

    wf = lpad.get_wf_by_fw_id(fw_ids[0]) #assumes only one workflow...should fix

    fizzmap = dict()
    fws = [lpad.get_fw_by_id(i) for i in fw_ids]

    for fw in fws:
        fw.parents = []

    for i,fw in enumerate(fws):
        for fw2 in fws:
            if fw2.fw_id in wf.links[fw.fw_id]:
                fw2.parents.append(fw)

    completed_fws = [lpad.get_fw_by_id(idnum) for idnum in fw_ids if lpad.get_fw_by_id(idnum).state == "COMPLETED"]
    complete_inds = [i for i,idnum in enumerate(fw_ids) if lpad.get_fw_by_id(idnum).state == "COMPLETED"]

    fizzfws = [lpad.get_fw_by_id(idnum) for idnum in fw_ids if lpad.get_fw_by_id(idnum).state == "FIZZLED"]

    for fizzfw in fizzfws:
        fw_id = fizzfw.fw_id
        trace,task = get_fw_traceback_task(fizzfw)
        new_fw = debug_fizzled(fizzfw,trace,task)
        fizzmap[fw_ids.index(fw_id)] = new_fw

    for i in fizzmap.keys():
        fw_old = fws[i]
        fw_new = fizzmap[i]
        for fw in fws:
            if fw_old in fw.parents:
                fw.parents.remove(fw_old)
                fw.parents.append(fw_new)
        fws[i] = fw_new

    for ind in sorted(complete_inds)[::-1]:
        fws.pop(ind)

    wf = Workflow(fws, name="restart")
    lpad.reset('', require_password=False)
    lpad.add_wf(wf)

    if queue:
        rapidfirequeue(lpad,fworker,qadapter)
    else:
        rapidfire(lpad)
