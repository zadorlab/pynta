from ase import Atoms
from ase.optimize import *
from ase.constraints import *
from ase.io import write, read
from ase.io.trajectory import Trajectory
from ase.calculators.socketio import SocketIOCalculator
from ase.vibrations import Vibrations
from acat.adsorption_sites import SlabAdsorptionSites
from sella import Sella, Constraints, IRC
from importlib import import_module
from fireworks import *
from fireworks.core.rocket_launcher import rapidfire
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.queue.queue_launcher import rapidfire as rapidfirequeue
from fireworks.utilities.fw_serializers import load_object_from_file
from fireworks.core.fworker import FWorker
import fireworks.fw_config
from pynta.ts import TS
from pynta.transitionstate import get_unique_optimized_adsorbates,determine_TS_construction,get_unique_TS_structs,generate_constraints_harmonic_parameters,get_surface_forming_bond_pairings
from pynta.symmetry import get_unique_sym, get_unique_sym_structs
from pynta.vib import AfterTS
from pynta.io import IO
from pynta.penalty_fun import AdsorbatePlacer
from xtb.ase.calculator import XTB
import multiprocessing as mp
import json
import copy
import sys
import shutil
import time
import logging
from copy import deepcopy
from pathlib import Path

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
                      run_kwargs={},constraints=[],parents=[],out_path=None,ignore_errors=False):
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
    d["ignore_errors"] = ignore_errors
    t1 = MolecularOptimizationTask(d)
    directory = os.path.dirname(xyz)
    if out_path is None: out_path = os.path.join(directory,label+".xyz")
    t2 = FileTransferTask({'files': [{'src': label+'.xyz', 'dest': out_path}], 'mode': 'copy', 'ignore_errors' : ignore_errors})
    return Firework([t1,t2],parents=parents,name=label+"opt")

@explicit_serialize
class MolecularOptimizationTask(OptimizationTask):
    required_params = ["software","label"]
    optional_params = ["software_kwargs","opt_method",
        "opt_kwargs","run_kwargs", "constraints","sella","order","socket","ignore_errors"]
    def run_task(self, fw_spec):
        errors = []
        software_kwargs = deepcopy(self["software_kwargs"]) if "software_kwargs" in self.keys() else dict()
        socket = self["socket"] if "socket" in self.keys() else False
        if socket:
            unixsocket = "ase_"+self["software"].lower()+"_"+self["label"]+"_"+self["xyz"].replace("/","_").replace(".","_")
            socket_address = os.path.join("/tmp","ipi_"+unixsocket)
            if "command" in software_kwargs.keys() and "{unixsocket}" in software_kwargs["command"]:
                software_kwargs["command"] = software_kwargs["command"].format(unixsocket=unixsocket)

        software = name_to_ase_software(self["software"])(**software_kwargs)

        opt_kwargs = deepcopy(self["opt_kwargs"]) if "opt_kwargs" in self.keys() else dict()
        opt_method = name_to_ase_opt(self["opt_method"]) if "opt_method" in self.keys() else BFGS
        run_kwargs = deepcopy(self["run_kwargs"]) if "run_kwargs" in self.keys() else dict()
        sella = self["sella"] if "sella" in self.keys() else False
        order = self["order"] if "order" in self.keys() else 0
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else False

        label = self["label"]
        xyz = self['xyz']
        suffix = os.path.split(xyz)[-1].split(".")[-1]

        try:
            if suffix == "xyz":
                sp = read(xyz)
            elif suffix == "traj": #take last point on trajectory
                sp = Trajectory(xyz)[-1]
            else:
                raise ValueError("xyz input not understood")
        except Exception as e:
            if not ignore_errors:
                raise e
            else:
                errors.append(e)

        if socket and os.path.exists(socket_address):
            os.unlink(socket_address)

        sp.calc = SocketIOCalculator(software,log=sys.stdout,unixsocket=unixsocket) if socket else software

        constraints = deepcopy(self["constraints"]) if "constraints" in self.keys() else []

        if not sella:
            assert order == 0
            for c in constraints:
                if isinstance(c,dict):
                    constraint = construct_constraint(c)
                    sp.set_constraint(constraint)
                elif c == "freeze half slab":
                    sp.set_constraint(FixAtoms([
                        atom.index for atom in sp if atom.position[2] < sp.cell[2, 2] / 2.
                    ]))

            opt_kwargs["trajectory"] = label+".traj"

            opt = opt_method(sp,**opt_kwargs)

            try:
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
            cons = Constraints(sp)
            for c in constraints:
                if isinstance(c,dict):
                    add_sella_constraint(cons,c)
                elif c == "freeze half slab":
                    for atom in sp:
                        if atom.position[2] < sp.cell[2, 2] / 2.:
                            cons.fix_translation(atom.index)

            opt = Sella(sp,constraints=cons,trajectory=label+".traj",order=order)
            try:
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

        if not opt.converged():
            e = ValueError
            if not ignore_errors:
                raise e
            else:
                errors.append(e)

        if len(errors) == 0:
            write(label+".xyz", sp)
        else:
            return FWAction(stored_data={"error": errors},exit=True)

        return FWAction()

@explicit_serialize
class MolecularOptimizationFailTask(OptimizationTask):
    required_params = ["software","label"]
    optional_params = ["software_kwargs","opt_method",
        "opt_kwargs","run_kwargs"]
    def run_task(self, fw_spec):
        print(fw_spec)
        software_kwargs = deepcopy(self["software_kwargs"]) if "software_kwargs" in self.keys() else dict()
        software = name_to_ase_software(self["software"])(**software_kwargs)

        opt_kwargs = deepcopy(self["opt_kwargs"]) if "opt_kwargs" in self.keys() else dict()
        opt_method = name_to_ase_opt(self["opt_method"]) if "opt_method" in self.keys() else BFGS
        run_kwargs = deepcopy(self["run_kwargs"]) if "run_kwargs" in self.keys() else dict()

        label = self["label"]
        xyz = self['xyz']
        sp = read(xyz)

        sp.calc = software
        opt = opt_method(sp,trajectory=label+".traj")
        opt.run(fmax=0.02,steps=2)

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
        software = name_to_ase_software(self["software"])
        label = self["label"]
        software_kwargs = deepcopy(self["software_kwargs"]) if "software_kwargs" in self.keys() else dict()
        energy_kwargs = deepcopy(self["energy_kwargs"]) if "energy_kwargs" in self.keys() else dict()
        ignore_errors = deepcopy(self["ignore_errors"]) if "ignore_errors" in self.keys() else False

        try:
            sp = read(xyz)
            sp.calc = software(**software_kwargs)
            en = sp.get_potential_energy(**energy_kwargs)
            with open(label+'_energy.json', 'a') as file:
                wr = json.dump(en,file)
        except Exception as e:
            if not ignore_errors:
                raise e
            else:
                return FWAction(stored_data={"error": e}, exit=True)

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
        software = name_to_ase_software(self["software"])(**software_kwargs)
        ignore_errors = deepcopy(self["ignore_errors"]) if "ignore_errors" in self.keys() else False
        socket = self["socket"] if "socket" in self.keys() else False
        if socket:
            unixsocket = "ase_"+self["software"].lower()+"_"+self["label"]+"_"+self["xyz"].replace("/","_").replace(".","_")
            socket_address = os.path.join("/tmp","ipi_"+unixsocket)
            if "command" in software_kwargs.keys() and "{unixsocket}" in software_kwargs["command"]:
                software_kwargs["command"] = software_kwargs["command"].format(unixsocket=unixsocket)

        if socket and os.path.exists(socket_address):
            os.unlink(socket_address)

        try:
            sp = read(xyz)
            sp.calc = SocketIOCalculator(software,log=sys.stdout,unixsocket=unixsocket) if socket else software

            constraints = deepcopy(self["constraints"]) if "constraints" in self.keys() else []
            for c in constraints:
                if isinstance(c,dict):
                    constraint = construct_constraint(c)
                    sp.set_constraint(constraint)
                elif c == "freeze half slab":
                    sp.set_constraint(FixAtoms([
                        atom.index for atom in sp if atom.position[2] < sp.cell[2, 2] / 2.
                    ]))
                    indices = [atom.index for atom in sp if atom.position[2] > sp.cell[2, 2] / 2.]

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

        return FWAction()


@explicit_serialize
class MolecularTSEstimate(FiretaskBase):
    required_params = ["rxn","ts_path","slab_path","adsorbates_path","rxns_file","repeats","path","metal","facet",
                        "name_to_adjlist_dict", "gratom_to_molecule_atom_maps",
                        "gratom_to_molecule_surface_atom_maps","opt_obj_dict",
                                "vib_obj_dict","IRC_obj_dict","nslab"]
    optional_params = ["out_path","spawn_jobs","nprocs",]
    def run_task(self, fw_spec):
        gratom_to_molecule_atom_maps = {sm: {int(k):v for k,v in d.items()} for sm,d in self["gratom_to_molecule_atom_maps"].items()}
        gratom_to_molecule_surface_atom_maps = {sm: {int(k):v for k,v in d.items()} for sm,d in self["gratom_to_molecule_surface_atom_maps"].items()}
        out_path = self["out_path"] if "out_path" in self.keys() else ts_path
        spawn_jobs = self["spawn_jobs"] if "spawn_jobs" in self.keys() else False
        nprocs = self["nprocs"] if "nprocs" in self.keys() else 1

        ts_path = self["ts_path"]
        rxn = self["rxn"]
        index = rxn["index"]
        metal = self["metal"]
        facet = self["facet"]
        nslab = self["nslab"]

        slab_path = self["slab_path"]
        # slab = read(slab_path)
        # cas = SlabAdsorptionSites(slab,facet,allow_6fold=False,composition_effect=False,
        #                     label_sites=True,
        #                     surrogate_metal=metal)
        #
        # adsorbates_path = self["adsorbates_path"]
        #
        #
        # reactants = Molecule().from_adjacency_list(rxn["reactants"])
        # products = Molecule().from_adjacency_list(rxn["products"])
        #
        # reactant_names = rxn["reactant_names"]
        # product_names = rxn["product_names"]
        #
        # mol_dict = {name: Molecule().from_adjacency_list(adj) for name,adj in self["name_to_adjlist_dict"].items()}
        #
        # reactant_mols = [mol_dict[name] for name in reactant_names]
        # product_mols = [mol_dict[name] for name in product_names]
        #
        # adsorbates = get_unique_optimized_adsorbates(rxn,adsorbates_path)
        #
        # forward,species_names = determine_TS_construction(reactant_names,
        #             reactant_mols,product_names,product_mols)
        #
        # ordered_adsorbates = [adsorbates[name] for name in species_names]
        #
        # rnum_surf_sites = [len(mol.get_surface_sites()) for i,mol in enumerate(reactant_mols)]
        # pnum_surf_sites = [len(mol.get_surface_sites()) for i,mol in enumerate(product_mols)]
        #
        # if forward:
        #     num_surf_sites = rnum_surf_sites
        # else:
        #     num_surf_sites = pnum_surf_sites
        #
        # if forward:
        #     reverse_names = product_names
        # else:
        #     reverse_names = reactant_names
        #
        # tsstructs = get_unique_TS_structs(adsorbates,species_names,cas,nslab,num_surf_sites,mol_dict,
        #                          gratom_to_molecule_atom_maps,gratom_to_molecule_surface_atom_maps,
        #                          facet,metal)
        #
        # constraint_lists,atom_bond_potential_lists,site_bond_potential_lists,site_bond_dict_list,site_fixed_bond_dict_list = generate_constraints_harmonic_parameters(
        #                                     tsstructs,adsorbates,slab,reactants,
        #                                      products,rxn["reaction_family"],template_reversed=(not forward),
        #                                     ordered_names=species_names,reverse_names=reverse_names,
        #                                     mol_dict=mol_dict,gratom_to_molecule_atom_maps=gratom_to_molecule_atom_maps,
        #                                     gratom_to_molecule_surface_atom_maps=gratom_to_molecule_surface_atom_maps,
        #                                     nslab=nslab,facet=facet,metal=metal,cas=cas)
        #
        # out_tsstructs,new_atom_bond_potential_lists,new_site_bond_potential_lists,new_constraint_lists,site_bond_potential_check_lists = get_surface_forming_bond_pairings(
        #                     tsstructs,atom_bond_potential_lists,site_bond_potential_lists,constraint_lists,site_bond_dict_list,
        #                     site_fixed_bond_dict_list,cas)
        #
        #
        # xyzs = []
        # for j,tsstruct in enumerate(out_tsstructs):
        #     os.makedirs(os.path.join(ts_path,str(j)))
        #     write(os.path.join(ts_path,str(j),"xtb_init.xyz"),tsstruct)
        #     sp = run_harmonically_forced_xtb_sella(out_tsstructs[j],new_atom_bond_potential_lists[j],new_site_bond_potential_lists[j],
        #                    nslab=nslab,constraints=new_constraint_lists[j],
        #                    site_bond_potential_check_lists=site_bond_potential_check_lists[j])
        #     if sp:
        #         write(os.path.join(ts_path,str(j),"xtb.xyz"),sp)
        #         xyzs.append(os.path.join(ts_path,str(j),"xtb.xyz"))
        #
        # if spawn_jobs:
        #     ctask = MolecularCollect({"xyzs":xyzs,"check_symm":True,"fw_generators": ["optimize_firework",["vibrations_firework","IRC_firework"]],
        #         "fw_generator_dicts": [self["opt_obj_dict"],[self["vib_obj_dict"],self["IRC_obj_dict"]]],
        #             "out_names": ["opt.xyz",["vib.json","irc.traj"]],"future_check_symms": [True,False], "label": "TS"+str(rxn_no)+"_"+rxn_name})
        #     cfw = Firework([ctask],name="TS"+str(rxn_no)+"_"+rxn_name+"_collect")
        #     newwf = Workflow([cfw],name='rxn_'+str(rxn_no)+str(rxn_name))
        #     return FWAction(detours=newwf) #using detour allows us to inherit children from the original collect to the subsequent collects
        # else:
        #     return FWAction()
        return FWAction()

def collect_firework(xyzs,check_symm,fw_generators,fw_generator_dicts,out_names,future_check_symms,parents=[],label=""):
    task = MolecularCollect({"xyzs": xyzs, "check_symm": check_symm, "fw_generators": fw_generators,
        "fw_generator_dicts": fw_generator_dicts, "out_names": out_names, "future_check_symms": future_check_symms, "label": label})
    return Firework([task],parents=parents,name=label+"collect")

@explicit_serialize
class MolecularCollect(CollectTask):
    required_params = ["xyzs","check_symm","fw_generators","fw_generator_dicts","out_names","future_check_symms","label"]
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
                d["out_path"] = os.path.join(os.path.split(xyz)[0],self["out_names"][0][i])
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
                    "out_names": out_names[1:],"future_check_symms": future_check_symms[1:],"label": self["label"]})
            cfw = Firework([task],parents=fws,name=self["label"]+"collect")
            newwf = Workflow(fws+[cfw],name=self["label"]+"collect"+str(-len(self["fw_generators"])))
            return FWAction(detours=newwf) #using detour allows us to inherit children from the original collect to the subsequent collects
        else:
            return FWAction(detours=fws)

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
        socket=False,software_kwargs={},opt_kwargs={},run_kwargs={},constraints=[],parents=[],ignore_errors=False):
        if out_path is None: out_path = os.path.join(directory,label+"_irc.traj")
        t1 = MolecularIRC(xyz=xyz,label=label,software=software,
            socket=socket,software_kwargs=software_kwargs,opt_kwargs=opt_kwargs,run_kwargs=run_kwargs,
            constraints=constraints,ignore_errors=ignore_errors)
        t2 = FileTransferTask({'files': [{'src': label+'_irc.traj', 'dest': out_path}], 'mode': 'copy', 'ignore_errors' : ignore_errors})
        fw = Firework([t1,t2],parents=[],name=label+"_IRC")
        return fw

@explicit_serialize
class MolecularIRC(FiretaskBase):
    required_params = ["xyz","label"]
    optional_params = ["software","socket",
            "software_kwargs", "opt_kwargs", "run_kwargs", "constraints", "ignore_errors"]
    def run_task(self, fw_spec):
        errors = []
        software_kwargs = deepcopy(self["software_kwargs"]) if "software_kwargs" in self.keys() else dict()
        socket = self["socket"] if "socket" in self.keys() else False
        if socket:
            unixsocket = "ase_"+self["software"].lower()+"_"+self["label"]+"_"+self["xyz"].replace("/","_").replace(".","_")
            socket_address = os.path.join("/tmp","ipi_"+unixsocket)
            if "command" in software_kwargs.keys() and "{unixsocket}" in software_kwargs["command"]:
                software_kwargs["command"] = software_kwargs["command"].format(unixsocket=unixsocket)

        software = name_to_ase_software(self["software"])(**software_kwargs)

        opt_kwargs = deepcopy(self["opt_kwargs"]) if "opt_kwargs" in self.keys() else dict()
        run_kwargs = deepcopy(self["run_kwargs"]) if "run_kwargs" in self.keys() else dict()
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else False

        label = self["label"]
        xyz = self['xyz']
        suffix = os.path.split(xyz)[-1].split(".")[-1]

        try:
            if suffix == "xyz":
                sp = read(xyz)
            elif suffix == "traj": #take last point on trajectory
                sp = Trajectory(xyz)[-1]
            else:
                raise ValueError("xyz input not understood")
        except Exception as e:
            if not ignore_errors:
                raise e
            else:
                errors.append(e)

        if socket and os.path.exists(socket_address):
            os.unlink(socket_address)

        sp.calc = SocketIOCalculator(software,log=sys.stdout,unixsocket=unixsocket) if socket else software

        constraints = deepcopy(self["constraints"]) if "constraints" in self.keys() else []

        cons = Constraints(sp)
        for c in constraints:
            if isinstance(c,dict):
                add_sella_constraint(cons,c)
            elif c == "freeze half slab":
                for atom in sp:
                    if atom.position[2] < sp.cell[2, 2] / 2.:
                        cons.fix_translation(atom.index)

        opt = IRC(sp,constraints=cons,trajectory=label+"_irc.traj",dx=0.1,eta=1e-4,gamma=0.4)
        try:
            run_kwargs["direction"] = "forward"
            opt.run(**run_kwargs)
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

        if not opt.converged():
            e = ValueError
            if not ignore_errors:
                raise e
            else:
                errors.append(e)

        if len(errors) == 0:
            pass
        else:
            return FWAction(stored_data={"error": errors},exit=True)

        return FWAction()

def run_harmonically_forced_xtb(atoms,atom_bond_potentials,site_bond_potentials,nslab,method="GFN1-xTB",
                               site_bond_potential_check_lists=None,constraints=[]):
    cons = Constraints(atoms)

    for c in constraints:
        if isinstance(c,dict):
            add_sella_constraint(cons,c)
        elif c == "freeze slab":
            for i,atom in enumerate(atoms): #freeze the slab
                if i < nslab:
                    cons.fix_translation(atom.index)
        else:
            raise ValueError("Constraint {} not understood".format(c))


    hfxtb = HarmonicallyForcedXTB(method="GFN1-xTB",
                              atom_bond_potentials=atom_bond_potentials,
                             site_bond_potentials=site_bond_potentials)

    atoms.calc = hfxtb

    opt = Sella(atoms,constraints=cons,trajectory="xtbharm.traj",order=0)

    opt.run(fmax=0.02)

    site_bond_correct = True

    # for atm_ind,site_bond_p_list in site_bond_potential_check_lists.items():
    #     used_p = [sbp for sbp in site_bond_potentials if sbp["ind"] == atm_ind][0]
    #     Estar = get_energy_site_bond(atoms,**used_p)
    #     print("Estar {}".format(atm_ind))
    #     print(Estar)
    #     for sbp in site_bond_p_list:
    #         E = get_energy_site_bond(atoms,**sbp)
    #         print(E)
    #         if E-Estar < 0.25:#E < Estar:
    #             print(sbp)
    #             site_bond_correct = False
    #             break
    #
    # print(site_bond_correct)

    if site_bond_correct:
        return atoms
    else:
        return None

def run_parallel_gfn1xtb_opt(inputs,nprocs):
    with mp.Pool(nprocs) as pool:
        errors = pool.map(run_gfn1xtb_opt,inputs)

def run_gfn1xtb_opt(inputs):
    xyz,xyzout,label,slab_path,bonds,av_dists_tuple,repeats = inputs

    adsorbed = read(xyz)
    slab = read(slab_path)
    big_slab = slab * repeats
    nbig_slab = len(big_slab)
    ts_estimate = adsorbed[nbig_slab:]
    traj_path = label+".traj"
    adsplacer = AdsorbatePlacer(
        big_slab, ts_estimate, bonds, av_dists_tuple,
        trajectory=traj_path,
    )

    adsplacer.ads_ref.set_calculator(XTB(method="GFN1-xTB"))

    try:
        opt = adsplacer.optimize()
        write(xyzout,read(traj_path))
        return None
    except Exception as e:
         return e

def TSxTBOpt_firework(xyz,slab_path,bonds,repeats,av_dists_tuple,out_path=None,label="",parents=[],ignore_errors=False):
    d = {"xyz": xyz, "slab_path": slab_path, "bonds": bonds, "repeats": repeats, "av_dists_tuple": av_dists_tuple,
        "label": label, "ignore_errors": ignore_errors}
    t1 = MolecularTSxTBOpt(d)
    directory = os.path.dirname(xyz)
    if out_path is None: out_path = os.path.join(directory,label+".traj")
    t2 = FileTransferTask({'files': [{'src': label+'.traj', 'dest': out_path}], 'mode': 'copy', "ignore_errors": ignore_errors})
    return Firework([t1,t2],parents=parents,name=label+"TSxTBopt")

@explicit_serialize
class MolecularTSxTBOpt(OptimizationTask):
    required_params = ["xyz","slab_path","bonds","repeats","av_dists_tuple"]
    optional_params = ["label","ignore_errors"]

    def run_task(self, fw_spec):
        label = self["label"] if "label" in self.keys() else "xtb"
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else False
        try:
            adsorbed = read(self["xyz"])
            slab = read(self["slab_path"])
            big_slab = slab * self["repeats"]
            nbig_slab = len(big_slab)
            ts_estimate = adsorbed[nbig_slab:]

            adsplacer = AdsorbatePlacer(
                big_slab, ts_estimate, self["bonds"], self["av_dists_tuple"],
                trajectory=label+".traj"
            )

            adsplacer.ads_ref.set_calculator(XTB(method="GFN1-xTB"))
            opt = adsplacer.optimize()
        except Exception as e:
            if not ignore_errors:
                raise e
            else:
                return FWAction(stored_data={"error": e}, exit=True)

        return FWAction()

def name_to_ase_software(software_name):
    if software_name == "XTB":
        module = import_module("xtb.ase.calculator")
        return getattr(module, software_name)
    else:
        module = import_module("ase.calculators."+software_name.lower())
        return getattr(module, software_name)

def name_to_ase_opt(opt_name):
    module = import_module("ase.optimize")
    return getattr(module, opt_name)

def get_task_index(task_dict,task_list):
    for i,d in enumerate(task_list):
        if d['_fw_name'] == task_dict['_fw_name']:
            return i
    else:
        raise IndexError

def restart_opt_firework(task,task_list):
    traj_file = task["label"]+".traj"
    shutil.copy(traj_file,os.path.join(os.path.split(task["xyz"])[0],traj_file))
    d = deepcopy(task.as_dict())
    d["xyz"] = os.path.join(os.path.split(task["xyz"])[0],traj_file)
    new_task = reconstruct_task(d,orig_task=task)
    return reconstruct_firework(new_task,task,task_list,full=True)

def reconstruct_task(task_dict,orig_task=None):
    name = task_dict["_fw_name"]
    if orig_task:
        fcn = orig_task.__class__
    else:
        fcn = globals()[name]
    d = copy.deepcopy(task_dict)
    return fcn(d)

def reconstruct_firework(new_task,old_task,task_list,full=True):
    task_index = get_task_index(old_task.as_dict(),task_list)
    tasks = []
    for i,d in enumerate(task_list):
        if i == task_index:
            tasks.append(new_task)
        elif full or i > task_index:
            tasks.append(reconstruct_task(d))
    return Firework(tasks)

def construct_constraint(d):
    constraint_dict = copy.deepcopy(d)
    constructor = getattr("ase.constraints",constraint_dict["type"])
    del constraint_dict["type"]
    return constructor(**constraint_dict)

def add_sella_constraint(cons,d):
    constraint_dict = copy.deepcopy(d)
    constructor = getattr(cons,constraint_dict["type"])
    del constraint_dict["type"]
    constructor(**constraint_dict)
    return

def get_fizzled_fws(lpad):
    return [lpad.get_fw_by_id(i) for i in lpad.get_fw_ids_in_wfs() if lpad.get_fw_by_id(i).state == "FIZZLED"]

def get_completed_fws(lpad):
    return [lpad.get_fw_by_id(i) for i in lpad.get_fw_ids_in_wfs() if lpad.get_fw_by_id(i).state == "COMPLETED"]

def get_fw_traceback_task(fizzfw):
    trace = fizzfw.to_dict()["launches"][0]["action"]["stored_data"]["_exception"]["_stacktrace"]
    task_dict = fizzfw.to_dict()["launches"][0]["action"]["stored_data"]["_task"]
    task_list = fizzfw.to_dict()["spec"]["_tasks"]
    task_index = get_task_index(task_dict,task_list)
    task = fizzfw.tasks[task_index]
    return trace,task

def debug_fizzled(fw,trace,task):
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
