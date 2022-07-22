from ase import Atoms
from ase.optimize import *
from ase.constraints import *
from ase.io import write, read
from ase.io.trajectory import Trajectory
from ase.calculators.socketio import SocketIOCalculator
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from ase.vibrations import Vibrations
from sella import Sella, Constraints
from importlib import import_module
from fireworks import *
from fireworks.core.rocket_launcher import rapidfire
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.queue.queue_launcher import rapidfire as rapidfirequeue
from fireworks.utilities.fw_serializers import load_object_from_file
from fireworks.core.fworker import FWorker
import fireworks.fw_config
from pynta.ts import TS
from pynta.vib import AfterTS
from pynta.io import IO
from xtb.ase.calculator import XTB
import json
import copy
import sys
import shutil
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
        software_kwargs = self["software_kwargs"] if "software_kwargs" in self.keys() else dict()
        socket = self["socket"] if "socket" in self.keys() else False
        if socket:
            unixsocket = "ase_"+self["software"].lower()+"_"+self["label"]+"_"+self["xyz"].replace("/","_").replace(".","_")
            if "command" in self["software_kwargs"].keys() and "{unixsocket}" in self["software_kwargs"]["command"]:
                self["software_kwargs"]["command"] = self["software_kwargs"]["command"].format(unixsocket=unixsocket)

        software = name_to_ase_software(self["software"])(**software_kwargs)

        opt_kwargs = self["opt_kwargs"] if "opt_kwargs" in self.keys() else dict()
        opt_method = name_to_ase_opt(self["opt_method"]) if "opt_method" in self.keys() else BFGS
        run_kwargs = self["run_kwargs"] if "run_kwargs" in self.keys() else dict()
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

        if socket and os.path.isfile(os.path.join("/tmp","ipi_"+unixsocket)):
            os.unlink(os.path.join("/tmp","ipi_"+unixsocket))

        sp.calc = SocketIOCalculator(software,log=sys.stdout,unixsocket=unixsocket) if socket else software

        constraints = self["constraints"] if "constraints" in self.keys() else []

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
            except ConnectionResetError as e:
                if socket:
                    sp.calc.close()
                    fw = restart_opt_firework(self,fw_spec["_tasks"])
                    return FWAction(detours=[fw])
                else:
                    if not ignore_errors:
                        raise e
                    else:
                        errors.append(e)
            except Exception as e:
                if not ignore_errors:
                    raise e
                else:
                    errors.append(e)
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
        software_kwargs = self["software_kwargs"] if "software_kwargs" in self.keys() else dict()
        software = name_to_ase_software(self["software"])(**software_kwargs)

        opt_kwargs = self["opt_kwargs"] if "opt_kwargs" in self.keys() else dict()
        opt_method = name_to_ase_opt(self["opt_method"]) if "opt_method" in self.keys() else BFGS
        run_kwargs = self["run_kwargs"] if "run_kwargs" in self.keys() else dict()

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
        software_kwargs = self["software_kwargs"] if "software_kwargs" in self.keys() else dict()
        energy_kwargs = self["energy_kwargs"] if "energy_kwargs" in self.keys() else dict()
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else False

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

def vibrations_firework(xyz,software,label,software_kwargs={},parents=[],out_path=None,constraints=[],ignore_errors=False):
    d = {"xyz" : xyz, "software" : software, "label" : label}
    if software_kwargs: d["software_kwargs"] = software_kwargs
    if constraints: d["constraints"] = constraints
    d["ignore_errors"] = ignore_errors
    t1 = MolecularVibrationsTask(d)
    directory = os.path.dirname(xyz)
    if out_path is None: out_path = os.path.join(directory,label+"_vib.json")
    t2 = FileTransferTask({'files': [{'src': label+'_vib.json', 'dest': os.path.join(os.path.split(out_path)[0],label+'_vib.json')},{'src':'vib.0.traj',
        'dest': out_path}], 'mode': 'copy', 'ignore_errors': ignore_errors})
    return Firework([t1,t2],parents=parents,name=label+"vib")

@explicit_serialize
class MolecularVibrationsTask(VibrationTask):
    required_params = ["software","label"]
    optional_params = ["software_kwargs","constraints","ignore_errors"]
    def run_task(self, fw_spec):
        indices = None
        xyz = self['xyz']
        software = name_to_ase_software(self["software"])
        label = self["label"]
        software_kwargs = self["software_kwargs"] if "software_kwargs" in self.keys() else dict()
        ignore_errors = self["ignore_errors"] if "ignore_errors" in self.keys() else False

        try:
            sp = read(xyz)
            sp.calc = software(**software_kwargs)

            constraints = self["constraints"] if "constraints" in self.keys() else []
            for c in constraints:
                if isinstance(c,dict):
                    constraint = construct_constraint(c)
                    sp.set_constraint(constraint)
                elif c == "freeze half slab":
                    sp.set_constraint(FixAtoms([
                        atom.index for atom in sp if atom.position[2] < sp.cell[2, 2] / 2.
                    ]))
                    indices = [atom.index for atom in atoms if atom.position[2] > atoms.cell[2, 2] / 2.]

            vib = Vibrations(sp,indices=indices)
            vib.run()
            vib.write_mode(n=0)
            vibdata = vib.get_vibrations()
            vibdata.write(label+"_vib.json")
        except Exception as e:
            if not ignore_errors:
                raise e
            else:
                return FWAction(stored_data={"error": e}, exit=True)

        return FWAction()

@explicit_serialize
class MolecularTSEstimate(FiretaskBase):
    required_params = ["rxn","ts_path","slab_path","adsorbates_path","rxns_file","repeats","path","metal"]
    optional_params = ["out_path","scfactor","scfactor_surface","scaled1","scaled2","spawn_jobs","opt_obj_dict",
            "vib_obj_dict","TSnudge_obj_dict"]

    def run_task(self, fw_spec):
        out_path = self["out_path"] if "out_path" in self.keys() else ts_path
        scfactor = self["scfactor"] if "scfactor" in self.keys() else 1.4
        scfactor_surface = self["scfactor"] if "scfactor" in self.keys() else 1.0
        scaled1 = self["scaled1"] if "scaled1" in self.keys() else True
        scaled2 = self["scaled2"] if "scaled2" in self.keys() else True
        spawn_jobs = self["spawn_jobs"] if "spawn_jobs" in self.keys() else False

        ts_path = self["ts_path"]
        rxn = self["rxn"]
        slab = self["slab_path"]
        metal = self["metal"]

        ts = TS(
            self["slab_path"],
            slab,
            ts_path,
            self["rxns_file"],
            self["repeats"],
            self["path"])

        r_name_list, p_name_list = IO.get_reactants_and_products(rxn)
        rxn_name = IO.get_rxn_name(rxn)
        species_dict = IO().get_species_dict(self["rxns_file"])
        rxn_no = rxn["index"]
        species_list = species_dict['rxn' + str(rxn_no)]
        reacting_atoms = IO.get_all_reacting_atoms(self["rxns_file"])['rxn_' + str(rxn_no)]

        ts.TS_placer(
            ts_path,
            scfactor,
            rxn,
            rxn_name,
            r_name_list,
            p_name_list,
            reacting_atoms)

        ts.filtered_out_equiv_ts_estimates(
            ts_path,
            rxn_name)

        #create the xtb-penalty and maybe ts opt fireworks
        # get a dictionary with average distances for all species in
        # species_list, e.g. {'CO2': 4.14.., 'H': 1.665..., 'O': 1.847...}
        sp_surf_av_dists = TS.get_av_dist_dict(
            species_list, metal, self["adsorbates_path"], scfactor_surface,
            scaled1)

        # convert sp_surf_av_dists dict to a tuple and take into accout that
        # the same type of species can be included into penalty function
        # calculations many times, e.g. ['C', 'H', 'O', 'O'], so the av_dist
        # for the 'O' have to be specified twice (order not important)
        av_dists_tuple = TS.get_av_dists_tuple(
            reacting_atoms, sp_surf_av_dists)

        # get all .xyz files with TS estimates
        ts_estimates_xyz_files = []
        ts_est = Path(ts_path).glob('*.xyz')
        for ts in ts_est:
            ts_estimates_xyz_files.append(str(ts))

        # sort it in increasing order
        ts_estimates_xyz_files = sorted(ts_estimates_xyz_files)

        # take a first file and use it as a template to get info about
        # surface atom and adsorbate atoms indices
        tmp_ts_atom = read(ts_estimates_xyz_files[0])
        tmp_adsorbate = tmp_ts_atom[ts.nslab:]

        # convert tag indicies into index indicies
        new_reacting_idx = []
        for reacting_atom_idx in reacting_atoms.values():
            for atom in tmp_adsorbate:
                if atom.tag == reacting_atom_idx:
                    new_reacting_idx.append(atom.index)

        # create surface_atoms_idx dict with all surface atoms and theirs
        # corresponding indicies
        surface_atoms_idxs = {
            atom.symbol + '_' + str(atom.index): atom.index
            for atom in tmp_ts_atom if atom.symbol == metal_atom}

        optfws = []
        xyzs = []
        # Loop through all .xyz files
        for prefix, xyz_file in enumerate(ts_estimates_xyz_files):
            bonds = self.get_bonds_penalty(
                new_reacting_idx,
                surface_atoms_idxs,
                xyz_file)

            # set up some variables
            prefix = str(prefix).zfill(3)
            calc_dir = os.path.join(ts_path, str(prefix))
            os.makedirs(calc_dir, exist_ok=True)

            shutil.copy(xyz_file,calc_dir)
            if spawn_jobs:
                xyz = os.path.join(calc_dir,os.path.basename(xyz_file))
                fwxtb = TSxTBOpt_firework(xyz,self["slab_path"],bonds,self["repeats"],av_dist_tuple,label=str(prefix),parents=[])
                xyzs.append(os.path.join(calc_dir,str(prefix)+".traj"))
                optfws.append(fwxtb)

        if spawn_jobs:
            ["xyzs","check_symm","fw_generators","fw_generator_dicts","out_names","future_check_symms"]
            ctask = MolecularCollect({"xyzs":xyzs,"check_symm":True,"fw_generators": ["optimize_firework","vibrations_firework","TSnudge_firework"],
                "fw_generator_dicts": [self["opt_obj_dict"], self["vib_obj_dict"],
                self["TSnudge_obj_dict"]],
                    "out_names": ["final.xyz","vib.0.traj",""],"future_check_symms": [True,False,False], "label": "TS"+str(rxn_no)+"_"+rxn_name})
            cfw = Firework([ctask],parents=optfws,name="TS"+str(rxn_no)+"_"+rxn_name+"_collect")

        if spawn_jobs:
            return FWAction(additions=optfws,detours=[cfw]) #using detour allows us to inherit children from the original collect to the subsequent collects
        else:
            return FWAction()

def collect_firework(xyzs,check_symm,fw_generators,fw_generator_dicts,out_names,future_check_symms,parents=[],label=""):
    task = MolecularCollect({"xyzs": xyzs, "check_symm": check_symm, "fw_generators": fw_generators,
        "fw_generator_dicts": fw_generator_dicts, "out_names": out_names, "future_check_symms": future_check_symms, "label": label})
    return Firework([task],parents=parents,name=label+"collect")

@explicit_serialize
class MolecularCollect(CollectTask):
    required_params = ["xyzs","check_symm","fw_generators","fw_generator_dicts","out_names","future_check_symms","label"]
    def run_task(self, fw_spec):
        import logging
        logging.error(self["xyzs"])
        xyzs = [xyz for xyz in self["xyzs"] if os.path.isfile(xyz)] #if the associated task errored a file may not be present
        if self["check_symm"]:
            xyzs = get_unique_sym(xyzs) #only unique structures
        if len(xyzs) == 0:
            raise ValueError("No xyzs to collect")
        fw_generator = globals()[self["fw_generators"][0]]
        fws = []
        out_xyzs = []
        for xyz in xyzs:
            d = deepcopy(self["fw_generator_dicts"][0])
            d["xyz"] = xyz
            d["out_path"] = os.path.join(os.path.split(xyz)[0],self["out_names"][0])
            d["label"] = self["out_names"][0]
            d["ignore_errors"] = True
            out_xyzs.append(d["out_path"])
            fw = fw_generator(**d)
            if not isinstance(fw,list):
                fw = [fw]
            fws.extend(fw)

        if len(self["fw_generators"]) > 1:
            task = MolecularCollect({"xyzs": out_xyzs,"check_symm": self["future_check_symms"][0],
                    "fw_generators": self["fw_generators"][1:],"fw_generator_dicts": self["fw_generator_dicts"][1:],
                    "out_names": self["out_names"][1:],"future_check_symms": self["future_check_symms"][1:],"label": self["label"]})
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

def TSxTBOpt_firework(xyz,slab_path,bonds,repeats,av_dist_tuple,out_path=None,label="",parents=[],ignore_errors=False):
    d = {"xyz": xyz, "slab_path": slab_path, "bonds": bonds, "repeats": repeats, "av_dist_tuple": av_dist_tuple,
        "label": label, "ignore_errors": ignore_errors}
    t1 = MolecularTSxTBOpt(d)
    directory = os.path.dirname(xyz)
    if out_path is None: out_path = os.path.join(directory,label+".traj")
    t2 = FileTransferTask({'files': [{'src': label+'.traj', 'dest': out_path}], 'mode': 'copy', "ignore_errors": ignore_errors})
    return Firework([t1,t2],parents=parents,name=label+"TSxTBopt")

class MolecularTSxTBOpt(OptimizationTask):
    required_params = ["xyz","slab_path","bonds","repeats","av_dist_tuple"]
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
        import logging
        logging.error("couldn't find index")
        logging.error(task_dict)
        logging.error(task_list)
        raise IndexError

def restart_opt_firework(task,task_list):
    import logging
    logging.error(task)
    logging.error(task_list)
    traj_file = task["label"]+".traj"
    shutil.copy(traj_file,os.path.join(os.path.split(task["xyz"])[0],traj_file))
    d = deepcopy(task.as_dict())
    d["xyz"] = os.path.join(os.path.split(task["xyz"])[0],traj_file)
    new_task = reconstruct_task(d,orig_task=task)
    return reconstruct_firework(new_task,task,task_list,full=True)

def reconstruct_task(task_dict,orig_task=None):
    import logging
    logging.error("reconstructing task")
    logging.error(task_dict)
    logging.error(orig_task)
    name = task_dict["_fw_name"]
    logging.error(name)
    if orig_task:
        fcn = orig_task.__class__
    else:
        fcn = globals()[name]
    d = copy.deepcopy(task_dict)
    return fcn(d)

def reconstruct_firework(new_task,old_task,task_list,full=True):
    import logging
    logging.error("in reconstruct_firework")
    logging.error(new_task)
    logging.error(old_task)
    logging.error(task_list)
    task_index = get_task_index(old_task.as_dict(),task_list)
    logging.error(task_index)
    tasks = []
    for i,d in enumerate(task_list):
        if i == task_index:
            logging.error("found specified task")
            logging.error(i)
            logging.error(d)
            tasks.append(new_task)
        elif full or i > task_index:
            logging.error("found other task")
            logging.error(i)
            logging.error(d)
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

def get_unique_sym(geoms):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of paths to .xyz or .traj files to compare

    Returns
    ________
    idx_list : list(str)
        a list with prefixes of all symmetrically distinct sites

    '''
    comparator = SymmetryEquivalenceCheck()

    good_adsorbates_atom_obj_list = []
    geos_out = []

    for geom in geoms:
        adsorbate_atom_obj = read(geom)
        adsorbate_atom_obj.pbc = True
        comparision = comparator.compare(
            adsorbate_atom_obj, good_adsorbates_atom_obj_list)

        if comparision is False:
            good_adsorbates_atom_obj_list.append(adsorbate_atom_obj)
            geos_out.append(geom)

    return geos_out
