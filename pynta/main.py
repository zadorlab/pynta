from pynta.tasks import *
from pynta.mol import get_adsorbate, generate_unique_site_additions, generate_adsorbate_guesses, get_name,generate_unique_placements
from molecule.molecule import Molecule
import ase.build
from ase.io import read, write
from ase import Atoms, Atom
from acat.adsorption_sites import SlabAdsorptionSites
from acat.adsorbate_coverage import SlabAdsorbateCoverage
from acat.settings import site_heights, adsorbate_molecule
import os
import time
import yaml
from copy import deepcopy
import numpy as np
from pynta.calculator import get_lattice_parameter
from fireworks import LaunchPad, Workflow
from fireworks.queue.queue_launcher import rapidfire as rapidfirequeue
from fireworks.features.multi_launcher import launch_multiprocess
from fireworks.utilities.fw_serializers import load_object_from_file
from fireworks.core.rocket_launcher import rapidfire
from fireworks.core.fworker import FWorker
import fireworks.fw_config
import logging

class Pynta:
    def __init__(self,path,rxns_file,surface_type,metal,label,launchpad_path=None,fworker_path=None,
        vacuum=8.0,repeats=(3,3,4),slab_path=None,software="Espresso",machine=None,socket=False,queue=False,njobs_queue=0,a=None,
        software_kwargs={'kpts': (3, 3, 1), 'tprnfor': True, 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-6, 'mixing_mode': 'local-TF',
                            "pseudopotentials": {"Cu": 'Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',"H": 'H.pbe-kjpaw_psl.1.0.0.UPF',"O": 'O.pbe-n-kjpaw_psl.1.0.0.UPF',"C": 'C.pbe-n-kjpaw_psl.1.0.0.UPF',"N": 'N.pbe-n-kjpaw_psl.1.0.0.UPF',
                            }, },
        software_kwargs_gas=None,
        TS_opt_software_kwargs=None,
        lattice_opt_software_kwargs={'kpts': (25,25,25), 'ecutwfc': 70, 'degauss':0.02, 'mixing_mode': 'plain'},
        reset_launchpad=False,queue_adapter_path=None,num_jobs=25,max_num_hfsp_opts=None,#max_num_hfsp_opts is mostly for fast testing
        Eharmtol=3.0,Eharmfiltertol=30.0,Ntsmin=5,frozen_layers=2,fmaxopt=0.05,fmaxirc=0.1,fmaxopthard=0.05):

        self.surface_type = surface_type
        if launchpad_path:
            launchpad = LaunchPad.from_file(launchpad_path)
        else:
            launchpad = LaunchPad()

        if reset_launchpad:
            launchpad.reset('', require_password=False)
        self.launchpad = launchpad
        self.slab_path = slab_path
        self.vacuum = vacuum
        self.a = a
        self.software = software
        self.machine = machine #define which machine you are using
        self.socket = socket
        self.repeats = repeats
        self.path = os.getcwd() if path is None else path
        self.facet = metal + surface_type
        self.fws = []
        self.metal = metal
        self.adsorbate_fw_dict = dict()
        self.software_kwargs = software_kwargs

        if software_kwargs_gas:
            self.software_kwargs_gas = software_kwargs_gas
        else:
            self.software_kwargs_gas = deepcopy(software_kwargs)
            self.software_kwargs_gas["kpts"] = 'gamma'
            self.software_kwargs_gas["smearing"] = 'gauss'
            self.software_kwargs_gas["degauss"] = 0.005
            self.software_kwargs_gas["mixing_beta"] = 0.2
            self.software_kwargs_gas["mixing_ndim"] = 10

        self.software_kwargs_TS = deepcopy(software_kwargs)
        if TS_opt_software_kwargs:
            for key,val in TS_opt_software_kwargs.items():
                self.software_kwargs_TS[key] = val

        self.lattice_opt_software_kwargs = deepcopy(software_kwargs)
        if lattice_opt_software_kwargs:
            for key,val in lattice_opt_software_kwargs.items():
                self.lattice_opt_software_kwargs[key] = val

        self.queue = queue
        self.fworker = None
        self.qadapter = None
        if fworker_path:
            self.fworker = FWorker.from_file(fworker_path)
        else:
            self.fworker = FWorker()

        self.rxns_file = rxns_file
        with open(self.rxns_file,'r') as f:
            self.rxns_dict = yaml.safe_load(f)
        self.slab = read(self.slab_path) if self.slab_path else None
        self.njobs_queue = njobs_queue
        self.num_jobs = num_jobs
        self.label = label
        if queue:
            self.qadapter = load_object_from_file(queue_adapter_path)
        if self.slab_path is None:
            self.nslab = int(np.prod(np.array(self.repeats)))
        else:
            self.nslab = len(read(self.slab_path))
        self.layers = self.repeats[2]
        self.frozen_layers = frozen_layers
        self.freeze_ind = int((self.nslab/self.layers)*self.frozen_layers)
        self.mol_dict = None
        self.Eharmtol = Eharmtol
        self.Eharmfiltertol = Eharmfiltertol
        self.Ntsmin = Ntsmin
        self.max_num_hfsp_opts = max_num_hfsp_opts
        self.fmaxopt = fmaxopt
        self.fmaxirc = fmaxirc
        self.fmaxopthard = fmaxopthard

    def generate_slab(self,skip_launch=False):
        """
        generates and optimizes a small scale slab that can be scaled to a large slab as needed
        optimization occurs through fireworks and this process waits until the optimization is completed
        """
        slab_type = getattr(ase.build,self.surface_type)
        #optimize the lattice constant
        if self.a is None:
            a = get_lattice_parameter(self.metal,self.surface_type,self.software,self.lattice_opt_software_kwargs)
            print("computed lattice constant of: {} Angstroms".format(a))
            self.a = a
        else:
            a = self.a
        #construct slab with optimial lattice constant
        slab = slab_type(self.metal,self.repeats,a,self.vacuum)
        slab.pbc = (True, True, False)
        write(os.path.join(self.path,"slab_init.xyz"),slab)
        self.slab_path = os.path.join(self.path,"slab.xyz")
        if self.software != "XTB":
            fwslab = optimize_firework(os.path.join(self.path,"slab_init.xyz"),self.software,"slab",
                opt_method="BFGSLineSearch",socket=self.socket,software_kwargs=self.software_kwargs,
                run_kwargs={"fmax" : self.fmaxopt},out_path=os.path.join(self.path,"slab.xyz"),constraints=["freeze up to {}".format(self.freeze_ind)],priority=1000)
            wfslab = Workflow([fwslab], name=self.label+"_slab")
            self.launchpad.add_wf(wfslab)
            if skip_launch:
                return
            self.launch(single_job=True)
            while not os.path.exists(self.slab_path): #wait until slab optimizes, this is required anyway and makes the rest of the code simpler
                time.sleep(1)
            self.slab = read(self.slab_path)
        else: #testing
            self.slab = slab
            write(self.slab_path,slab)

    def analyze_slab(self):
        full_slab = self.slab
        cas = SlabAdsorptionSites(full_slab, self.surface_type,allow_6fold=False,composition_effect=False,
                        label_sites=True,
                        surrogate_metal=self.metal)

        self.cas = cas

        unique_site_lists,unique_site_pairs_lists,single_site_bond_params_lists,double_site_bond_params_lists = generate_unique_placements(full_slab,cas)

        self.single_site_bond_params_lists = single_site_bond_params_lists
        self.single_sites_lists = unique_site_lists
        self.double_site_bond_params_lists = double_site_bond_params_lists
        self.double_sites_lists = unique_site_pairs_lists

    def generate_mol_dict(self):
        """
        generates all unique Molecule objects based on the reactions and generates a dictionary
        mapping smiles to each unique Molecule object
        also updates self.rxns_dict with addtional useful information for each reaction
        """
        mols = []

        for r in self.rxns_dict:
            r["reactant_mols"] = []
            r["product_mols"] = []
            react = Molecule().from_adjacency_list(r["reactant"])
            prod = Molecule().from_adjacency_list(r["product"])
            react.clear_labeled_atoms()
            prod.clear_labeled_atoms()
            for mol in react.split():
                mol.multiplicity = mol.get_radical_count() + 1
                if not mol.is_surface_site():
                    mols.append(mol)
                    r["reactant_mols"].append(mol)
            for mol in prod.split():
                mol.multiplicity = mol.get_radical_count() + 1
                if not mol.is_surface_site():
                    mols.append(mol)
                    r["product_mols"].append(mol)

        unique_mols = []
        for mol in mols:
            for m in unique_mols:
                if mol.is_isomorphic(m):
                    break
            else:
                unique_mols.append(mol)

        for mol in unique_mols:
            mol.multiplicity = mol.get_radical_count() + 1

        mol_dict = {get_name(mol):mol for mol in unique_mols}
        self.mol_dict = mol_dict
        self.name_to_adjlist_dict = {sm:mol.to_adjacency_list() for sm,mol in mol_dict.items()}


        for r in self.rxns_dict:
            r["reactant_names"] = []
            r["product_names"] = []

            for i,rmol in enumerate(r["reactant_mols"]):
                for sm,mol in mol_dict.items():
                    if mol is rmol or mol.is_isomorphic(rmol,save_order=True):
                        r["reactant_mols"][i] = mol
                        r["reactant_names"].append(sm)
                        break
                else:
                    print("rmol")
                    print(rmol.to_adjacency_list())
                    raise ValueError

            for i,rmol in enumerate(r["product_mols"]):
                for sm,mol in mol_dict.items():
                    if mol is rmol or mol.is_isomorphic(rmol,save_order=True):
                        r["product_mols"][i] = mol
                        r["product_names"].append(sm)
                        break
                else:
                    print("rmol")
                    print(rmol.to_adjacency_list())
                    raise ValueError

            r["reactant_mols"] = [x.to_adjacency_list() for x in r["reactant_mols"]]
            r["product_mols"] = [x.to_adjacency_list() for x in r["product_mols"]]

    def generate_initial_adsorbate_guesses(self,skip_structs=False):
        """
        Generates initial guess geometries for adsorbates and gas phase species
        Generates maps connecting the molecule objects with these adsorbates
        """
        cas = self.cas
        structures = dict()
        gratom_to_molecule_atom_maps = dict()
        gratom_to_molecule_surface_atom_maps = dict()
        for sm,mol in self.mol_dict.items():
            print(sm)
            surf_sites = mol.get_surface_sites()

            ads,mol_to_atoms_map = get_adsorbate(mol)
            if len(surf_sites) == 0:
                if not skip_structs:
                    ads.pbc = (True,True,False)
                    ads.center(vacuum=10)
                    structures[sm] = [ads]
                gratom_to_molecule_atom_maps[sm] = {val:key for key,val in mol_to_atoms_map.items()}
                gratom_to_molecule_surface_atom_maps[sm] = dict()
            else:
                if not skip_structs:
                    #check if this adsorbate is already calculated
                    if os.path.exists(os.path.join(self.path,"Adsorbates",sm)): #assume initial guesses already generated
                        structures[sm] = None
                    else:
                        structs = generate_adsorbate_guesses(mol,ads,self.slab,cas,mol_to_atoms_map,self.metal,
                                           self.single_site_bond_params_lists,self.single_sites_lists,
                                           self.double_site_bond_params_lists,self.double_sites_lists,
                                           self.Eharmtol,self.Eharmfiltertol,self.Ntsmin)
                        structures[sm] = structs


                gratom_to_molecule_atom_maps[sm] = {val:key for key,val in mol_to_atoms_map.items()}

                adatoms = []
                surf_index_atom_map = dict()
                for i,atm in enumerate(mol.atoms):
                    if atm.is_bonded_to_surface():
                        surf_index_atom_map[mol_to_atoms_map[i]] = i

                gratom_to_molecule_surface_atom_maps[sm] = surf_index_atom_map

        self.gratom_to_molecule_atom_maps = gratom_to_molecule_atom_maps
        self.gratom_to_molecule_surface_atom_maps = gratom_to_molecule_surface_atom_maps
        if not skip_structs:
            self.adsorbate_structures = structures

    def setup_adsorbates(self,initial_guess_finished=False):
        """
        Attaches each adsorbate structure to the slab and sets up fireworks to
        first optimize each possible geometry and then run vibrational calculations on each unique final geometry
        """
        if not initial_guess_finished:
            adsorbate_dict = dict()
            for sp_symbol, adsorbate in self.adsorbate_structures.items():
                adsorbate_dict[sp_symbol] = dict()
                if adsorbate is not None:
                    for prefix, structure in enumerate(adsorbate):
                        adsorbate_dict[sp_symbol][prefix] = structure
                else: #read from Adsorbates directory
                    prefixes = os.listdir(os.path.join(self.path,"Adsorbates",sp_symbol))
                    for prefix in prefixes:
                        if prefix == "info.json":
                            continue
                        adsorbate_dict[sp_symbol][prefix] = read(os.path.join(self.path,"Adsorbates",sp_symbol,str(prefix),str(prefix)+"_init.xyz"))

            big_slab = self.slab
            nsmall_slab = len(self.slab)
            for adsname,adsorbate in adsorbate_dict.items():
                xyzs = []
                optfws = []
                optfws2 = []
                mol = self.mol_dict[adsname]

                #check if this adsorbate is already calculated
                exists = False
                for prefix,structure in adsorbate.items():
                    if os.path.exists(os.path.join(self.path,"Adsorbates",adsname,str(prefix),str(prefix)+".xyz")):
                        exists = True

                if exists: #if this species already has a completed opt job in any prefix driectory skip it
                    self.adsorbate_fw_dict[adsname] = []
                    continue

                for prefix,structure in adsorbate.items():
                    if len(mol.get_surface_sites()) > 0:
                        big_slab_ads = structure
                        software_kwargs = deepcopy(self.software_kwargs)
                        target_site_num = len(mol.get_surface_sites())
                        if self.software != "XTB":
                            constraints = ["freeze up to {}".format(self.freeze_ind)]
                        else:
                            constraints = ["freeze up to "+str(self.nslab)]
                    else: #gas phase
                        big_slab_ads = structure
                        target_site_num = None #no slab so can't run site analysis
                        software_kwargs = deepcopy(self.software_kwargs_gas)
                        constraints = []
                        if len(big_slab_ads) == 1 and self.software == "Espresso": #monoatomic species
                            software_kwargs["command"] = software_kwargs["command"].replace("< PREFIX.pwi > PREFIX.pwo","-ndiag 1 < PREFIX.pwi > PREFIX.pwo")
                    try:
                        os.makedirs(os.path.join(self.path,"Adsorbates",adsname,str(prefix)))
                    except:
                        pass
                    write(os.path.join(self.path,"Adsorbates",adsname,str(prefix),str(prefix)+"_init.xyz"),big_slab_ads)
                    sp_dict = {"name":adsname, "adjlist":mol.to_adjacency_list(),"atom_to_molecule_atom_map": self.gratom_to_molecule_atom_maps[adsname],
                            "gratom_to_molecule_surface_atom_map": self.gratom_to_molecule_surface_atom_maps[adsname], "nslab": self.nslab}
                    with open(os.path.join(self.path,"Adsorbates",adsname,"info.json"),'w') as f:
                        json.dump(sp_dict,f)
                    xyz = os.path.join(self.path,"Adsorbates",adsname,str(prefix),str(prefix)+".xyz")
                    xyzs.append(xyz)
                    fwopt = optimize_firework(os.path.join(self.path,"Adsorbates",adsname,str(prefix),str(prefix)+"_init.xyz"),
                        self.software,"weakopt_"+str(prefix),
                        opt_method="MDMin",opt_kwargs={'dt': 0.05},socket=self.socket,software_kwargs=software_kwargs,
                        run_kwargs={"fmax" : 0.5, "steps" : 70},parents=[],constraints=constraints,
                        ignore_errors=True, metal=self.metal, facet=self.surface_type, target_site_num=target_site_num, priority=3)
                    fwopt2 = optimize_firework(os.path.join(self.path,"Adsorbates",adsname,str(prefix),"weakopt_"+str(prefix)+".xyz"),
                        self.software,str(prefix),
                        opt_method="QuasiNewton",socket=self.socket,software_kwargs=software_kwargs,
                        run_kwargs={"fmax" : self.fmaxopt, "steps" : 70},parents=[fwopt],constraints=constraints,
                        ignore_errors=True, metal=self.metal, facet=self.surface_type, target_site_num=target_site_num, priority=3, fmaxhard=self.fmaxopthard,
                        allow_fizzled_parents=True)
                    optfws.append(fwopt)
                    optfws.append(fwopt2)
                    optfws2.append(fwopt2)

                vib_obj_dict = {"software": self.software, "label": adsname, "software_kwargs": software_kwargs,
                    "constraints": ["freeze up to "+str(self.nslab)]}

                cfw = collect_firework(xyzs,True,["vibrations_firework"],[vib_obj_dict],["vib.json"],[],parents=optfws2,label=adsname)
                self.adsorbate_fw_dict[adsname] = optfws2
                logging.error(self.adsorbate_fw_dict.keys())
                self.fws.extend(optfws+[cfw])
        else:
            ads_path = os.path.join(self.path,"Adsorbates")
            for ad in os.listdir(ads_path):
                xyzs = []
                optfws = []
                optfws2 = []
                if ad in self.mol_dict.keys():
                    mol = self.mol_dict[ad]
                else:
                    continue #the species is not in the target reactions so skip it
                target_site_num = len(mol.get_surface_sites())
                ad_path = os.path.join(ads_path,ad)
                completed = False
                for prefix in os.listdir(ad_path):
                    if os.path.exists(os.path.join(ad_path,prefix,prefix+".xyz")):
                        completed = True
                if completed:
                    self.adsorbate_fw_dict[ad] = []
                    continue
                for prefix in os.listdir(ad_path):
                    if prefix.split(".")[-1] == "json":
                        continue
                    prefix_path = os.path.join(ad_path,prefix)
                    if target_site_num == 0:
                        software_kwargs = deepcopy(self.software_kwargs_gas)
                        constraints = []
                        if len(mol.atoms) == 1 and self.software == "Espresso": #monoatomic species
                            software_kwargs["command"] = software_kwargs["command"].replace("< PREFIX.pwi > PREFIX.pwo","-ndiag 1 < PREFIX.pwi > PREFIX.pwo")
                    else:
                        software_kwargs = deepcopy(self.software_kwargs)
                        if self.software != "XTB":
                            constraints = ["freeze up to {}".format(self.freeze_ind)]
                        else:
                            constraints = ["freeze up to "+str(self.nslab)]
                    xyz = os.path.join(prefix_path,str(prefix)+".xyz")
                    init_path = os.path.join(prefix_path,prefix+"_init.xyz")
                    assert os.path.exists(init_path), init_path
                    xyzs.append(xyz)
                    fwopt = optimize_firework(init_path,
                        self.software,"weakopt_"+str(prefix),self.machine,
                        opt_method="MDMin",opt_kwargs={'dt': 0.05},socket=self.socket,software_kwargs=software_kwargs,
                        run_kwargs={"fmax" : 0.5, "steps" : 70},parents=[],constraints=constraints,
                        ignore_errors=True, metal=self.metal, facet=self.surface_type, target_site_num=target_site_num, priority=3)
                    fwopt2 = optimize_firework(os.path.join(self.path,"Adsorbates",ad,str(prefix),"weakopt_"+str(prefix)+".xyz"),
                        self.software,str(prefix),self.machine,
                        opt_method="QuasiNewton",socket=self.socket,software_kwargs=software_kwargs,
                        run_kwargs={"fmax" : self.fmaxopt, "steps" : 70},parents=[fwopt],constraints=constraints,
                        ignore_errors=True, metal=self.metal, facet=self.surface_type, target_site_num=target_site_num, priority=3)
                    optfws.append(fwopt)
                    optfws.append(fwopt2)
                    optfws2.append(fwopt2)

                vib_obj_dict = {"software": self.software, "label": ad, "machine":self.machine, "software_kwargs": software_kwargs,
                    "constraints": ["freeze up to "+str(self.nslab)]}

                cfw = collect_firework(xyzs,True,["vibrations_firework"],[vib_obj_dict],["vib.json"],[True,False],parents=optfws2,label=ad,allow_fizzled_parents=False)
                self.adsorbate_fw_dict[ad] = optfws2
                logging.error(self.adsorbate_fw_dict.keys())
                self.fws.extend(optfws+[cfw])

    def setup_transition_states(self,adsorbates_finished=False):
        """
        Sets up fireworks to generate and filter a set of TS estimates, optimize each unique TS estimate,
        and run vibrational and IRC calculations on the each unique final transition state
        Note the vibrational and IRC calculations are launched at the same time
        """
        if self.software != "XTB":
            opt_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs_TS,
                "run_kwargs": {"fmax" : self.fmaxopt, "steps" : 70},"constraints": ["freeze up to {}".format(self.freeze_ind)],"sella":True,"order":1,}
        else:
            opt_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs_TS,
                "run_kwargs": {"fmax" : 0.02, "steps" : 70},"constraints": ["freeze up to "+str(self.nslab)],"sella":True,"order":1,}
        vib_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                "constraints": ["freeze up to "+str(self.nslab)]}
        IRC_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                "run_kwargs": {"fmax" : self.fmaxirc, "steps" : 70},"constraints":["freeze up to "+str(self.nslab)]}
        for i,rxn in enumerate(self.rxns_dict):
            ts_path = os.path.join(self.path,"TS"+str(i))
            os.makedirs(ts_path)
            ts_task = MolecularTSEstimate({"rxn": rxn,"ts_path": ts_path,"slab_path": self.slab_path,"adsorbates_path": os.path.join(self.path,"Adsorbates"),
                "rxns_file": self.rxns_file,"path": self.path,"metal": self.metal,"facet": self.surface_type, "out_path": ts_path,
                "spawn_jobs": True, "opt_obj_dict": opt_obj_dict, "vib_obj_dict": vib_obj_dict,
                    "IRC_obj_dict": IRC_obj_dict, "nprocs": 48, "name_to_adjlist_dict": self.name_to_adjlist_dict,
                    "gratom_to_molecule_atom_maps":{sm: {str(k):v for k,v in d.items()} for sm,d in self.gratom_to_molecule_atom_maps.items()},
                    "gratom_to_molecule_surface_atom_maps":{sm: {str(k):v for k,v in d.items()} for sm,d in self.gratom_to_molecule_surface_atom_maps.items()},
                    "nslab":self.nslab,"Eharmtol":self.Eharmtol,"Eharmfiltertol":self.Eharmfiltertol,"Ntsmin":self.Ntsmin,
                    "max_num_hfsp_opts":self.max_num_hfsp_opts})
            reactants = rxn["reactant_names"]
            products = rxn["product_names"]
            parents = []
            if not adsorbates_finished:
                for m in reactants+products:
                    parents.extend(self.adsorbate_fw_dict[m])
            fw = Firework([ts_task],parents=parents,name="TS"+str(i)+"est",spec={"_allow_fizzled_parents": True,"_priority": 10})
            self.fws.append(fw)

    def launch(self,single_job=False):
        """
        Call appropriate rapidfire function
        """
        if self.queue:
            rapidfirequeue(self.launchpad,self.fworker,self.qadapter,njobs_queue=self.njobs_queue,nlaunches="infinite")
        elif not self.queue and (self.num_jobs == 1 or single_job):
            rapidfire(self.launchpad,self.fworker,nlaunches="infinite")
        else:
            launch_multiprocess(self.launchpad,self.fworker,"INFO","infinite",self.num_jobs,5)

    def execute(self,generate_initial_ad_guesses=True,calculate_adsorbates=True,
                calculate_transition_states=True,launch=True):
        """
        generate and launch a Pynta Fireworks Workflow
        if generate_initial_ad_guesses is true generates initial guesses, otherwise assumes they are already there
        if calculate_adsorbates is true generates firework jobs for adsorbates, otherwise assumes they are not needed
        if calculate_transition_states is true generates fireworks jobs for transition states, otherwise assumes they are not needed
        if launch is true launches the fireworks workflow in infinite mode...this generates a process that will continue to spawn jobs
        if launch is false the Fireworks workflow is added to the launchpad, where it can be launched separately using fireworks commands
        """

        if not calculate_adsorbates: #default handling
            generate_initial_ad_guesses = False

        if self.slab_path is None: #handle slab
            self.generate_slab()

        self.analyze_slab()
        self.generate_mol_dict()
        self.generate_initial_adsorbate_guesses(skip_structs=(not generate_initial_ad_guesses))

        #adsorbate optimization
        if calculate_adsorbates:
            self.setup_adsorbates(initial_guess_finished=(not generate_initial_ad_guesses))

        if calculate_transition_states:
            #setup transition states
            self.setup_transition_states(adsorbates_finished=(not calculate_adsorbates))

        wf = Workflow(self.fws, name=self.label)
        self.launchpad.add_wf(wf)

        if launch:
            self.launch()


    def execute_from_initial_ad_guesses(self):
        if self.slab_path is None: #handle slab
            self.generate_slab()

        self.analyze_slab()
        self.generate_mol_dict()
        self.generate_initial_adsorbate_guesses(skip_structs=True)

        #adsorbate optimization
        self.setup_adsorbates(initial_guess_finished=True)

        #setup transition states
        self.setup_transition_states()

        wf = Workflow(self.fws, name=self.label)
        self.launchpad.add_wf(wf)


        self.launch()
