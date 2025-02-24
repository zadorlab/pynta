from pynta.tasks import *
from pynta.mol import get_adsorbate, generate_unique_site_additions, generate_adsorbate_guesses, get_name,generate_unique_placements
from pynta.coveragedependence import *
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
from pynta.calculator import get_lattice_parameters
from fireworks import LaunchPad, Workflow
from fireworks.queue.queue_launcher import rapidfire as rapidfirequeue
from fireworks.features.multi_launcher import launch_multiprocess
from fireworks.utilities.fw_serializers import load_object_from_file
from fireworks.core.rocket_launcher import rapidfire
from fireworks.core.fworker import FWorker
import fireworks.fw_config
import logging

#logger
logger = logging.getLogger(__name__)
logging.basicConfig(filename='pynta.log', level=logging.INFO)


class Pynta:
    def __init__(self,path,rxns_file,surface_type,metal,label,launchpad_path=None,fworker_path=None,
        vacuum=8.0,repeats=(3,3,4),slab_path=None,software="Espresso", pbc=(True,True,False),socket=False,queue=False,njobs_queue=0,a=None,
        software_kwargs={'kpts': (3, 3, 1), 'tprnfor': True, 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-6, 'mixing_mode': 'local-TF',
                            "pseudopotentials": {"Cu": 'Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',"H": 'H.pbe-kjpaw_psl.1.0.0.UPF',"O": 'O.pbe-n-kjpaw_psl.1.0.0.UPF',"C": 'C.pbe-n-kjpaw_psl.1.0.0.UPF',"N": 'N.pbe-n-kjpaw_psl.1.0.0.UPF',
                            }, },
        software_kwargs_gas=None,
        TS_opt_software_kwargs=None,
        irc_mode="fixed", #choose irc mode: 'skip', 'relaxed', 'fixed'
        lattice_opt_software_kwargs={'kpts': (25,25,25), 'ecutwfc': 70, 'degauss':0.02, 'mixing_mode': 'plain'},
        reset_launchpad=False,queue_adapter_path=None,num_jobs=25,max_num_hfsp_opts=None,#max_num_hfsp_opts is mostly for fast testing
        Eharmtol=3.0,Eharmfiltertol=30.0,Ntsmin=5,frozen_layers=2,fmaxopt=0.05,fmaxirc=0.1,fmaxopthard=0.05,c=None,
        surrogate_metal=None):

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
        self.pbc = pbc
        self.c = c
        self.software = software
        self.socket = socket
        self.repeats = repeats
        self.path = os.getcwd() if path is None else path
        self.facet = metal + surface_type
        self.fws = []
        self.metal = metal
        if surrogate_metal is None:
            self.surrogate_metal = metal
        else:
            self.surrogate_metal = surrogate_metal
        self.adsorbate_fw_dict = dict()
        self.software_kwargs = software_kwargs
        self.irc_mode = irc_mode

        if software.lower() == 'vasp':
            self.pbc = (True,True,True)

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

        logger.info('Pynta class is initiated')

    def generate_slab(self,skip_launch=False):
        """
        generates and optimizes a small scale slab that can be scaled to a large slab as needed
        optimization occurs through fireworks and this process waits until the optimization is completed
        """
        slab_type = getattr(ase.build,self.surface_type)
        #optimize the lattice constant
        if self.a is None:
            a = get_lattice_parameters(self.metal,self.surface_type,self.software,self.lattice_opt_software_kwargs)
            print("computed lattice constants of: {} Angstroms".format(a))
            if isinstance(a,float):
                self.a = a
            else:
                self.a = a[0]
                self.c = a[1]
        
        logger.info('Construct slab with optimal lattice constant')
        #construct slab with optimial lattice constant
        if self.c:
            slab = slab_type(symbol=self.metal,size=self.repeats,a=self.a,vacuum=self.vacuum,c=self.c)
        else:
            slab = slab_type(symbol=self.metal,size=self.repeats,a=self.a,vacuum=self.vacuum)
        slab.pbc = self.pbc
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
                        surrogate_metal=self.surrogate_metal)

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
                    ads.pbc = self.pbc
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
                        vib_constraints = ["freeze up to "+str(self.nslab)]
                    else: #gas phase
                        big_slab_ads = structure
                        target_site_num = None #no slab so can't run site analysis
                        software_kwargs = deepcopy(self.software_kwargs_gas)
                        constraints = []
                        vib_constraints = []
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
                    "constraints": vib_constraints}

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
                        vib_constraints = []
                        if len(mol.atoms) == 1 and self.software == "Espresso": #monoatomic species
                            software_kwargs["command"] = software_kwargs["command"].replace("< PREFIX.pwi > PREFIX.pwo","-ndiag 1 < PREFIX.pwi > PREFIX.pwo")
                    else:
                        software_kwargs = deepcopy(self.software_kwargs)
                        if self.software != "XTB":
                            constraints = ["freeze up to {}".format(self.freeze_ind)]
                        else:
                            constraints = ["freeze up to "+str(self.nslab)]
                        vib_constraints = ["freeze up to "+str(self.nslab)]
                    xyz = os.path.join(prefix_path,str(prefix)+".xyz")
                    init_path = os.path.join(prefix_path,prefix+"_init.xyz")
                    assert os.path.exists(init_path), init_path
                    xyzs.append(xyz)
                    fwopt = optimize_firework(init_path,
                        self.software,"weakopt_"+str(prefix),
                        opt_method="MDMin",opt_kwargs={'dt': 0.05},socket=self.socket,software_kwargs=software_kwargs,
                        run_kwargs={"fmax" : 0.5, "steps" : 70},parents=[],constraints=constraints,
                        ignore_errors=True, metal=self.metal, facet=self.surface_type, target_site_num=target_site_num, priority=3)
                    fwopt2 = optimize_firework(os.path.join(self.path,"Adsorbates",ad,str(prefix),"weakopt_"+str(prefix)+".xyz"),
                        self.software,str(prefix),
                        opt_method="QuasiNewton",socket=self.socket,software_kwargs=software_kwargs,
                        run_kwargs={"fmax" : self.fmaxopt, "steps" : 70},parents=[fwopt],constraints=constraints,
                        ignore_errors=True, metal=self.metal, facet=self.surface_type, target_site_num=target_site_num, priority=3)
                    optfws.append(fwopt)
                    optfws.append(fwopt2)
                    optfws2.append(fwopt2)

                vib_obj_dict = {"software": self.software, "label": ad, "software_kwargs": software_kwargs,
                    "constraints": vib_constraints}

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

        #logging.info
        logger.info(f"================= IRC mode is: {self.irc_mode} =======================")
        #pass through 
        
        for i,rxn in enumerate(self.rxns_dict):
            #if irc_mode is "fixed" freeze all slab and conduct MolecularTSEstimate. 
            if self.irc_mode == "fixed":
                IRC_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                    "run_kwargs": {"fmax" : self.fmaxopt, "steps" : 70},"constraints": ["freeze up to "+str(self.nslab)]}

            elif self.irc_mode == "relaxed":
                IRC_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                    "run_kwargs": {"fmax" : self.fmaxopt, "steps" : 70},"constraints": ["freeze up to {}".format(self.freeze_ind)]}
        # if irc_mode = "skip" : do not conduct IRC
            else:
                logger.info("Skip IRC: IRC is not conducted")
                IRC_obj_dict = {}
                pass

            ts_path = os.path.join(self.path,"TS"+str(i))
            os.makedirs(ts_path)
            ts_task = MolecularTSEstimate({"rxn": rxn,"ts_path": ts_path,"slab_path": self.slab_path,"adsorbates_path": os.path.join(self.path,"Adsorbates"),
                    "rxns_file": self.rxns_file,"path": self.path,"metal": self.metal,"facet": self.surface_type, "out_path": ts_path, "irc_mode": self.irc_mode,
                    "spawn_jobs": True, "opt_obj_dict": opt_obj_dict, "vib_obj_dict": vib_obj_dict, "IRC_obj_dict": IRC_obj_dict,
                    "nprocs": 48, "name_to_adjlist_dict": self.name_to_adjlist_dict,
                    "gratom_to_molecule_atom_maps":{sm: {str(k):v for k,v in d.items()} for sm,d in self.gratom_to_molecule_atom_maps.items()},
                    "gratom_to_molecule_surface_atom_maps":{sm: {str(k):v for k,v in d.items()} for sm,d in self.gratom_to_molecule_surface_atom_maps.items()},
                    "nslab":self.nslab,"Eharmtol":self.Eharmtol,"Eharmfiltertol":self.Eharmfiltertol,"Ntsmin":self.Ntsmin,
                    "max_num_hfsp_opts":self.max_num_hfsp_opts, "surrogate_metal":self.surrogate_metal})
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

class CoverageDependence:
    def __init__(self,path,metal,surface_type,repeats,pynta_run_directory,software,software_kwargs,label,sites,site_adjacency,coad_stable_sites,adsorbates=[],transition_states=dict(),coadsorbates=[],
                 max_dist=3.0,frozen_layers=2,fmaxopt=0.05,Ncalc_per_iter=6,TS_opt_software_kwargs=None,launchpad_path=None,
                 fworker_path=None,queue=False,njobs_queue=0,reset_launchpad=False,queue_adapter_path=None,
                 num_jobs=25,surrogate_metal=None,concern_energy_tol=None,max_iters=np.inf):
        self.path = path
        self.metal = metal
        self.repeats = repeats
        self.nslab = np.product(repeats)
        self.pynta_run_directory = pynta_run_directory
        self.pairs_directory = os.path.join(self.path,"pairs")
        self.slab_path = os.path.join(self.pynta_run_directory,"slab.xyz")
        self.adsorbates_path = os.path.join(self.pynta_run_directory,"Adsorbates")
        self.coad_stable_sites = coad_stable_sites
        self.software = software
        self.software_kwargs = software_kwargs
        self.surface_type = surface_type
        self.facet = metal + surface_type
        self.sites = sites
        self.site_adjacency = site_adjacency
        self.Ncalc_per_iter = Ncalc_per_iter
        self.software_kwargs_TS = deepcopy(software_kwargs)
        self.concern_energy_tol = concern_energy_tol
        if TS_opt_software_kwargs:
            for key,val in TS_opt_software_kwargs.items():
                self.software_kwargs_TS[key] = val
        
        if adsorbates != []:
            raise ValueError("Not implemented yet")
        self.adsorbates = adsorbates
        self.transition_states = transition_states
        self.coadsorbates = coadsorbates
        self.max_dist = max_dist
        self.frozen_layers = frozen_layers
        self.layers = self.repeats[2]
        self.freeze_ind = int((self.nslab/self.layers)*self.frozen_layers)
        self.fmaxopt = fmaxopt
        self.label = label
        self.max_iters = max_iters
        
        if launchpad_path:
            launchpad = LaunchPad.from_file(launchpad_path)
        else:
            launchpad = LaunchPad()

        if reset_launchpad:
            launchpad.reset('', require_password=False)
        self.launchpad = launchpad
        
        self.fworker_path = fworker_path
        if fworker_path:
            self.fworker = FWorker.from_file(fworker_path)
        else:
            self.fworker = FWorker()
            
        self.queue = queue
        self.njobs_queue = njobs_queue
        self.reset_launchpad = reset_launchpad
        if queue:
            self.qadapter = load_object_from_file(queue_adapter_path)
        self.num_jobs = num_jobs
        
        if surrogate_metal is None:
            self.surrogate_metal = metal
        else:
            self.surrogate_metal = surrogate_metal 
            
        self.pairs_fws = []
        self.fws = []
        
    def setup_pairs_calculations(self):
        tsdirs = [os.path.join(self.pynta_run_directory,t,ind) for t,ind in self.transition_states.items()]
        outdirs_ad,outdirs_ts = setup_pair_opts_for_rxns(self.path,tsdirs,self.coadsorbates,self.surrogate_metal,self.surface_type,max_dist=self.max_dist)
        
        for d in outdirs_ad:
            fwopt = optimize_firework(d,
                            self.software,"weakopt",
                            opt_method="MDMin",opt_kwargs={'dt': 0.05,"trajectory": "weakopt.traj"},software_kwargs=self.software_kwargs,order=0,
                            run_kwargs={"fmax" : 0.5, "steps" : 30},parents=[],
                              constraints=["freeze up to {}".format(self.freeze_ind)],
                            ignore_errors=True, metal=self.metal, facet=self.surface_type, priority=3)
            fwopt2 = optimize_firework(os.path.join(os.path.split(d)[0],"weakopt.xyz"),
                            self.software,"out",
                            opt_method="QuasiNewton",opt_kwargs={"trajectory": "out.traj"},software_kwargs=self.software_kwargs,order=0,
                            run_kwargs={"fmax" : self.fmaxopt, "steps" : 70},parents=[fwopt],
                              constraints=["freeze up to {}".format(self.freeze_ind)],
                            ignore_errors=True, metal=self.metal, facet=self.surface_type, priority=2)
        
            fwvib = vibrations_firework(os.path.join(os.path.split(d)[0],"out.xyz"),
                                        self.software,"vib",software_kwargs=self.software_kwargs,parents=[fwopt2],
                                        constraints=["freeze up to "+str(self.nslab)])
            self.pairs_fws.append(fwopt)
            self.pairs_fws.append(fwopt2)
            self.pairs_fws.append(fwvib)
        
        for d in outdirs_ts:
            fwopt = optimize_firework(d,
                            self.software,"out", sella=True, 
                            opt_kwargs={"trajectory": "out.traj"},software_kwargs=self.software_kwargs_TS,
                            order=1,
                            run_kwargs={"fmax" : self.fmaxopt, "steps" : 70},parents=[],
                              constraints=["freeze up to {}".format(self.freeze_ind)],
                            ignore_errors=True, metal=self.metal, facet=self.surface_type, priority=3)
            
            fwvib = vibrations_firework(os.path.join(os.path.split(d)[0],"out.xyz"),
                                        self.software,"vib",software_kwargs=self.software_kwargs,parents=[fwopt],
                                        constraints=["freeze up to "+str(self.nslab)])
            self.pairs_fws.append(fwopt)
            self.pairs_fws.append(fwvib)
            
        self.fws.extend(self.pairs_fws)
    
    def setup_active_learning_loop(self):
        admol_name_path_dict = {k: os.path.join(self.pynta_run_directory,k,v,"opt.xyz") for k,v in self.transition_states.items()}
        admol_name_structure_dict = dict()
        ads = self.adsorbates
        allowed_structure_site_structures = generate_allowed_structure_site_structures(os.path.join(self.pynta_run_directory,"Adsorbates"),self.sites,self.site_adjacency,self.nslab,max_dist=np.inf)

        for ts in self.transition_states.keys():
            info_path = os.path.join(self.pynta_run_directory,ts,"info.json")
            with open(info_path,'r') as f:
                info = json.load(f)
            reactants = Molecule().from_adjacency_list(info["reactants"])
            products = Molecule().from_adjacency_list(info["products"])
            keep_binding_vdW_bonds_in_reactants=False
            keep_vdW_surface_bonds_in_reactants=False
            mol = reactants
            for bd in mol.get_all_edges():
                if bd.order == 0:
                    if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
                        keep_binding_vdW_bonds_in_reactants = True
                        m = mol.copy(deep=True)
                        b = m.get_bond(m.atoms[mol.atoms.index(bd.atom1)],m.atoms[mol.atoms.index(bd.atom2)])
                        m.remove_bond(b)
                        out = m.split()
                        if len(out) == 1: #vdW bond is not only thing connecting adsorbate to surface
                            keep_vdW_surface_bonds_in_reactants = True
            keep_binding_vdW_bonds_in_products=False
            keep_vdW_surface_bonds_in_products=False
            mol = products
            for bd in mol.get_all_edges():
                if bd.order == 0:
                    if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
                        keep_binding_vdW_bonds_in_products = True
                        m = mol.copy(deep=True)
                        b = m.get_bond(m.atoms[mol.atoms.index(bd.atom1)],m.atoms[mol.atoms.index(bd.atom2)])
                        m.remove_bond(b)
                        out = m.split()
                        if len(out) == 1: #vdW bond is not only thing connecting adsorbate to surface
                            keep_vdW_surface_bonds_in_products = True
            
            keep_binding_vdW_bonds = keep_binding_vdW_bonds_in_reactants and keep_binding_vdW_bonds_in_products
            keep_vdW_surface_bonds = keep_vdW_surface_bonds_in_reactants and keep_vdW_surface_bonds_in_products
            
            atoms = read(admol_name_path_dict[ts])
            st,_,_ = generate_TS_2D(atoms, info_path,  self.metal, self.surface_type, self.sites, self.site_adjacency, self.nslab,
                     max_dist=np.inf, allowed_structure_site_structures=allowed_structure_site_structures,
                     keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
            admol_name_structure_dict[ts] = st
            with open(info_path,"r") as f:
                info = json.load(f)
                for name in info["species_names"]+info["reverse_names"]:
                    if name not in ads:
                        ads.append(name)
        
        for ad in ads:
            p = os.path.join(self.pynta_run_directory,"Adsorbates",ad)
            with open(os.path.join(p,"info.json")) as f:
                info = json.load(f)
                
            mol = Molecule().from_adjacency_list(info["adjlist"])
            
            if mol.contains_surface_site():
                keep_binding_vdW_bonds=False 
                keep_vdW_surface_bonds=False
                for bd in mol.get_all_edges():
                    if bd.order == 0:
                        if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
                            keep_binding_vdW_bonds = True
                            m = mol.copy(deep=True)
                            b = m.get_bond(m.atoms[mol.atoms.index(bd.atom1)],m.atoms[mol.atoms.index(bd.atom2)])
                            m.remove_bond(b)
                            out = m.split()
                            if len(out) == 1: #vdW bond is not only thing connecting adsorbate to surface
                                keep_vdW_surface_bonds = True
                        
                ad_xyz = get_best_adsorbate_xyz(p,self.sites,self.nslab)
                admol_name_path_dict[ad] = ad_xyz 
                atoms = read(ad_xyz)
                st,_,_ = generate_adsorbate_2D(atoms, self.sites, self.site_adjacency, self.nslab, max_dist=np.inf, allowed_structure_site_structures=allowed_structure_site_structures,
                                               keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
                admol_name_structure_dict[ad] = st
                
        calculation_directories = [] #identify pairs directories
        for p1 in os.listdir(os.path.join(self.path,"pairs")):
            for p2 in os.listdir(os.path.join(self.path,"pairs",p1)):
                calculation_directories.append(os.path.join(self.path,"pairs",p1,p2))
        
        
        fw = train_covdep_model_firework(self.path,admol_name_path_dict,admol_name_structure_dict,self.sites,self.site_adjacency,
                                self.pynta_run_directory, self.metal, self.surface_type, self.slab_path, calculation_directories, self.coadsorbates[0], 
                                self.coad_stable_sites, self.software, self.software_kwargs, self.software_kwargs_TS, self.freeze_ind, self.fmaxopt,
                                parents=self.fws, max_iters=self.max_iters,
                                Ncalc_per_iter=self.Ncalc_per_iter,iter=0,concern_energy_tol=self.concern_energy_tol,ignore_errors=True)

        self.fws.append(fw)
        
        
    def execute(self,run_pairs=True,run_active_learning=False):
        if run_pairs:
            self.setup_pairs_calculations()
        if run_active_learning:
            self.setup_active_learning_loop()
        wf = Workflow(self.fws, name=self.label)
        self.launchpad.add_wf(wf)