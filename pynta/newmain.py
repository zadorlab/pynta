from pynta.tasks import *
from pynta.io import IO
from pynta.adsorbates import Adsorbates
from pynta.molecule import get_adsorbate, generate_unique_site_additions, generate_adsorbate_guesses
from pynta.excatkit.adsorption import Builder
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
from fireworks import LaunchPad, Workflow
from fireworks.queue.queue_launcher import rapidfire as rapidfirequeue
from fireworks.utilities.fw_serializers import load_object_from_file
from fireworks.core.rocket_launcher import rapidfire
from fireworks.core.fworker import FWorker
import fireworks.fw_config

class Pynta:
    def __init__(self,path,launchpad_path,fworker_path,rxns_file,surface_type,metal,label,a=3.6,vaccum=8.0,
        repeats=[(3,3,1),(1,1,4)],slab_path=None,software="Espresso",socket=False,queue=False,njobs_queue=0,
        software_kwargs={'kpts': (3, 3, 1), 'tprnfor': True, 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-11, 'mixing_mode': 'local-TF',
                            "pseudopotentials": {"Cu": 'Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',"H": 'H.pbe-kjpaw_psl.1.0.0.UPF',"O": 'O.pbe-n-kjpaw_psl.1.0.0.UPF',"C": 'C.pbe-n-kjpaw_psl.1.0.0.UPF',"N": 'N.pbe-n-kjpaw_psl.1.0.0.UPF',
                            }, },
        software_kwargs_gas=None,
        reset_launchpad=False,queue_adapter_path=None,nprocs=48,opt_time_limit_hrs=12.0,
        Eharmtol=3.0,Eharmfiltertol=30.0,Ntsmin=5):

        self.surface_type = surface_type
        launchpad = LaunchPad.from_file(launchpad_path)
        if reset_launchpad:
            launchpad.reset('', require_password=False)
        self.launchpad = launchpad
        self.slab_path = slab_path
        self.vaccum = vaccum
        self.a = a
        self.software = software
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

        self.queue = queue
        self.fworker = None
        self.qadapter = None
        self.fworker = FWorker.from_file(fworker_path)
        self.rxns_file = rxns_file
        with open(self.rxns_file,'r') as f:
            self.rxns_dict = yaml.safe_load(f)
        self.slab = read(self.slab_path) if self.slab_path else None
        self.njobs_queue = njobs_queue
        self.label = label
        if queue:
            self.qadapter = load_object_from_file(queue_adapter_path)
        self.nprocs = nprocs
        self.nslab = int(np.prod(np.array(self.repeats[0])*np.array(self.repeats[1])))
        self.mol_dict = None
        self.Eharmtol = Eharmtol
        self.Eharmfiltertol = Eharmfiltertol
        self.Ntsmin = Ntsmin
        self.opt_time_limit_hrs = opt_time_limit_hrs

    def generate_slab(self):
        """
        generates and optimizes a small scale slab that can be scaled to a large slab as needed
        optimization occurs through fireworks and this process waits until the optimization is completed
        """
        slab_type = getattr(ase.build,self.surface_type)
        slab = slab_type(self.metal,self.repeats[1],self.a,self.vaccum)
        slab.pbc = (True, True, False)
        write(os.path.join(self.path,"slab_init.xyz"),slab)
        self.slab_path = os.path.join(self.path,"slab.xyz")
        fwslab = optimize_firework(os.path.join(self.path,"slab_init.xyz"),self.software,"slab",
            opt_method="BFGSLineSearch",socket=self.socket,software_kwargs=self.software_kwargs,
            run_kwargs={"fmax" : 0.01},out_path=os.path.join(self.path,"slab.xyz"))
        wfslab = Workflow([fwslab], name=self.label+"_slab")
        self.launchpad.add_wf(wfslab)
        self.rapidfire()
        while not os.path.exists(self.slab_path): #wait until slab optimizes, this is required anyway and makes the rest of the code simpler
            time.sleep(1)

        self.slab = read(self.slab_path)

    def analyze_slab(self):
        full_slab = self.slab * self.repeats[0]
        cas = SlabAdsorptionSites(full_slab, self.surface_type,allow_6fold=False,composition_effect=False,
                        label_sites=True,
                        surrogate_metal=self.metal)

        self.cas = cas

        single_geoms,single_site_bond_params_lists,single_sites_lists = generate_unique_site_additions(full_slab,cas,len(full_slab),site_bond_params_list=[],sites_list=[])

        double_geoms_full = []
        double_site_bond_params_lists_full = []
        double_sites_lists_full = []
        for i in range(len(single_geoms)):
            double_geoms,double_site_bond_params_lists,double_sites_lists = generate_unique_site_additions(single_geoms[i],
                                                                cas,len(full_slab),single_site_bond_params_lists[i],single_sites_lists[i])
            double_geoms_full.extend(double_geoms)
            double_site_bond_params_lists_full.extend(double_site_bond_params_lists)
            double_sites_lists_full.extend(double_sites_lists)

        self.single_site_bond_params_lists = single_site_bond_params_lists
        self.single_sites_lists = single_sites_lists
        self.double_site_bond_params_lists = double_site_bond_params_lists_full
        self.double_sites_lists = double_sites_lists_full

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

        mol_dict = {mol.to_smiles():mol for mol in unique_mols}
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

    def generate_initial_adsorbate_guesses(self):
        """
        Generates initial guess geometries for adsorbates and gas phase species
        Generates maps connecting the molecule objects with these adsorbates
        """
        grslab = get_grslab(self.slab_path)
        ads_builder = Builder(grslab)
        structures = dict()
        gratom_to_molecule_atom_maps = dict()
        gratom_to_molecule_surface_atom_maps = dict()
        for sm,mol in self.mol_dict.items():
            gratom,surf_indexes,atom_map,surf_index_atom_map = molecule_to_gratoms(mol)
            gratom_to_molecule_atom_maps[sm] = atom_map
            gratom_to_molecule_surface_atom_maps[sm] = surf_index_atom_map
            try:
                if len(surf_indexes) > 0:
                    structs = ads_builder.add_adsorbate(gratom,index=-1,bonds=surf_indexes)
                    structures[sm] = structs
                else:
                    structures[sm] = [gratom]
            except IndexError:
                print('sp_gratoms, sp_gratoms.edges, sp.gratoms.tags')
                print(gratom, gratom.edges,
                                  gratom.get_tags())

        self.gratom_to_molecule_atom_maps = gratom_to_molecule_atom_maps
        self.gratom_to_molecule_surface_atom_maps = gratom_to_molecule_surface_atom_maps
        self.adsorbate_structures = structures

    def setup_adsorbates(self):
        """
        Attaches each adsorbate structure to the slab and sets up fireworks to
        first optimize each possible geometry and then run vibrational calculations on each unique final geometry
        """
        adsorbate_dict = dict()
        for sp_symbol, adsorbate in self.adsorbate_structures.items():
            adsorbate_dict[sp_symbol] = dict()
            for prefix, structure in enumerate(adsorbate):
                adsorbate_dict[sp_symbol][prefix] = structure

        big_slab = self.slab * self.repeats[0]
        nsmall_slab = len(self.slab)
        for adsname,adsorbate in adsorbate_dict.items():
            xyzs = []
            optfws = []
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
                    big_slab_ads = big_slab + structure[nsmall_slab:]
                    software_kwargs = deepcopy(self.software_kwargs)
                    target_site_num = len(mol.get_surface_sites())
                else: #gas phase
                    big_slab_ads = structure
                    target_site_num = None #no slab so can't run site analysis
                    software_kwargs = deepcopy(self.software_kwargs_gas)
                    if len(big_slab_ads) == 1 and self.software == "Espresso": #monoatomic species
                        software_kwargs["command"] = software_kwargs["command"].replace("< PREFIX.pwi > PREFIX.pwo","-ndiag 1 < PREFIX.pwi > PREFIX.pwo")
                try:
                    os.makedirs(os.path.join(self.path,"Adsorbates",adsname,str(prefix)))
                except:
                    pass
                write(os.path.join(self.path,"Adsorbates",adsname,str(prefix),str(prefix)+"_init.xyz"),big_slab_ads)
                xyz = os.path.join(self.path,"Adsorbates",adsname,str(prefix),str(prefix)+".xyz")
                xyzs.append(xyz)
                fwopt = optimize_firework(os.path.join(self.path,"Adsorbates",adsname,str(prefix),str(prefix)+"_init.xyz"),
                    self.software,str(prefix),
                    opt_method="QuasiNewton",socket=self.socket,software_kwargs=software_kwargs,
                    run_kwargs={"fmax" : 0.01, "steps" : 70},parents=[],constraints=["freeze half slab"], time_limit_hrs=self.opt_time_limit_hrs,
                    fmaxhard=0.05, ignore_errors=True, metal=self.metal, facet=self.surface_type, target_site_num=target_site_num)
                optfws.append(fwopt)

            vib_obj_dict = {"software": self.software, "label": adsname, "software_kwargs": software_kwargs,
                "constraints": ["freeze all "+self.metal]}

            cfw = collect_firework(xyzs,True,[["vibrations_firework"]],[[vib_obj_dict]],[["vib.json"]],[[False]],parents=optfws,label=adsname)
            self.adsorbate_fw_dict[adsname] = optfws
            logging.error(self.adsorbate_fw_dict.keys())
            self.fws.extend(optfws+[cfw])


    def setup_transition_states(self,adsorbates_finished=False):
        """
        Sets up fireworks to generate and filter a set of TS estimates, optimize each unique TS estimate,
        and run vibrational and IRC calculations on the each unique final transition state
        Note the vibrational and IRC calculations are launched at the same time
        """
        opt_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                "run_kwargs": {"fmax" : 0.01, "steps" : 70},"constraints": ["freeze half slab"],"sella":True,"order":1,
                "fmaxhard": 0.05, "time_limit_hrs": self.opt_time_limit_hrs}
        vib_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                "constraints": ["freeze all "+self.metal]}
        IRC_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                "run_kwargs": {"fmax" : 0.1, "steps" : 1000},"constraints":["freeze all "+self.metal]}
        for i,rxn in enumerate(self.rxns_dict):
            ts_path = os.path.join(self.path,"TS"+str(i))
            os.makedirs(ts_path)
            ts_task = MolecularTSEstimate({"rxn": rxn,"ts_path": ts_path,"slab_path": self.slab_path,"adsorbates_path": os.path.join(self.path,"Adsorbates"),
                "rxns_file": self.rxns_file,"repeats": self.repeats[0],"path": self.path,"metal": self.metal,"facet": self.surface_type, "out_path": ts_path,
                "spawn_jobs": True, "opt_obj_dict": opt_obj_dict, "vib_obj_dict": vib_obj_dict,
                    "IRC_obj_dict": IRC_obj_dict, "nprocs": self.nprocs, "name_to_adjlist_dict": self.name_to_adjlist_dict,
                    "gratom_to_molecule_atom_maps":{sm: {str(k):v for k,v in d.items()} for sm,d in self.gratom_to_molecule_atom_maps.items()},
                    "gratom_to_molecule_surface_atom_maps":{sm: {str(k):v for k,v in d.items()} for sm,d in self.gratom_to_molecule_surface_atom_maps.items()},
                    "nslab":self.nslab,"Eharmtol":self.Eharmtol,"Eharmfiltertol":self.Eharmfiltertol,"Ntsmin":self.Ntsmin})
            reactants = rxn["reactant_names"]
            products = rxn["product_names"]
            parents = []
            if not adsorbates_finished:
                for m in reactants+products:
                    parents.extend(self.adsorbate_fw_dict[m])
            fw = Firework([ts_task],parents=parents,name="TS"+str(i)+"est")
            self.fws.append(fw)

    def rapidfire(self):
        """
        Call appropriate rapidfire function
        """
        if self.queue:
            rapidfirequeue(self.launchpad,self.fworker,self.qadapter,njobs_queue=self.njobs_queue)
        else:
            rapidfire(self.launchpad,fworker=self.fworker)

    def execute(self):
        if self.slab_path is None: #handle slab
            self.generate_slab()

        self.analyze_slab()
        self.generate_mol_dict()
        self.generate_initial_adsorbate_guesses()

        #adsorbate optimization
        self.setup_adsorbates()

        #setup transition states
        self.setup_transition_states()

        wf = Workflow(self.fws, name=self.label)
        self.launchpad.add_wf(wf)

        while True: #ensures lanuches continue throughout the calculation process
            self.rapidfire()
