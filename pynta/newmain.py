from pynta.tasks import *
from pynta.io import IO
from pynta.adsorbates import Adsorbates
from pynta.molecule import *
import ase.build
from ase.io import read, write
import os
import time
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
                            }},
        reset_launchpad=False,queue_adapter_path=None,nprocs=48):

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
        self.rxns = IO.open_yaml_file(rxns_file)
        self.metal = metal
        self.adsorbate_fw_dict = dict()
        self.software_kwargs = software_kwargs
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
        self.nslab = np.prod(np.array(self.repeats[0])*np.array(self.repeats[1]))
        self.mol_dict = None


    def generate_slab(self):
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

    def generate_mol_dict(self):
        mols = []
        for r in self.rxn_dict:
            react = Molecule().from_adjacency_list(r["reactant"])
            prod = Molecule().from_adjacency_list(r["product"])
            react.clear_labeled_atoms()
            prod.clear_labeled_atoms()
            for mol in react.split()+prod.split():
                if mol.contains_surface_site() and not mol.is_surface_site():
                    mols.append(mol)

        unique_mols = []
        for mol in mols:
            for m in unique_mols:
                if mol.is_isomorphic(m):#,save_order=True):
                    break
            else:
                unique_mols.append(mol)

        mol_dict = {mol.to_smiles():mol for mol in unique_mols}
        self.mol_dict = mol_dict

    def generate_initial_adsorbate_guesses(self):
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
                structs = ads_builder.add_adsorbate(gratom,index=-1,bonds=surf_indexes)
                structures[sm] = structs
            except IndexError:
                print('sp_gratoms, sp_gratoms.edges, sp.gratoms.tags')
                print(gratom, gratom.edges,
                                  gratom.get_tags())

        self.gratom_to_molecule_atom_maps = gratom_to_molecule_atom_maps
        self.gratom_to_molecule_surface_atom_maps = gratom_to_molecule_surface_atom_maps
        self.adsorbate_structures = structures

    def setup_adsorbates(self):
        nslab = self.nslab

        adsorbate_dict = dict()
        for sp_symbol, adsorbate in self.adsorbate_structures.items():
            adsorbate_dict[sp_symbol] = dict()
            for prefix, structure in enumerate(adsorbate):
                adsorbate_dict[sp_symbol][prefix] = structure
        
        big_slab = self.slab * self.repeats[0]
        for adsname,adsorbate in adsorbate_dict.items():
            xyzs = []
            optfws = []
            for prefix,structure in adsorbate.items():
                big_slab_ads = big_slab + structure[nslab:]
                os.makedirs(os.path.join(self.path,"Adsorbates",adsname,str(prefix)))
                write(os.path.join(self.path,"Adsorbates",adsname,str(prefix),str(prefix)+"_init.xyz"),big_slab_ads)
                xyz = os.path.join(self.path,"Adsorbates",adsname,str(prefix),str(prefix)+".xyz")
                xyzs.append(xyz)
                fwopt = optimize_firework(os.path.join(self.path,"Adsorbates",adsname,str(prefix),str(prefix)+"_init.xyz"),
                    self.software,str(prefix),
                    opt_method="QuasiNewton",socket=self.socket,software_kwargs=self.software_kwargs,
                    run_kwargs={"fmax" : 0.01, "steps" : 70},parents=[],constraints=["freeze half slab"], ignore_errors=True)
                optfws.append(fwopt)

            vib_obj_dict = {"software": self.software, "label": str(prefix), "software_kwargs": self.software_kwargs,
                "constraints": ["freeze half slab"]}

            cfw = collect_firework(xyzs,False,[["vibrations_firework"]],[[vib_obj_dict]],[["vib.json"]],[[False]],parents=optfws,label=adsname)
            self.adsorbate_fw_dict[adsname] = optfws
            logging.error(self.adsorbate_fw_dict.keys())
            self.fws.extend(optfws+[cfw])


    def setup_transition_states(self,adsorbates_finished=False):
        opt_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                "run_kwargs": {"fmax" : 0.01, "steps" : 70},"constraints": ["freeze half slab"],"sella":True,"order":1}
        vib_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                "constraints": ["freeze half slab"]}
        IRC_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                "constraints":["freeze half slab"]}
        for i,rxn in enumerate(self.rxns):
            ts_path = os.path.join(self.path,"TS"+str(i))
            os.makedirs(ts_path)
            ts_task = MolecularTSEstimate({"rxn": rxn,"ts_path": ts_path,"slab_path": self.slab_path,"adsorbates_path": os.path.join(self.path,"Adsorbates"),
                "rxns_file": self.rxns_file,"repeats": self.repeats[0],"path": self.path,"metal": self.metal,"out_path": ts_path,
                    "scfactor": 1.4,"scfactor_surface": 1.0,
                    "scaled1": True, "scaled2": False, "spawn_jobs": True, "opt_obj_dict": opt_obj_dict, "vib_obj_dict": vib_obj_dict,
                    "IRC_obj_dict": IRC_obj_dict, "nprocs": self.nprocs})
            reactants,products = IO.get_reactants_and_products(rxn)
            parents = []
            if not adsorbates_finished:
                for m in reactants+products:
                    parents.extend(self.adsorbate_fw_dict[m])
            fw = Firework([ts_task],parents=parents,name="TS"+str(i)+"est")
            self.fws.append(fw)

    def rapidfire(self):
        if self.queue:
            rapidfirequeue(self.launchpad,self.fworker,self.qadapter,njobs_queue=self.njobs_queue)
        else:
            rapidfire(self.launchpad,fworker=self.fworker)

    def execute(self):
        if self.slab_path is None: #handle slab
            self.generate_slab()

        #adsorbate optimization
        self.setup_adsorbates()

        #setup transition states
        self.setup_transition_states()

        wf = Workflow(self.fws, name=self.label)
        self.launchpad.add_wf(wf)

        while True:
            self.rapidfire()
