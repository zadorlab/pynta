from pynta.tasks import *
from fireworks import LaunchPad
import ase.build
from ase.io import read, write
import os
import time
from pynta.io import IO

class Pynta:
    def __init__(self,path,rxns_file,surface_type,metal,a=3.6,vaccum=8.0,
    repeats=(3,3,1),slab_path=None,software="espresso",socket=False,
    software_kwargs={'kpts': (3, 3, 1), 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-11, 'mixing_mode': 'local-TF',
                            "pseudopotentials": {"Cu": 'Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',"H": 'H.pbe-kjpaw_psl.1.0.0.UPF',"O": 'O.pbe-n-kjpaw_psl.1.0.0.UPF',"C": 'C.pbe-n-kjpaw_psl.1.0.0.UPF',"N": 'N.pbe-n-kjpaw_psl.1.0.0.UPF',
                            }}):
        self.surface_type = surface_type
        launchpad = LaunchPad()
        launchpad.reset('', require_password=False)
        self.launchpad = launchpad
        self.slab_path = slab_path
        self.vaccum = vaccum
        self.a = a
        self.software = software
        self.repeats = repeats
        self.path = os.getcwd() if path is None else path
        self.facet = metal + surface_type
        self.fws = []
        self.rxns = IO.open_yaml_file(rxns_file)
        self.metal = metal
        self.adsorbate_fw_dict = dict()
        self.software_kwargs = software_kwargs

    def generate_slab(self):
        slab_type = getattr(ase.build,self.surface_type)
        slab = slab_type(self.metal,self.repeats,self.a,self.vaccum)
        slab.pbc = (True, True, False)
        write(os.path.join(self.path,"slab_init.xyz"),slab)
        self.slab_path = os.path.join(self.path,"slab.xyz")
        fwslab = optimize_firework(os.path.join(self.path,"slab_init.xyz"),self.software,"slab",
            opt_method="BFGSLineSearch",socket=self.socket,software_kwargs=self.software_kwargs,
            run_kwargs={"fmax" : 0.01})
        return fwslab

    def setup_adsorbates(self):
        put_adsorbates = Adsorbates(self.path, self.slab_path, self.repeats, self.rxns_file)
        adsorbate_dict = put_adsorbates.adjacency_to_3d()
        for adsname,adsorbate in adsorbate_dict.items():
            xyzs = []
            optfws = []
            for prefix,structure in adsorbate.items():
                os.makedirs(os.path.join(self.path,"Adsorbates",adsname,prefix))
                write(os.path.join(self.path,"Adsorbates",adsname,prefix,prefix+"_init.xyz"),structure)
                xyz = os.path.join(self.path,"Adsorbates",adsname,prefix,prefix+".xyz")
                xyzs.append(xyz)
                fwopt = optimize_firework(os.path.join(self.path,"Adsorbates",adsname,prefix,prefix+"_init.xyz"),
                    self.software,prefix,
                    opt_method="QuasiNewton",socket=self.socket,software_kwargs=self.software_kwargs,
                    run_kwargs={"fmax" : 0.01, "steps" : 70},parents=ads_parents,constraints=["freeze slab"])

                fwvib = vibrations_firework(os.path.join(self.path,"Adsorbates",adsname,prefix,prefix+".xyz"),self.software,
                    prefix,software_kwargs=self.software_kwargs,parents=ads_parents+[fwopt],constraints=["freeze slab"])
                optfws.append(fwopt)

            vib_obj_dict = {"software": self.software, "label": prefix, "software_kwargs": self.software_kwargs,
                "constraints": ["freeze slab"]}
            ctask = MolecularCollect(xyzs,False,[vibrations_firework], [vib_obj_dict],
                    ["vib.json"],[False])
            cfw = Firework([ctask],parents=optfws)
            self.adsorbate_fw_dict[adsname] = cfw
            self.fws.extend(optfws+[cfw])


    def setup_transition_states(self):
        opt_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                "run_kwargs": {"fmax" : 0.01, "steps" : 70},"constraints": ["freeze slab"],"sella":True,"order":1}
        vib_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                "constraints": constraints}
        TSnudge_obj_dict = {"software":self.software,"label":"prefix","socket":self.socket,"software_kwargs":self.software_kwargs,
                "run_kwargs":{"fmax" : 0.01, "steps" : 70},"constraints":["freeze slab"],"opt_method":"QuasiNewton"}
        for i,rxn in enumerate(self.rxns):
            ts_path = os.path.join(self.path,"TS"+str(i))
            os.makedirs(ts_path)

            ts_task = MolecularTSEstimate(rxn,ts_path,self.slab_path,os.path.join(self.path,"Adsorbates"),
                self.rxns_file,self.repeats,self.path.self.metal,out_path=ts_path,scfactor=1.4,scfactor_surface=1.0,
                    scaled1=True,scaled2=False,spawn_jobs=True,opt_obj_dict=opt_obj_dict,vib_obj_dict=vib_obj_dict,
                    TSnudge_obj_dict=TSnudge_obj_dict)
            reactants,products = IO.get_reactants_and_products(rxn)
            parents = []
            for m in reactants+products:
                parents.append(self.adsorbate_fw_dict[m])
            fw = Firework([ts_task],parents=parents)
            self.fws.append(fw)

    def execute(self):
        if self.slab_path is None: #handle slab
            fwslab = self.generate_slab()
            wfslab = Workflow([fwslab], name="slab")
            lpad.add_wf(wfslab)
            rapidfire(lpad)
            while not os.path.exists(self.slab_path): #wait until slab optimizes, this is required anyway and makes the rest of the code simpler
                time.sleep(1)

        slab = read(self.slab_path)

        #adsorbate optimization
        self.setup_adsorbates()

        #setup transition states
        self.setup_transition_states()

        wf = Workflow(self.fws, name="pynta")
        lpad.add_wf(wf)
        rapidfire(lpad)
