import unittest
from nose.plugins.attrib import attr
import os
import shutil
from pynta.main import Pynta, CoverageDependence
from pynta.utils import clean_pynta_path
from ase.io import read
from fireworks import LaunchPad, Workflow
from fireworks.core.rocket_launcher import rapidfire, launch_rocket
from acat.adsorption_sites import SlabAdsorptionSites
import numpy as np
import logging

@attr('functional')
class MainTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """A function that is run ONCE before all unit tests in this class."""
        fpath = os.path.abspath(os.path.dirname(__file__))
        pynta_path = os.path.split(fpath)[0]
        cls.path = os.path.join(pynta_path,"test","pyntatest")
        clean_pynta_path(cls.path)
        cls.launchpath = os.path.join(cls.path,"launches")
        cls.covdep_path = os.path.join(pynta_path,"test","covdeptest")
        cls.launchpath_covdep = os.path.join(cls.covdep_path,"launches")

    @classmethod
    def tearDownClass(cls):
        """A function that is run ONCE after all unit tests in this class."""
        clean_pynta_path(cls.path,save_initial_guess=False) #clean Pynta stuff

        #remove launch files
        fnames = os.listdir(cls.launchpath)
        for n in fnames:
            pv = os.path.join(cls.launchpath,n)
            if os.path.isfile(pv) and n.split(".")[-1] != ".yaml":
                os.remove(pv)
            else:
                shutil.rmtree(pv)
        
        # remove covdep files
        if os.path.exists(os.path.join(cls.covdep_path,"Iterations")):
            shutil.rmtree(os.path.join(cls.covdep_path,"Iterations"))
        if os.path.exists(os.path.join(cls.covdep_path,"Configurations")):
            shutil.rmtree(os.path.join(cls.covdep_path,"Configurations"))
        if os.path.exists(os.path.join(cls.covdep_path,"pairs_datums.json")):
            os.remove(os.path.join(cls.covdep_path,"pairs_datums.json"))
        
    def test_workflow(self):
        returndir = os.getcwd()
        os.chdir(self.launchpath)

        lpad_name=os.path.join(os.path.split(self.path)[0],"my_launchpad.yaml")
        lpad = LaunchPad.from_file(lpad_name)
        wfnames = [lpad.get_wf_summary_dict(wf_id)['name'] for wf_id in lpad.get_wf_ids()]
        name = "Pynta_functional_test_"
        ind = 0
        while name+str(ind) in wfnames:
            ind += 1
        name = name+str(ind)
        
        pyn = Pynta(path=self.path,
                rxns_file=os.path.join(self.path,"rxn_test.yaml"),
                software="XTB",
                surface_type="fcc111",metal="Cu",socket=False,queue=False,a=3.61,
                repeats=(3,3,4),label=name,num_jobs=1,max_num_hfsp_opts=2,
                software_kwargs={"method": "GFN1-xTB","verbosity":0},
                software_kwargs_gas={"method": "GFN1-xTB","verbosity":0},
               TS_opt_software_kwargs={},
               lattice_opt_software_kwargs={},
               slab_path=os.path.join(self.path,"slab.xyz"),
               Eharmtol=1.0, Nharmmin=2,nprocs_harm=4,
               launchpad_path=lpad_name,
               calculate_thermodynamic_references=False,
               )

        pyn.analyze_slab()
        pyn.generate_mol_dict()
        pyn.generate_atom_maps()

        #adsorbate optimization
        pyn.setup_adsorbates()
        
        #setup transition states
        pyn.setup_transition_states()

        wf = Workflow(pyn.fws, name=pyn.label)
        pyn.launchpad.add_wf(wf)

        wf_id = [wf_id for wf_id in pyn.launchpad.get_wf_ids() if pyn.launchpad.get_wf_summary_dict(wf_id)['name']==name][0]

        launch_rocket(pyn.launchpad,fworker=pyn.fworker)
        state = pyn.launchpad.get_wf_summary_dict(wf_id)["state"]

        while state == "RUNNING" or state == "READY" or state == "RESERVED":
            launch_rocket(pyn.launchpad,fworker=pyn.fworker)
            state = pyn.launchpad.get_wf_summary_dict(wf_id)["state"]

        self.assertTrue(pyn.launchpad.get_wf_summary_dict(wf_id)["state"]=="COMPLETED")
        
        self.assertTrue(os.path.exists(os.path.join(self.path,"thermo_library.py")))
        self.assertTrue(os.path.exists(os.path.join(self.path,"reaction_library","reactions.py")))
        self.assertTrue(os.path.exists(os.path.join(self.path,"reaction_library","dictionary.txt")))
        
        os.chdir(returndir)
        
    def test_cov_dep_workflow(self):
        returndir = os.getcwd()
        os.chdir(self.launchpath_covdep)

        lpad_name = os.path.join(os.path.split(self.covdep_path)[0],"my_launchpad.yaml")
        lpad = LaunchPad.from_file(lpad_name)
        wfnames = [lpad.get_wf_summary_dict(wf_id)['name'] for wf_id in lpad.get_wf_ids()]
        name = "CovDep_functional_test_"
        ind = 0
        while name+str(ind) in wfnames:
            ind += 1
        name = name+str(ind)
        
        pynta_path = os.path.join(self.covdep_path,"pyntaxtb") 
        slab = read(os.path.join(pynta_path,"slab.xyz"))
        metal = "Pt"
        facet = "fcc111"
        
        cas = SlabAdsorptionSites(slab, facet,allow_6fold=False,composition_effect=False,
                                    label_sites=True,
                                    surrogate_metal=metal)
        
        sites = cas.get_sites()
        site_adjacency = cas.get_neighbor_site_list()

        covdep = CoverageDependence(self.covdep_path,metal,facet,(3,3,4),pynta_path,software="XTB",software_kwargs={"method": "GFN1-xTB","verbosity":0},
                                    label=name,sites=sites,site_adjacency=site_adjacency,coad_stable_sites={"O=[Pt]":["fcc"],"[Pt]":["fcc"]},adsorbates=["CO[Pt]"],transition_states={"TS0":"9"},
                                    coadsorbates=["O=[Pt]","[Pt]"],
                        frozen_layers=4,fmaxopt=0.05,Ncalc_per_iter=1,launchpad_path=lpad_name,max_iters=1,
                                max_dist=1.0,imag_freq_max=np.inf,num_jobs=1)
        
        covdep.execute(run_pairs=False, run_active_learning=True, launch=False)
        
        wf_id = [wf_id for wf_id in covdep.launchpad.get_wf_ids() if covdep.launchpad.get_wf_summary_dict(wf_id)['name']==name][0]

        launch_rocket(covdep.launchpad,fworker=covdep.fworker)
        state = covdep.launchpad.get_wf_summary_dict(wf_id)["state"]

        while state == "RUNNING" or state == "READY" or state == "RESERVED":
            launch_rocket(covdep.launchpad,fworker=covdep.fworker)
            state = covdep.launchpad.get_wf_summary_dict(wf_id)["state"]

        self.assertTrue(covdep.launchpad.get_wf_summary_dict(wf_id)["state"]=="COMPLETED")
        
        self.assertTrue(os.path.exists(os.path.join(self.covdep_path,"Iterations","1","Ncoad_config_TS0_[Pt].json")))
        self.assertTrue(os.path.exists(os.path.join(self.covdep_path,"Iterations","1","regressor.json")))
        
        os.chdir(returndir)

if __name__ == '__main__':
    unittest.main()
