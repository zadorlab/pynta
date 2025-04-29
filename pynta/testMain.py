import unittest
from nose.plugins.attrib import attr
import os
import shutil
from pynta.main import Pynta
from pynta.utils import clean_pynta_path
from fireworks import LaunchPad, Workflow
from fireworks.core.rocket_launcher import rapidfire, launch_rocket


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

        os.chdir(returndir)

if __name__ == '__main__':
    unittest.main()
