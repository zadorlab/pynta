import unittest
import os
from pynta.calculator import *
from ase.io import read

class UtilsTest(unittest.TestCase):

    def test_get_lattice_parameters(self):
        a = get_lattice_parameters('Cu','fcc111','EMT',dict())
        self.assertAlmostEqual(a,3.5898294220923366)

    def test_run_harmonically_forced(self):
        fpath = os.path.abspath(os.path.dirname(__file__))
        atoms = read(os.path.join(os.path.split(fpath)[0],"pynta/testing_data/CH2_Ru0001_dissoc_init.xyz"))
        harm_dict = {"harmonic energy": 0.003004935452450757, "harmonic force": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-0.03715139805848659, -0.5989027528842809, 0.7340574578789104], [0.0, 0.0, 0.0], [0.03579299529037104, 0.9053324017503849, -0.204660094328337]], "atom_bond_potentials": [{"ind1": 36, "ind2": 38, "k": 100.0, "deq": 1.5781941646206141}], "site_bond_potentials": [{"ind": 36, "site_pos": [6.82467097, 3.94022995, 15.819851064573047], "k": 100.0, "deq": 0.0}, {"ind": 38, "site_pos": [6.82448852, 5.51655804, 15.469391200279716], "k": 100.0, "deq": 0.17858443277482008}], "molecule_to_atom_maps": [{"0": 0, "1": 1, "2": 2}], "ase_to_mol_num": {"37": 0, "36": 0, "38": 0}}
        atomsout,Eharm,Fharm = run_harmonically_forced(atoms,harm_dict["atom_bond_potentials"],harm_dict["site_bond_potentials"],36,
        molecule_to_atom_maps=harm_dict["molecule_to_atom_maps"],ase_to_mol_num=harm_dict["ase_to_mol_num"],harm_f_software="TBLite",
                                harm_f_software_kwargs={"method": "GFN1-xTB","verbosity": 0},constraints=["freeze all Ru"])
        
        self.assertAlmostEqual(Eharm,0.00306794366556015)
        
if __name__ == '__main__':
    unittest.main()
