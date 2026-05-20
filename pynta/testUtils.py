import unittest
import os
from pynta.utils import *
from ase.calculators.mixing import SumCalculator
from ase.calculators.espresso import Espresso
from ase.calculators.emt import EMT

class UtilsTest(unittest.TestCase):

    def test_name_to_ase_software(self):
        soft = name_to_ase_software("Espresso")

    def test_name_to_ase_opt(self):
        opt = name_to_ase_opt("MDMin")
    
    def test_to_ase_software(self):
        software_kwargs=[{'kpts': (4, 4, 1), 'tprnfor': True, 'occupations': 'smearing',
                            'smearing':  'gauss', 'input_dft': 'BEEF-vdW',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-6, 'mixing_mode': 'local-TF',
                            "pseudopotentials": {"Pt": 'Pt.pbe-n-kjpaw_psl.1.0.0.UPF',"H": 'H.pbe-kjpaw_psl.1.0.0.UPF',"O": 'O.pbe-n-kjpaw_psl.1.0.0.UPF',"C": 'C.pbe-n-kjpaw_psl.1.0.0.UPF',"N": 'N.pbe-n-kjpaw_psl.1.0.0.UPF'},
                            "profile": {"type": "EspressoProfile",
                                        "command": 'srun /opt/custom/espresso/6.6_nostress/bin/pw.x < PREFIX.pwi > PREFIX.pwo',
                                        "pseudo_dir": ".",
                            }
                            },dict()]
        soft = to_ase_software(["Espresso","EMT"],software_kwargs)
        assert isinstance(soft,SumCalculator)
        assert isinstance(soft.mixer.calcs[0],Espresso)
        assert isinstance(soft.mixer.calcs[1],EMT)

if __name__ == '__main__':
    unittest.main()
