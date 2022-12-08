import unittest
import os
from pynta.utils import *

class UtilsTest(unittest.TestCase):

    def test_name_to_ase_software(self):
        soft = name_to_ase_software("Espresso")

    def test_name_to_ase_opt(self):
        opt = name_to_ase_opt("MDMin")

if __name__ == '__main__':
    unittest.main()
