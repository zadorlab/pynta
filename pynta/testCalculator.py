import unittest
import os
from pynta.calculator import *

class UtilsTest(unittest.TestCase):

    def test_get_lattice_parameters(self):
        a = get_lattice_parameters('Cu','fcc111','EMT',dict())
        self.assertAlmostEqual(a,3.5898294220923366)

if __name__ == '__main__':
    unittest.main()
