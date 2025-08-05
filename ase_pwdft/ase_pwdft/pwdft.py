""" This module defines an ASE interface to PWDFT

https://github.com/ebylaska/PWDFT
"""
import os
import numpy as np

from ase_pwdft.pwdftio.pwdftwriter import write_pwdft_in
from ase_pwdft.pwdftio.pwdftreader import read_pwdft_out
from ase.units import Hartree
from ase.calculators.calculator import FileIOCalculator


class PWDFT(FileIOCalculator):
    """ Class for doing PWDFT calculations.

        calc = PWDFT(label='pwdft', xc='LDA', ecut=70)
    """

    implemented_properties = ['energy', 'forces']
    command = 'pwdft < PREFIX.nwxi > PREFIX.nwxo'
    accepts_bandpath_keyword = True  # To Check
    discard_results_on_any_change = True  # To Check
    echo = True
    # Disable by default
    twodhcurve = False
    lmbfgs = False

    default_parameters = dict()

    def __init__(self, restart=None,
                 ignore_bad_restart_file=FileIOCalculator._deprecated,
                 label='pwdft', atoms=None, **kwargs):
        """ Construct PWDFT-Calculator object"""

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)
        self.calc = None

    def write_input(self, atoms, properties, system_changes):
        """ Write input parameters to files"""
        # Prepare perm and scratch directories
        perm = os.path.abspath(self.parameters.get('perm', 'perm'))
        scratch = os.path.abspath(self.parameters.get('scratch', 'perm'))
        # Ensure directories exist
        if not os.path.exists(perm):
            os.makedirs(perm, exist_ok=True)
        if not os.path.exists(scratch):
            os.makedirs(scratch, exist_ok=True)

        # Determine if system is periodic (surface or crystal)
        pbc = getattr(atoms, 'pbc', None)
        enable_twodhcurve = False
        enable_lmbfgs = False
        if pbc is not None:
            npbc = sum(pbc) if hasattr(pbc, '__iter__') else int(pbc)
            if npbc == 2 or npbc == 3:
                enable_twodhcurve = True
                enable_lmbfgs = True

        with open(self.label + '.nwxi', 'w') as fd:
            write_pwdft_in(
                fd, atoms, properties,
                twod_hcurve=enable_twodhcurve,
                lmbfgs=enable_lmbfgs,
                **self.parameters)

    def read_results(self):
        output = read_pwdft_out(self.label + '.nwxo')
        self.calc = output.calc
        self.results = output.calc.results
