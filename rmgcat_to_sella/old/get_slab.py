from ase.build import fcc111
from ase.io import read, write

from gpaw import GPAW, PW

import os

os.environ['GPAW_SETUP_PATH']
# avaiable options for slab build are: fcc100, fcc110, bcc100, bcc110, bcc111, fcc111, hcp0001, hcp10m10, diamond100, diamond111

def get_slab(symbol, a, vacuum, size, name, ext):
	slab = fcc111(symbol, size, a, vacuum, orthogonal=False, periodic=True)
	
	#setting up calculator
	calc = GPAW(mode = 'pw', kpts=(4, 4, 4))

	slab.set_calculator(calc)
	ener = slab.get_potential_energy()
	force = slab.get_forces()
	dot = '.'
	write(name + dot + ext, slab)
