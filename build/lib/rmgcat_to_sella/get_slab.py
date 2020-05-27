import os
import shutil
import getpass

from ase.build import fcc111, fcc211
from ase.io import read, write
from ase.calculators.espresso import Espresso
from ase.calculators.socketio import SocketIOCalculator

from ase.optimize import LBFGS

from gpaw import GPAW, PW

from sella import Sella

import os

# os.environ['GPAW_SETUP_PATH']
# avaiable options for slab build are: fcc100, fcc110, bcc100, bcc110, bcc111, fcc111, hcp0001, hcp10m10, diamond100, diamond111

def get_slab_fcc_111(symbol, a, vacuum, size, name, ext):
	slab = fcc111(symbol, size, a, vacuum, orthogonal=False, periodic=True)
	
	#setting up calculator
	# calc = GPAW(xc='PBE', mode = 'pw', kpts=(4, 4, 4))
	# slab.set_calculator(calc)
	# dyn = LBFGS(slab, trajectory = 'slab_Cu.traj')
	# dyn.run(fmax=0.01)

	unixsocket = name
	socketpath = f'/tmp/ipi_{unixsocket}'
	if os.path.exists(socketpath):
		os.remove(socketpath)

	# jobdir = os.path.join(unixsocket, 'opt')
	if os.path.exists(unixsocket):
		shutil.rmtree(unixsocket)
	os.makedirs(unixsocket)

	label = os.path.join(unixsocket, 'Cu_111_slab_opt')

	espresso = Espresso(command='mpirun -np 8 /Users/mgierad/00_SANDIA_WORK/03_codes/build/q-e-qe-6.4.1/bin/pw.x -inp PREFIX.pwi --ipi {{unixsocket}}:UNIX > PREFIX.pwo'
								.format(unixsocket=unixsocket),
						label=label,
						pseudopotentials=dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
											  ),
						pseudo_dir='/Users/mgierad/00_SANDIA_WORK/03_codes/build/q-e-qe-6.4.1/pseudoPOT',
						kpts=(3, 3, 1),
						occupations='smearing',
						smearing='marzari-vanderbilt',
						degauss=0.01,  # Rydberg
						ecutwfc=40,  # Rydberg
						nosym=True,  # Allow symmetry breaking during optimization
						conv_thr=1e-11,
						mixing_mode='local-TF',
						)

	with SocketIOCalculator(espresso, unixsocket=unixsocket) as calc:
		slab.calc = calc
		opt = Sella(slab, order=0, delta0=1e-2, trajectory=label + '.traj')
		opt.run(fmax=0.01)


	slab.get_potential_energy()
	slab.get_forces()
	write(name + '.' + ext, slab)


def get_slab_fcc_211(symbol, a, vacuum, size, name, ext):
	slab = fcc211(symbol, size, a, vacuum, orthogonal=True)
	
	'''
	Using gipaw - fast calculations
	'''
	#setting up calculator
	# calc = GPAW(xc='PBE', mode = 'pw')
	# slab.set_calculator(calc)
	# dyn = LBFGS(slab, trajectory = 'slab_Cu.traj')
	# dyn.run(fmax=0.01)
	# slab.get_potential_energy()
	# slab.get_forces()

	unixsocket = name
	socketpath = f'/tmp/ipi_{unixsocket}'
	if os.path.exists(socketpath):
		os.remove(socketpath)

	# jobdir = os.path.join(unixsocket, 'opt')
	if os.path.exists(unixsocket):
		shutil.rmtree(unixsocket)
	os.makedirs(unixsocket)

	label = os.path.join(unixsocket, 'Cu_211_slab_opt')

	espresso = Espresso(command='mpirun -np 8 /Users/mgierad/00_SANDIA_WORK/03_codes/build/q-e-qe-6.4.1/bin/pw.x -inp PREFIX.pwi --ipi {{unixsocket}}:UNIX > PREFIX.pwo'
								.format(unixsocket=unixsocket),
						label=label,
						pseudopotentials=dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
											  ),
						pseudo_dir='/Users/mgierad/00_SANDIA_WORK/03_codes/build/q-e-qe-6.4.1/pseudoPOT',
						kpts=(3, 3, 1),
						occupations='smearing',
						smearing='marzari-vanderbilt',
						degauss=0.01,  # Rydberg
						ecutwfc=40,  # Rydberg
						nosym=True,  # Allow symmetry breaking during optimization
						conv_thr=1e-11,
						mixing_mode='local-TF',
						)

	with SocketIOCalculator(espresso, unixsocket=unixsocket) as calc:
		slab.calc = calc
		opt = Sella(slab, order=0, delta0=1e-2, trajectory=label + '.traj')
		opt.run(fmax=0.01)


	ener = slab.get_potential_energy()
	force = slab.get_forces()
	write(name + '.' + ext, slab)
