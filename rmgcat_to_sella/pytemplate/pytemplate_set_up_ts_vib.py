#!/usr/bin/env python3

import os

import datetime
from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO

from ase.io import read, write
from ase.constraints import FixAtoms
from ase.vibrations import Vibrations

from sella import Sella

from numpy import floor

geom = '{geom}'
prefix = geom[:2]
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
nimages = {nimages}

start = datetime.datetime.now()

with open(geom + '_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()


# the first and nimages/2 are the same structures - no displacement.
# for inmages = 16 it would be like this
# start...max...start...min...
# floor(nimages/4) is a maximum in one direction (max displacement),
# nimages - floor(nimages/4) is a minimum in the another one (min displacement)
index_forward = int(floor(nimages/4))
index_reverse = int(nimages - index_forward)

atoms = read(os.path.join(prefix, geom + '.xyz'))

# freeze all surface atoms
atoms.set_constraint(FixAtoms(
    [atom.index for atom in atoms if atom.symbol == 'Cu']))
# vibrate only adsorbed species
indices = [atom.index for atom in atoms if atom.symbol != 'Cu']

extra_calc_keywords = dict(
    pseudopotentials={pseudopotentials},
    pseudo_dir='{pseudo_dir}',
    label=geom
)

atoms.calc = EspressoBalsamSocketIO(
    workflow='QE_Socket',
    job_kwargs=balsam_exe_settings,
    **calc_keywords
)

atoms.calc.set(**extra_calc_keywords)

# start vibrations calculations
vib = Vibrations(atoms, indices=indices)
vib.run()
vib.summary()
vib.clean()

# write the first vibration mode to vib.0.traj file (default) - imaginary freq
vib.write_mode(0, nimages=nimages)

# should be one traj file, though
for traj in os.listdir(os.getcwd()):
    if traj.startswith('vib') and traj.endswith('traj'):
        # get forward displacement
        write('forward.xyz', read(traj, index=index_forward))
        # get reverse displacement
        write('reverse.xyz', read(traj, index=index_reverse))

end = datetime.datetime.now()

with open(geom + '_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()
