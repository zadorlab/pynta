from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.vibrations import Vibrations
import os
from ase.io import read, write
from numpy import floor

nimages = 30
# the first and nimages/2 are the same structures - no displacement.
# for inmages = 16 it would be like this
# start...max...start...min...
# floor(nimages/4) is a maximum in one direction (max displacement),
# nimages - floor(nimages/4) is a minimum in the another one (min displacement)
index_forward = int(floor(nimages/4))
index_reverse = int(nimages - index_forward)

geom = read('{traj}')
geom.calc = (EMT())
print(geom)
vib = Vibrations(geom)
vib.run()
vib.summary()
vib.clean()

# write the last vibration mode to traj file
vib.write_mode(0, nimages=nimages)

# should be one traj file, though
for traj in os.listdir(os.getcwd()):
    if traj.endswith('traj'):
        # get forward displacement
        write('forward.xyz', read(traj, index=index_forward))
        # get reverse displacement
        write('reverse.xyz', read(traj, index=index_reverse))

f = read('forward.xyz')
r = read('reverse.xyz')
# print(f.get_distance(0, 1), r.get_distance(0, 1))
