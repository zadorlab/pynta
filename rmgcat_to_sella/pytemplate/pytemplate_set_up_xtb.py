#!/usr/bin/env python3

from xtb.ase.calculator import XTB
from rmgcat_to_sella.prepare_ts_with_xtb import AdsorbatePlacer
import os
import sys
from pathlib import Path
import datetime
from xtb import GFN1
from ase.io import read, write


from pathlib import Path
cwd = Path.cwd().as_posix()
path = Path(cwd).parents[2]
sys.path.append(str(path))


geom = '{geom}'
bonds = {bonds}
av_dists_tuple = {av_dists_tuple}
traj_path = '{traj_path}'
slab = '../../../{slabopt}'
repeats = {repeats}
prefix = '{prefix}'

# label = geom[:2]

adsorbed = read(geom)
slab = read(slab)
big_slab = slab * repeats
nbig_slab = len(big_slab)
ts_estimate = adsorbed[nbig_slab:]

start = datetime.datetime.now()
with open(prefix + '_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
f.close()

adsplacer = AdsorbatePlacer(bigSlab, TS_candidate, bonds, avDists,
                            GFN1(accuracy=0.01,
                                 max_iterations=1000),
                            trajectory=trajPath)
adsplacer.set_calculator(XTB(method="GFN1-xTB"))
opt = adsplacer.optimize()
# visualize end point of each trajectory
write(traj_path[:-5] + '_final.png', read(traj_path))
write(traj_path[:-5] + '_final.xyz', read(traj_path))

end = datetime.datetime.now()
with open(prefix + '_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
f.close()
