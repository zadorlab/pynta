#!/usr/bin/env python3

import sys
import os
import datetime
from pathlib import Path

from pynta.penalty_fun import AdsorbatePlacer
from xtb.ase.calculator import XTB

from ase.io import read, write

cwd = Path.cwd().as_posix()
path = Path(cwd).parents[2]
sys.path.append(str(path))


geom = '{geom}'
bonds = {bonds}
av_dists_tuple = {av_dists_tuple}
traj_path = '{traj_path}'
creation_dir = '{creation_dir}'
slab = os.path.join(creation_dir, '{slabopt}')
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

adsplacer = AdsorbatePlacer(
    big_slab, ts_estimate, bonds, av_dists_tuple,
    trajectory=traj_path
)

adsplacer.ads_ref.set_calculator(XTB(method="GFN1-xTB"))
opt = adsplacer.optimize()
# visualize initial point of each trajectory
write(geom[:-4] + '_initial.png', read(geom))
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
