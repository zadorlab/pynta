#!/usr/bin/env python3
#SBATCH -J {geom_name}_xtb
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1gb
#SBATCH -p day-long-cpu
#SBATCH -t 1-00:00:00
#SBATCH -e %x.err
#SBATCH -o %x.out

import os
import sys
from pathlib import Path
import datetime
from xtb import GFN1
from ase.io import read, write

submitDir = os.environ['SLURM_SUBMIT_DIR']
os.chdir(submitDir)
path = Path(submitDir).parents[2]
sys.path.append(str(path))

from rmgcat_to_sella.prepare_ts_with_xtb import AdsorbatePlacer


geom     = '{geom}'
bonds    = {bonds}
av_dists_tuple  = {av_dists_tuple}
traj_path = '{traj_path}'
slab     = '../../../{slabopt}'
repeats  = {repeats}
prefix   = '{prefix}'

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

adsplacer = AdsorbatePlacer(big_slab, ts_estimate, bonds, av_dists_tuple,
                            GFN1(accuracy=0.01,
                                max_iterations=1000),
                                trajectory=traj_path)
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
