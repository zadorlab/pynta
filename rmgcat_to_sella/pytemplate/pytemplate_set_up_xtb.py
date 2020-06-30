#!/usr/bin/env python3
#SBATCH -J {geomName}_xtb
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

submitDir = os.environ['SLURM_SUBMIT_DIR']
os.chdir(submitDir)
path = Path(submitDir).parents[2]
sys.path.append(str(path))
# import inputR2S

import datetime
from xtb import GFN1
from ase.io import read, write
from rmgcat_to_sella.prepare_ts_with_xtb import AdsorbatePlacer


geom = '{geom}'
bonds = {bonds}
avDists = {avDists}
trajPath = '{trajPath}'
slab = '../../../{slabopt}'
repeats = {repeats}
prefix = '{prefix}'

# label = geom[:2]

adsorbed = read(geom)
slab = read(slab)
bigSlab = slab * repeats
nbigSlab = len(bigSlab)
TS_candidate = adsorbed[nbigSlab:]

start = datetime.datetime.now()
with open(prefix + '_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
f.close()

adsplacer = AdsorbatePlacer(bigSlab, TS_candidate, bonds, avDists,
                            GFN1(accuracy=0.01,
                                 max_iterations=1000),
                            trajectory=trajPath)
opt = adsplacer.optimize()
# visualize end point of each trajectory
write(trajPath[:-5] + '_final.png', read(trajPath))
write(trajPath[:-5] + '_final.xyz', read(trajPath))

end = datetime.datetime.now()
with open(prefix + '_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
f.close()
