#!/usr/bin/env python3

import os
import sys
from pathlib import Path


from pathlib import Path
cwd=Path.cwd().as_posix()
path = Path(cwd).parents[2]
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

from xtb.ase.calculator import XTB
adsplacer = AdsorbatePlacer(bigSlab, TS_candidate, bonds, avDists,
                            GFN1(accuracy=0.01,
                                 max_iterations=1000),
                            trajectory=trajPath)
adsplacer.set_calculator(XTB(method="GFN1-xTB"))
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
