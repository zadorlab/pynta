#!/usr/bin/env python3
#SBATCH -J {geomName}_xtb
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p day-long-cpu
#SBATCH -t 1-00:00:00
#SBATCH -e %x.err
#SBATCH -o %x.out

from rmgcat_to_sella.prepare_ts_with_xtb import AdsorbatePlacer

from ase.io import read, write

from xtb import GFN1

import datetime

geom = '{geom}'
bonds = {bonds}
avDists = {avDists}
trajPath = '{trajPath}'
slab = '../../../Cu_111_slab_opt.xyz'
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

end = datetime.datetime.now()
with open(prefix + '_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
f.close()
