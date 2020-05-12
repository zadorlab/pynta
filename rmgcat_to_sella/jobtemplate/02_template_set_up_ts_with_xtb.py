#!/usr/bin/env python3
#SBATCH -J set_up_TS_with_xtb
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p day-long-cpu
#SBATCH -t 1-00:00:00
#SBATCH -e %x.err
#SBATCH -o %x.out

from rmgcat_to_sella.ts import genTSestimate, set_up_penalty_xtb

from ase.io import read

import os

slab           = read('{slabopt}')
repeats        = {repeats}
yamlfile       = '{yamlfile}'
facetpath      = '{facetpath}'
rotAngle       = {rotAngle}
scfactor       = {scfactor}
pytemplate_xtb = '{pytemplate_xtb}'
path           = os.path.join(facetpath, 'TS_estimate')
sp1 = '{sp1}'
sp2 = '{sp2}'

genTSestimate(slab, repeats, yamlfile, facetpath, rotAngle, scfactor)
set_up_penalty_xtb(path, pytemplate_xtb, repeats, sp1, sp2)

bashCommand = os.popen(
    "cd {facetpath}/TS_estimate/; for i in $(ls -d */); do cd $i; sbatch *py; cd ../ || exit; done > ../../submitted_02.txt; cd ../../")
print(bashCommand.read())
