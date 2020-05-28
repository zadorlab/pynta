#!/usr/bin/env python3
#SBATCH -J set_up_TS_with_xtb
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1gb
#SBATCH -p day-long-cpu
#SBATCH -t 1-00:00:00
#SBATCH -e %x.err
#SBATCH -o %x.out

import os
import sys
submitDir = os.environ['SLURM_SUBMIT_DIR']
os.chdir(submitDir)
sys.path.append(os.getcwd())

# import inputR2S

from rmgcat_to_sella.ts import genTSestimate, set_up_penalty_xtb, copyMinimasPrevCalculated

from ase.io import read

slab             = read('{slabopt}')
slabopt          = '{slabopt}'
repeats          = {repeats}
yamlfile         = '{yamlfile}'
facetpath        = '{facetpath}'
rotAngle         = {rotAngle}
scfactor         = {scfactor}
scfactor_surface = {scfactor_surface}
pytemplate_xtb   = '{pytemplate_xtb}'
path             = os.path.join(facetpath, 'TS_estimate')
sp1              = '{sp1}'
sp2              = '{sp2}'
checkMinimaDir   = os.path.dirname(os.getcwd())
dstDir           = os.path.join(facetpath, 'minima')
scaled1          = {scaled1}
scaled2          = {scaled2}

copyMinimasPrevCalculated(checkMinimaDir, sp1, sp2, dstDir, facetpath)
genTSestimate(slab, repeats, yamlfile, facetpath, rotAngle, scfactor)
set_up_penalty_xtb(path, pytemplate_xtb, repeats, slabopt, sp1, sp2, scfactor_surface, scaled1, scaled2)

bashCommand = os.popen(
    "cd {facetpath}/TS_estimate/; for i in $(ls -d */); do cd $i; sbatch *py; cd ../ || exit; done > ../../submitted_02.txt; cd ../../")
print(bashCommand.read())
