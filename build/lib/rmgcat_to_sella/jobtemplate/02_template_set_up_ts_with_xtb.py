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

from rmgcat_to_sella.ts import TS

slab = '{slab}'
repeats = {repeats}
yamlfile = '{yamlfile}'
facetpath = '{facetpath}'
rotAngle = {rotAngle}
scfactor = {scfactor}
scfactor_surface = {scfactor_surface}
pytemplate_xtb = '{pytemplate_xtb}'
species = ['{sp1}', '{sp2}']
current_dir = os.path.dirname(os.getcwd())
minima_dir = os.path.join(facetpath, 'minima')
scaled1 = {scaled1}
scaled2 = {scaled2}
ts_dir = 'TS_estimate'

ts = TS(facetpath, ts_dir, yamlfile)
ts.copy_minimas_prev_calculated(current_dir, species, minima_dir)
ts.prepare_ts_estimate(slab, repeats, scfactor, scfactor_surface, rotAngle,
                       pytemplate_xtb, species, scaled1, scaled2)

bashCommand = os.popen(
    "cd {facetpath}/TS_estimate/; for i in $(ls -d */); do cd $i; sbatch *py; cd ../ || exit; done > ../../submitted_02.txt; cd ../../")
print(bashCommand.read())
