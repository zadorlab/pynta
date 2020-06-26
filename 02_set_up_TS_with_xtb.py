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

from rmgcat_to_sella.ts import TS

from ase.io import read

slab             = 'Cu_100_slab_opt.xyz'
repeats          = (3, 4, 1)
yamlfile         = 'reactions.yaml'
facetpath        = 'Cu_100'
rotAngle         = 60
scfactor         = 1.4
scfactor_surface = 1.0
pytemplate_xtb   = '/Users/mgierad/.local/lib/python3.7/site-packages/rmgcat_to_sella-0.0.1-py3.7.egg/rmgcat_to_sella/pytemplate/pytemplate_set_up_xtb.py'
ts_estimate_path = os.path.join(facetpath, 'TS_estimate')
species          = ['O', 'H']
current_dir      = os.path.dirname(os.getcwd())
minima_dir       = os.path.join(facetpath, 'minima')
scaled1          = False
scaled2          = False

ts = TS(ts_estimate_path, slab, repeats, yamlfile, facetpath, rotAngle,
        scfactor, scfactor_surface, scaled1, scaled2)
ts.copy_minimas_prev_calculated(current_dir, species, minima_dir)
ts.prepare_ts_estimate()
ts.set_up_penalty_xtb(pytemplate_xtb, slabopt, species)

bashCommand = os.popen(
    "cd Cu_100/TS_estimate/; for i in $(ls -d */); do cd $i; sbatch *py; cd ../ || exit; done > ../../submitted_02.txt; cd ../../")
print(bashCommand.read())
