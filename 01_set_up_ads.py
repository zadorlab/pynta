#!/usr/bin/env python3
#SBATCH -J relax_surf_ads
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

from ase.io import read

from rmgcat_to_sella import adjacency_to_3d, create_relax_jobs

slab     = read('Cu_100_slab_opt.xyz')
slab.pbc = [True, True, False]
repeats  = (3, 4, 1)
pytemplate = '/Users/mgierad/.local/lib/python3.7/site-packages/rmgcat_to_sella-0.0.1-py3.7.egg/rmgcat_to_sella/pytemplate/pytemplate_relax_Cu_111_ads.py'
adjacency_to_3d('reactions.yaml', slab, repeats, 'Cu_100')
create_relax_jobs('Cu_100', pytemplate)

bashCommand = os.popen(
    "cd ./Cu_100/minima/; for i in $(ls | grep py); do sbatch $i; done > ../../submitted_01.txt; cd ../../")
print(bashCommand.read())
