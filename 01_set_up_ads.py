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

from rmgcat_to_sella.adsorbates import Adsorbates

facetpath        = 'Cu_100'
slab             = 'Cu_100_slab_opt.xyz'
repeats          = (3, 4, 1)
yamlfile         = 'reactions.yaml'
pytemplate       = '/Users/mgierad/.local/lib/python3.7/site-packages/rmgcat_to_sella-0.0.1-py3.7.egg/rmgcat_to_sella/pytemplate/pytemplate_relax_Cu_111_ads.py'
pseudopotentials = dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF', H='H.pbe-kjpaw_psl.1.0.0.UPF', O='O.pbe-n-kjpaw_psl.1.0.0.UPF', C='C.pbe-n-kjpaw_psl.1.0.0.UPF')
pseudo_dir       = '/home/mgierad/espresso/pseudo'

put_adsorbates = Adsorbates(facetpath, slab, repeats, yamlfile)
put_adsorbates.adjacency_to_3d()
put_adsorbates.create_relax_jobs(pytemplate, pseudopotentials, pseudo_dir)

bashCommand = os.popen(
    "cd ./Cu_100/minima/; for i in $(ls | grep py); do sbatch $i; done > ../../submitted_01.txt; cd ../../")
print(bashCommand.read())
