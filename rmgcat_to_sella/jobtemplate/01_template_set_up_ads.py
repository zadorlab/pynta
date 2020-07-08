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

facetpath        = '{facetpath}'
slab             = '{slabopt}'
repeats          = {repeats}
yamlfile         = '{yamlfile}'
pytemplate       = '{pytemplate}'
pseudopotentials = {pseudopotentials}
pseudo_dir       = '{pseudo_dir}'

put_adsorbates = Adsorbates(facetpath, slab, repeats, yamlfile)
put_adsorbates.adjacency_to_3d()
put_adsorbates.create_relax_jobs(pytemplate, pseudopotentials, pseudo_dir)

bashCommand = os.popen(
    "cd ./{facetpath}/minima/; for i in $(ls | grep py); do sbatch $i; done > ../../submitted_01.txt; cd ../../")
print(bashCommand.read())
