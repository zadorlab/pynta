#!/usr/bin/env python3
#SBATCH -J runTS
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

slab             = '{slab}'
repeats          = {repeats}
yamlfile         = '{yamlfile}'
facetpath        = '{facetpath}'
pytemplate       = '{pytemplate}'
ts_dir           = 'TS_estimate'
pseudopotentials = {pseudopotentials}
pseudo_dir       = '{pseudo_dir}'

ts = TS(facetpath, slab, ts_dir, yamlfile, repeats)
ts.create_unique_TS()
ts.create_TS_unique_job_files(pytemplate, pseudopotentials, pseudo_dir)

bashCommand = os.popen(
    "cd {facetpath}/TS_estimate_unique; for i in $(ls | grep py); do sbatch $i; done > ../../submitted_03.txt; cd ../../")
print(bashCommand.read)
