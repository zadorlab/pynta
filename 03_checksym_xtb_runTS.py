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

from rmgcat_to_sella.ts import create_unique_TS, create_TS_unique_job_files

facetpath  = 'Cu_100'
pytemplate = '/Users/mgierad/.local/lib/python3.7/site-packages/rmgcat_to_sella-0.0.1-py3.7.egg/rmgcat_to_sella/pytemplate/pytemplate_set_up_ts.py'
TSdir      = 'TS_estimate'
create_unique_TS(facetpath, TSdir)
create_TS_unique_job_files(facetpath, TSdir, pytemplate)

bashCommand = os.popen(
    "cd Cu_100/TS_estimate_unique; for i in $(ls | grep py); do sbatch $i; done > ../../submitted_03.txt; cd ../../")
print(bashCommand.read)
