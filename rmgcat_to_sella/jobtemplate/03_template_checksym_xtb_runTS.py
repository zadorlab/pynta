#!/usr/bin/env python3
#SBATCH -J runTS
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p day-long-cpu
#SBATCH -t 1-00:00:00
#SBATCH -e %x.err
#SBATCH -o %x.out

import os
import sys
submitDir = os.environ['SLURM_SUBMIT_DIR']
os.chdir(submitDir)
sys.path.append(os.getcwd())

import inputR2S

from rmgcat_to_sella.ts import create_unique_TS, create_TS_unique_job_files

facetpath  = '{facetpath}'
pytemplate = '{pytemplate}'
TSdir      = 'TS_estimate'
create_unique_TS(facetpath, TSdir)
create_TS_unique_job_files(facetpath, TSdir, pytemplate)

bashCommand = os.popen(
    "cd {facetpath}/TS_estimate_unique; for i in $(ls | grep py); do sbatch $i; done > ../../submitted_03.txt; cd ../../")
print(bashCommand.read)
