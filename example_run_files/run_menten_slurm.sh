#!/usr/bin/env python3
#SBATCH -J main_job_Cu_111
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p week-long-cpu
#SBATCH -t 7-00:00:00
#SBATCH -e %x.err
#SBATCH -o %x.out

import os
import sys
submitDir = os.environ['SLURM_SUBMIT_DIR']
os.chdir(submitDir)
sys.path.append(os.getcwd())
import time
import inputR2S
from rmgcat_to_sella.main import WorkFlow

workflow = WorkFlow()
workflow.gen_job_files()
workflow.execute()
