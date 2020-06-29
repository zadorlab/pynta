#!/usr/bin/env python3
#SBATCH -J runOptIRC
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

from rmgcat_to_sella.irc import IRC

yamlfile   = '{yamlfile}'
facetpath  = '{facetpath}'
pytemplate = '{pytemplate}'
ts_dir     = 'TS_estimate_unique'
irc_dir    = 'IRC'

irc = IRC(facetpath, ts_dir, yamlfile)
irc.opt_after_IRC(irc_dir, pytemplate)

bashCommand = os.popen(
    'cd {facetpath}/IRC; for i in $(ls -d */); do cd $i; for j in $(ls -d irc*/); do cd $j; sbatch *py; cd ../ || exit; done; cd ../ || exit; done > ../../submitted_05.txt; cd ../../')
print(bashCommand.read())
