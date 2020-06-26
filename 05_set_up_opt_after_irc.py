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

# import inputR2S

from rmgcat_to_sella.irc import optAfterIRC

facetpath  = 'Cu_100'
IRCPath    = 'IRC'
pytemplate = '/Users/mgierad/.local/lib/python3.7/site-packages/rmgcat_to_sella-0.0.1-py3.7.egg/rmgcat_to_sella/pytemplate/pytemplate_set_up_opt_irc.py'
path       = os.path.join(facetpath, IRCPath)

optAfterIRC(path, pytemplate)

bashCommand = os.popen(
    'cd Cu_100/IRC; for i in $(ls -d */); do cd $i; for j in $(ls -d irc*/); do cd $j; sbatch *py; cd ../ || exit; done; cd ../ || exit; done > ../../submitted_05.txt; cd ../../')
print(bashCommand.read())
