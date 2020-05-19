#!/usr/bin/env python3
#SBATCH -J runOptIRC
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

from rmgcat_to_sella.irc import optAfterIRC

facetpath  = '{facetpath}'
IRCPath    = 'IRC'
pytemplate = '{pytemplate}'
path       = os.path.join(facetpath, IRCPath)

optAfterIRC(path, pytemplate)

bashCommand = os.popen(
    'cd {facetpath}/IRC; for i in $(ls -d */); do cd $i; for j in $(ls -d irc*/); do cd $j; sbatch *py; cd ../ || exit; done; cd ../ || exit; done > ../../submitted_05.txt; cd ../../')
print(bashCommand.read())
