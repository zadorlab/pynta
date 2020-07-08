#!/usr/bin/env python3
#SBATCH -J runIRC
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

slab             = '{slab}'
repeats          = {repeats}
facetpath        = '{facetpath}'
pytemplate_f     = '{pytemplate_f}'
pytemplate_r     = '{pytemplate_r}'
yamlfile         = '{yamlfile}'
ts_dir           = 'TS_estimate'
pseudopotentials = {pseudopotentials}
pseudo_dir       = '{pseudo_dir}'

irc = IRC(facetpath, slab, repeats, ts_dir, yamlfile,
          pseudopotentials, pseudo_dir)
irc.set_up_irc(pytemplate_f, pytemplate_r)

bashCommand = os.popen(
    "cd {facetpath}/IRC/; for i in $(ls | grep 'py'); do sbatch $i; done > ../../submitted_04.txt; cd ../../")
print(bashCommand.read())
