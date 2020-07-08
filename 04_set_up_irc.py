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

slab             = 'Cu_100_slab_opt.xyz'
repeats          = (3, 4, 1)
facetpath        = 'Cu_100'
pytemplate_f     = '/Users/mgierad/.local/lib/python3.7/site-packages/rmgcat_to_sella-0.0.1-py3.7.egg/rmgcat_to_sella/pytemplate/pytemplate_set_up_irc_f.py'
pytemplate_r     = '/Users/mgierad/.local/lib/python3.7/site-packages/rmgcat_to_sella-0.0.1-py3.7.egg/rmgcat_to_sella/pytemplate/pytemplate_set_up_irc_r.py'
yamlfile         = 'reactions.yaml'
ts_dir           = 'TS_estimate'
pseudopotentials = dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF', H='H.pbe-kjpaw_psl.1.0.0.UPF', O='O.pbe-n-kjpaw_psl.1.0.0.UPF', C='C.pbe-n-kjpaw_psl.1.0.0.UPF')
pseudo_dir       = '/home/mgierad/espresso/pseudo'

irc = IRC(facetpath, slab, repeats, ts_dir, yamlfile,
          pseudopotentials, pseudo_dir)
irc.set_up_irc(pytemplate_f, pytemplate_r)

bashCommand = os.popen(
    "cd Cu_100/IRC/; for i in $(ls | grep 'py'); do sbatch $i; done > ../../submitted_04.txt; cd ../../")
print(bashCommand.read())
