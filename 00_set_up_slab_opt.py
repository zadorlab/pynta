#!/usr/bin/env python3
#SBATCH -J slab_opt
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem=4gb
#SBATCH -p day-long-cpu
#SBATCH -t 1-00:00:00
#SBATCH -e %x.err
#SBATCH -o %x.out

import os
import sys
submitDir = os.environ['SLURM_SUBMIT_DIR']
os.chdir(submitDir)
sys.path.append(os.getcwd())

from rmgcat_to_sella.get_slab import GetSlab

surface_type     = 'fcc100'
symbol           = 'Cu'
a                = 3.6
repeats_surface  = (1, 1, 4)
vacuum           = 8.0
slab_name        = 'Cu_100_slab_opt'
pseudopotentials = dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF', H='H.pbe-kjpaw_psl.1.0.0.UPF', O='O.pbe-n-kjpaw_psl.1.0.0.UPF', C='C.pbe-n-kjpaw_psl.1.0.0.UPF')

get_slab = GetSlab(surface_type, symbol, a, repeats_surface, vacuum,
                   slab_name, pseudopotentials)
get_slab.run_slab_opt()
