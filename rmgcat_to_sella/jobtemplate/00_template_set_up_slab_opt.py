#!/usr/bin/env python3
#SBATCH -J slab_opt
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1gb
#SBATCH -p day-long-cpu
#SBATCH -t 1-00:00:00
#SBATCH -e %x.err
#SBATCH -o %x.out

from rmgcat_to_sella.get_slab import GetSlab
import os
import sys

submitDir = os.environ['SLURM_SUBMIT_DIR']
os.chdir(submitDir)
sys.path.append(os.getcwd())

surface_type     = '{surface_type}'
symbol           = '{symbol}'
a                = {a}
repeats_surface  = {repeats_surface}
vacuum           = {vacuum}
slab_name        = '{slab_name}'
pseudopotentials = {pseudopotentials}


get_slab = GetSlab(surface_type, symbol, a, repeats_surface, vacuum,
                   slab_name, pseudopotentials)
get_slab.run_slab_opt()
