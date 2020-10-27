#!/usr/bin/env python3
from pathlib import Path
'''
####################################################
                    Basic Input
####################################################
'''
####################################################
# do you want to run surface optimization
optimize_slab = True
####################################################
# specify facet orientation, repeats of the slab+ads
# and repeats of the slab_opt unit cell
surface_types_and_repeats = {'fcc111': [(3, 3, 1), (1, 1, 4)],
                             'fcc100': [(3, 4, 1), (1, 1, 4)]}
####################################################
# surface atoms
symbol = 'Cu'
####################################################
# lattice constant
a = 3.6
####################################################
# vacuum in the z direction (Angstrem)
vacuum = 8.0
####################################################
# Quantum Espresso pseudopotantials and exe settings
# for DFT calculations
pseudo_dir = '/home/mgierad/espresso/pseudo'

pseudopotentials = "dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',"\
    + "H='H.pbe-kjpaw_psl.1.0.0.UPF',"\
    + "O='O.pbe-n-kjpaw_psl.1.0.0.UPF'," \
    + "C='C.pbe-n-kjpaw_psl.1.0.0.UPF')"

executable = '/home/mgierad/00_codes/build/q-e-qe-6.4.1/build/bin/pw.x'
####################################################
# Baslam settings
node_packing_count = 48
balsam_exe_settings = {'num_nodes': 1,  # nodes per each balsam job
                       'ranks_per_node': node_packing_count,  # cores per node
                       'threads_per_rank': 1
                       }
calc_keywords = {'kpts': (3, 3, 1),
                 'occupations': 'smearing',
                 'smearing': 'marzari-vanderbilt',
                 'degauss': 0.01,  # Rydberg
                 'ecutwfc': 40,  # Rydberg
                 'nosym': True,  # Allow symmetry breaking during optimization
                 'conv_thr': 1e-11,
                 'mixing_mode': 'local-TF'
                 }
####################################################
# Set up a working directory (this is default)
creation_dir = Path.cwd().as_posix()
####################################################
# filename of the .yaml file with reactions
yamlfile = 'reactions.yaml'
####################################################
# specify the scaling factor to scale the bond distance
# between two atoms taking part in the reaction
scfactor = 1.4
####################################################
# specify the scaling factor to scale the target distance
# i.e. the average bond distance between adsorbate and
# the nearest surface metal atom
scfactor_surface = 1.0
####################################################
# do you want to apply the scfactor_surface to the species 1?
scaled1 = False
####################################################
# do you want to apply scfactor_surface to the species 2?
scaled2 = False
####################################################
