'''
####################################################
                    Basic Input
####################################################
'''
####################################################
# specify the name of the main directory with calculations
facetpath        = 'Cu_100'
####################################################
# do you want to run surface optimization
optimize_slab    = True
####################################################
# specify name of the slab
slab_name        = 'Cu_100_slab_opt'
####################################################
# specify facet orientation
surface_type     = 'fcc100'
####################################################
# surface atoms
symbol           = 'Cu'
####################################################
# lattice constant
a                = 3.6
####################################################
# vacuum in the z direction (Angstrem)
vacuum           = 8.0
####################################################
# filename of the optimized surface slab
slabopt          = 'Cu_100_slab_opt.xyz'
####################################################
# Quantum Espresso pseudopotantials for DFT calculations
#pseudo_dir = '/home/mgierad/espresso/pseudo'
pseudo_dir = '/home/brossdh/espresso/pseudo'
pseudopotentials = "dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF', H='H.pbe-kjpaw_psl.1.0.0.UPF', O='O.pbe-n-kjpaw_psl.1.0.0.UPF', C='C.pbe-n-kjpaw_psl.1.0.0.UPF')"
executable = '\'/soft/applications/quantum_espresso/6.4/bin/pw.x\''
balsam_exe_settings = {'num_nodes': 1,
        'ranks_per_node': 16,
        'threads_per_rank': 4,
        'threads_per_core': 1
        }

calc_keywords= {'kpts':(3, 3, 1),
        'occupations':'smearing',
        'smearing':'marzari-vanderbilt',
        'degauss':0.01,  # Rydberg
        'ecutwfc':40,  # Rydberg
        'nosym':True,  # Allow symmetry breaking during optimization
        'conv_thr':1e-11,
        'mixing_mode':'local-TF'
        }

####################################################
# filename of the .yaml file with reactions
yamlfile         = 'reactions.yaml'
####################################################
# specify repeats of the surface in (x, y, z) direction
repeats_surface  = (1, 1, 4)
####################################################
# specify repeats of the surface in (x, y, z) direction
repeats          = (3, 4, 1)
####################################################
# specify the angle of TS estimate addut rotation
rotAngle         = 60
####################################################
# specify the scaling factor to scale the bond distance
# between two atoms taking part in the reaction
scfactor         = 1.4
####################################################
# specify the scaling factor to scale the target distance
# i.e. the average bond distance between adsorbate and
# the nearest surface metal atom
scfactor_surface = 1.0
####################################################
# species list
species_list     = ['O', 'H'] 
####################################################
# do you want to apply the scfactor_surface to the species 1?
scaled1          = False
####################################################
# do you want to apply scfactor_surface to the species 2?
scaled2          = False
####################################################
'''
####################################################
                    Scripts
####################################################
'''
slab_opt_script = '00_set_up_slab_opt.py'
SurfaceAdsorbateScript = '01_set_up_ads.py'
TSxtbScript = '02_set_up_TS_with_xtb.py'
TSScript = '03_checksym_xtb_runTS.py'
IRCScript = '04_set_up_irc.py'
IRCoptScript = '05_set_up_opt_after_irc.py'
