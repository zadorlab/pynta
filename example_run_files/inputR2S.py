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
pseudo_dir = '/home/mgierad/espresso/pseudo'
pseudopotentials = "dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF', H='H.pbe-kjpaw_psl.1.0.0.UPF', O='O.pbe-n-kjpaw_psl.1.0.0.UPF', C='C.pbe-n-kjpaw_psl.1.0.0.UPF')"
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
# species 1
sp1              = 'O'
####################################################
# species 2
sp2              = 'H'
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
