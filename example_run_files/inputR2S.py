'''
####################################################
                    Basic Input
####################################################
'''
facetpath              = 'Cu_111'
slabopt                = 'Cu_111_slab_opt.xyz'
yamlfile               = 'reactions.yaml'
repeats                = (3, 3, 1)
rotAngle               = 60
scfactor               = 1.4
scfactor_surface       = 1.0
sp1                    = 'O'
sp2                    = 'H'
scaled1                = False
scaled2                = False
'''
####################################################
                    Scripts
####################################################
'''
SurfaceScript          = '00_gen_surface.py'
SurfaceAdsorbateScript = '01_set_up_ads.py'
TSxtbScript            = '02_set_up_TS_with_xtb.py'
TSScript               = '03_checksym_xtb_runTS.py'
IRCScript              = '04_set_up_irc.py'
IRCoptScript           = '05_set_up_opt_after_irc.py'

