#!/usr/bin/env python3

from rmgcat_to_sella.get_slab import GetSlab

surface_type     = '{surface_type}'
symbol           = '{symbol}'
a                = {a}
repeats_surface  = {repeats_surface}
vacuum           = {vacuum}
slab_name        = '{slab_name}'
pseudopotentials = {pseudopotentials}
pseudo_dir       = '{pseudo_dir}'
workflow_name    = slab_name+'00'
executable       = {executable}
balsam_exe_settings = {balsam_exe_settings}
calc_keywords    = {calc_keywords}
creation_dir     = '{creation_dir}'

get_slab = GetSlab(surface_type, symbol, a, repeats_surface, vacuum,
                   slab_name, pseudopotentials, pseudo_dir, executable,
                   balsam_exe_settings,calc_keywords,creation_dir)
get_slab.run_slab_opt()
