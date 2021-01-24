#!/usr/bin/env python3
from pynta.get_slab import GetSlab

socket_calculator = '{socket_calculator}'
surface_type = '{surface_type}'
metal_atom = '{metal_atom}'
a = {a}
repeats_surface = {repeats_surface}
vacuum = {vacuum}
slab_name = '{slab_name}'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
workflow_name = slab_name + '00'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'

get_slab = GetSlab(socket_calculator, surface_type, metal_atom, a,
                   repeats_surface, vacuum, slab_name, pseudopotentials,
                   pseudo_dir, balsam_exe_settings, calc_keywords,
                   creation_dir)
get_slab.run_slab_opt()
