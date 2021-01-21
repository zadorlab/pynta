#!/usr/bin/env python3
import os
import datetime

from ase.io import read, write
from ase.constraints import FixAtoms
from sella import Sella

rxn_name = '{rxn_name}'
prefix = '{prefix}'
geom = '{ts_fname}'
facetpath = '{facetpath}'
socket_calculator = '{socket_calculator}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}

trajdir = os.path.join(prefix, prefix + '_' +
                       facetpath + '_' + rxn_name + '.traj')
label = os.path.join(prefix, prefix)

start = datetime.datetime.now()
with open(label + '_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()

ts_atom = read(os.path.join(prefix, geom))
# fix all atoms but not adsorbates
# ts_atom.set_constraint(FixAtoms([
#     atom.index for atom in ts_atom if atom.index < len(ts_atom) - 2
# ]))
# fix bottom half of the slab
ts_atom.set_constraint(FixAtoms([
    atom.index for atom in ts_atom if atom.position[2] < ts_atom.cell[2, 2] / 2.
]))

# update balsam_exe_settings with info about a new num_nodes
# balsam_exe_settings['num_nodes'] = {n_kpts}

extra_calc_keywords = dict(
    pseudopotentials={pseudopotentials},
    pseudo_dir='{pseudo_dir}',
    label=prefix
)

# kpts={repeats},
# jobs_args='-nk {n_kpts}',

balsamcalc_module = __import__('pynta.balsamcalc', fromlist=[
    socket_calculator])
sock_calc = getattr(balsamcalc_module, socket_calculator)

ts_atom.calc = sock_calc(
    workflow='QE_Socket',
    job_kwargs=balsam_exe_settings,
    **calc_keywords
)

ts_atom.calc.set(**extra_calc_keywords)

opt = Sella(ts_atom, order=1, delta0=1e-2, gamma=1e-3, trajectory=trajdir)
opt.run(fmax=0.01, steps=70)
ts_atom.calc.close()

write_dir = os.path.join(prefix, prefix + '_' + facetpath + '_' + rxn_name)
write(write_dir + '_ts_final.png', read(trajdir))
write(write_dir + '_ts_final.xyz', read(trajdir))

end = datetime.datetime.now()
with open(label + '_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()
