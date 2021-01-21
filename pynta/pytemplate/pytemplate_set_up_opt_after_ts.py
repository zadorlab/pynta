#!/usr/bin/env python3
import os

import datetime

from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton

geom = '{geom}'
prefix = geom[:2]
socket_calculator = '{socket_calculator}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
geom_prefix = os.path.join(prefix, geom)

start = datetime.datetime.now()

with open(geom_prefix + '_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()

atoms = read(os.path.join(prefix, geom + '.xyz'))
atoms.set_constraint(FixAtoms([
    atom.index for atom in atoms if atom.position[2] < atoms.cell[2, 2] / 2.
]))

# update balsam_exe_settings with info about a new num_nodes
# balsam_exe_settings['num_nodes'] = {n_kpts}

extra_calc_keywords = dict(
    pseudopotentials={pseudopotentials},
    pseudo_dir='{pseudo_dir}',
    label=geom
)

# kpts={repeats},
# jobs_args='-nk {n_kpts}',

balsamcalc_module = __import__('pynta.balsamcalc', fromlist=[
    socket_calculator])
sock_calc = getattr(balsamcalc_module, socket_calculator)

atoms.calc = sock_calc(
    workflow='QE_Socket',
    job_kwargs=balsam_exe_settings,
    **calc_keywords
)

atoms.calc.set(**extra_calc_keywords)

opt = QuasiNewton(atoms=atoms, trajectory=geom_prefix + '.traj')
opt.run(fmax=0.01, steps=70)
atoms.calc.close()

write_file = os.path.join(geom_prefix + '_final')
write(write_file + '.png', read(geom_prefix + '.traj'))
write(write_file + '.xyz', read(geom_prefix + '.traj'))

end = datetime.datetime.now()

with open(geom_prefix + '_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()
