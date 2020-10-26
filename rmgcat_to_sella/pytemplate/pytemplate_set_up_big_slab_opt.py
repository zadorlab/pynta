#!/usr/bin/env python3
import os

from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO

from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton

slab_name = '{slab_name}'
big_slab_name = '{facetpath}' + '_big_slab_opt'
repeats = {repeats}
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'

atoms = read(slab_name + '.xyz') * repeats
atoms.set_constraint(FixAtoms([
    atom.index for atom in atoms if atom.position[2] < atoms.cell[2, 2] / 2.
]))

extra_calc_keywords = dict(
    pseudopotentials={pseudopotentials},
    pseudo_dir='{pseudo_dir}',
    label=label
)

atoms.calc = EspressoBalsamSocketIO(
    workflow='QE_Socket',
    job_kwargs=balsam_exe_settings,
    **calc_keywords
)

atoms.calc.set(**extra_calc_keywords)

opt = QuasiNewton(atoms=atoms, trajectory=big_slab_name + '.traj')
# opt = Sella(atoms, order=0, delta0=1e-2, trajectory=jobdir + '.traj')
opt.run(fmax=0.01)
atoms.calc.close()

write_dir = os.path.join(big_slab_name + '.xyz')
write(write_dir, read(big_slab_name + '.traj'))
