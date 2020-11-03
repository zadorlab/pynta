#!/usr/bin/env python3

import os
import shutil

import datetime
from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO

from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton

geom = '{geom}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'

start = datetime.datetime.now()

with open(outdir + '_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()

atoms = read(jobdir + '.xyz')
atoms.set_constraint(FixAtoms([
    atom.index for atom in atoms if atom.position[2] < atoms.cell[2, 2] / 2.
]))

extra_calc_keywords = dict(
    pseudo_dir='{pseudo_dir}',
    pseudopotentials={pseudopotentials},
    label=label
)


atoms.calc = EspressoBalsamSocketIO(
    workflow='QE_Socket',
    job_kwargs=balsam_exe_settings,
    **calc_keywords
)

atoms.calc.set(**extra_calc_keywords)

opt = QuasiNewton(atoms=atoms, trajectory=jobdir + '.traj')
# opt = Sella(atoms, order=0, delta0=1e-2, trajectory=jobdir + '.traj')
opt.run(fmax=0.06)
atoms.calc.close()

png_write_dir = os.path.join(jobdir + '_final.png')
write(png_write_dir, read(jobdir + '.traj'))

end = datetime.datetime.now()

with open(outdir + '_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()
