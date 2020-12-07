#!/usr/bin/env python3

from balsam.launcher.dag import BalsamJob, add_dependency
import os
import shutil

import datetime
from pynta.balsamcalc import EspressoBalsamSocketIO

from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton

adsorbate = '{adsorbate}'
prefix = '{prefix}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
jobdir = os.path.join(adsorbate, prefix)
outdir = os.path.join(jobdir, prefix)

if os.path.exists(jobdir):
    shutil.rmtree(jobdir)
os.mkdir(jobdir)
label = os.path.join(jobdir, prefix)
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

# kpts={repeats},
# jobs_args='-nk {n_kpts}',

# # update balsam_exe_settings with info about a new num_nodes
# balsam_exe_settings['num_nodes'] = {n_kpts}

atoms.calc = EspressoBalsamSocketIO(
    workflow='QE_Socket',
    job_kwargs=balsam_exe_settings,
    **calc_keywords
)

atoms.calc.set(**extra_calc_keywords)

opt = QuasiNewton(atoms=atoms, trajectory=jobdir + '.traj')
opt.run(fmax=0.01, steps=70)
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
