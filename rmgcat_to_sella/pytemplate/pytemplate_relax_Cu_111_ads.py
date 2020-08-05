#!/usr/bin/env python3

import os
import shutil

from ase.io import read, write
from ase.constraints import FixAtoms

from sella import Sella


import datetime

adsorbate = '{adsorbate}'
prefix = '{prefix}'
executable='{executable}'
balsam_exe_settings={balsam_exe_settings}
calc_keywords={calc_keywords}
creation_dir='{creation_dir}'

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
atoms.set_constraint(FixAtoms([atom.index for atom in atoms if atom.position[2] < atoms.cell[2, 2] / 2.]))

from pathlib import Path
cwd = Path.cwd().as_posix()
from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO
EspressoBalsamSocketIO.exe = executable
job_kwargs=balsam_exe_settings.copy()
#job_kwargs.update([('user_workdir',cwd)])
QE_keywords=calc_keywords.copy()
QE_keywords.update([('pseudopotentials',{pseudopotentials}),('pseudo_dir','{pseudo_dir}'),('label',outdir)])
Calc = EspressoBalsamSocketIO(
    workflow='QE_Socket',
    job_kwargs=job_kwargs,
    **QE_keywords
    )

atoms.calc = Calc
from ase.optimize import BFGSLineSearch
opt = BFGSLineSearch(atoms=atoms,trajectory=jobdir+'.traj')
#opt = Sella(atoms, order=0, delta0=1e-2, trajectory=jobdir + '.traj')
opt.run(fmax=0.01)
atoms.calc.close()

pngWriteFile = os.path.join(jobdir + '_final.png')
write(pngWriteFile, read(jobdir + '.traj'))

end = datetime.datetime.now()

with open(outdir + '_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()

