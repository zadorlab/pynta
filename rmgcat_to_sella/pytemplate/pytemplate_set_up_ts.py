#!/usr/bin/env python3
import os
import shutil

from ase.io import read, write

from sella import Sella

from ase.constraints import FixAtoms

import datetime

rxn = '{rxn}'
prefix = '{prefix}'
trajdir = os.path.join(prefix, prefix + '_' + rxn + '.traj')
label = os.path.join(prefix, prefix)

start = datetime.datetime.now()
with open(label + '_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()




TS_est = read('{TS}')
# fix all atoms but not adsorbates
# TS_est.set_constraint(FixAtoms([atom.index for atom in TS_est if atom.index < len(TS_est) - 2]))
# fix bottom half of the slab
TS_est.set_constraint(FixAtoms([atom.index for atom in TS_est if atom.position[2] < TS_est.cell[2, 2] / 2.]))

from pathlib import Path
cwd = Path.cwd().as_posix()
from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO
EspressoBalsamSocketIO.exe = executable
job_kwargs=balsam_exe_settings.copy()
job_kwargs.update([('user_workdir',cwd)])
QE_keywords_slab=calc_keywords.copy()
QE_keywords.update([('pseudopotentials',{pseudopotentials}),'pseudo_dir','{pseudo_dir}',('label',label)])
Calc = EspressoBalsamSocketIO(
    workflow='QE_Socket',
    job_kwargs=job_kwargs,
    **QE_keywords_slab
    )

TS_est.calc = Calc
opt = Sella(TS_est, order=1, delta0=1e-2, gamma=1e-16, trajectory = trajdir)
opt.run(fmax=0.01)
TS_est.calc.close()

WriteDir = os.path.join(prefix, prefix + '_' + rxn)
write(WriteDir + '_ts_final.png', read(trajdir))
write(WriteDir + '_ts_final.xyz', read(trajdir))

end = datetime.datetime.now()
with open(label + '_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()
