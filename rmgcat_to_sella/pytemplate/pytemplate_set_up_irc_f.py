#!/usr/bin/env python3
import os
import shutil

from ase.io import read, write

from sella import IRC

from ase.constraints import FixAtoms

import datetime

rxn = '{rxn}'
prefix = '{prefix}'
trajdir = os.path.join(prefix, prefix + '_' + rxn + '_irc_f.traj')
# jobdir = os.path.join()
label = os.path.join(prefix, prefix)

start = datetime.datetime.now()
with open(label + '_irc_f_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()


from pathlib import Path
cwd = Path.cwd().as_posix()
from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO
EspressoBalsamSocketIO.exe = executable
job_kwargs=balsam_exe_settings.copy()
#job_kwargs.update([('user_workdir',cwd)])
QE_keywords=calc_keywords.copy()
QE_keywords.update([('pseudopotentials',{pseudopotentials}),('pseudo_dir','{pseudo_dir}'),('label',label)])
Calc = EspressoBalsamSocketIO(
    workflow='QE_Socket',
    job_kwargs=job_kwargs,
    **QE_keywords
    )



TS_geom = read('./{TS_xyz}')
TS_geom.set_constraint(FixAtoms(
    [atom.index for atom in TS_geom if atom.position[2] < TS_geom.cell[2, 2] / 2.]))

TS_geom.calc = Calc
opt = IRC(TS_geom, trajectory=trajdir, dx=0.1, eta=1e-4, gamma=0.4)
opt.run(fmax=0.1, steps=1000, direction='forward')

TS_geom.calc.close()

pngWriteDir_f = os.path.join(prefix, prefix + '_' + rxn + '_irc_f.png')
write(pngWriteDir_f, read(trajdir))

#####

end = datetime.datetime.now()
with open(label + '_irc_f_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()
