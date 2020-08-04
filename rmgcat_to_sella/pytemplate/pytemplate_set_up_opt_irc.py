#!/usr/bin/env python3
import os
import shutil

from ase.io import read, write

from sella import Sella

from ase.constraints import FixAtoms

import datetime

rxn = '{rxn}'
prefix = '{prefix}'
trajdir = os.path.join(prefix + '_' + rxn + '.traj')
# jobdir = os.path.join()
# label = os.path.join(prefix, prefix)

start = datetime.datetime.now()
with open(prefix + '_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()

# unixsocket = '_'.join([rxn, prefix])
# unixsocket = '{prefix}/{prefix}'.format(prefix=prefix)

from pathlib import Path
cwd = Path.cwd().as_posix()
from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO
EspressoBalsamSocketIO.exe = executable
job_kwargs=balsam_exe_settings.copy()
#job_kwargs.update([('user_workdir',cwd)])
QE_keywords=calc_keywords.copy()
QE_keywords.update([('pseudopotentials',{pseudopotentials}),'pseudo_dir','{pseudo_dir}',('label',prefix)])
Calc = EspressoBalsamSocketIO(
    workflow='QE_Socket',
    job_kwargs=job_kwargs,
    pseudopotentials=self.pseudopotentials,
    pseudo_dir=self.pseudo_dir,
    **QE_keywords
    )


geom_opt = read('{geom}')
geom_opt.set_constraint(FixAtoms(
    [atom.index for atom in geom_opt if atom.position[2] < geom_opt.cell[2, 2] / 2.]))

geom_opt.calc = Calc
opt = Sella(geom_opt, order=0, delta0=1e-2, trajectory=trajdir)
opt.run(fmax=0.01)
geom_opt.calc.close()

WriteDir = os.path.join(prefix + '_' + rxn + '_final')
write(WriteDir + '.png', read(trajdir))
write(WriteDir + '.xyz', read(trajdir))

end = datetime.datetime.now()
with open(prefix + '_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()
