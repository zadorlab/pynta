#!/usr/bin/env python3
import os
import datetime
from pathlib import Path

from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO

from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import BFGSLineSearch

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

cwd = Path.cwd().as_posix()
extra_calc_keywords = dict(
        pseudopotentials={pseudopotentials},
        pseudo_dir='{pseudo_dir}',
        label=prefix
        )

geom_opt = read('{geom}')
geom_opt.set_constraint(FixAtoms([
    x.index for x in geom_opt if x.position[2] < geom_opt.cell[2, 2] / 2.
]))

geom_opt.calc = EspressoBalsamSocketIO(
        workflow='QE_Socket',
        job_kwargs=balsam_exe_settings,
        **calc_keywords
        )

geom_opt.calc.set(**extra_calc_keywords)

opt = BFGSLineSearch(atoms=geom_opt, trajectory=trajdir)
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
