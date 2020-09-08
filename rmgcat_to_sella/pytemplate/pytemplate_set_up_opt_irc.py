#!/usr/bin/env python3
import os
import datetime

from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO

from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton

rxn = '{rxn}'
prefix = '{prefix}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
traj_dir = os.path.join(prefix + '_' + rxn + '.traj')
label = os.path.join(prefix, prefix)

start = datetime.datetime.now()
with open(prefix + '_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()

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

opt = QuasiNewton(atoms=geom_opt, trajectory=traj_dir)
opt.run(fmax=0.01)
geom_opt.calc.close()

write_dir = os.path.join(prefix + '_' + rxn + '_final')
write(write_dir + '.png', read(traj_dir))
write(write_dir + '.xyz', read(traj_dir))

end = datetime.datetime.now()
with open(prefix + '_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()
