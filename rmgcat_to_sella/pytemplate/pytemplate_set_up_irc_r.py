#!/usr/bin/env python3
import os
import datetime

from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO

from ase.io import read, write
from ase.constraints import FixAtoms
from sella import IRC

rxn = '{rxn}'
prefix = '{prefix}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
trajdir = os.path.join(prefix, prefix + '_' + rxn + '_irc_r.traj')
label = os.path.join(prefix, prefix)

start = datetime.datetime.now()
with open(label + '_irc_r_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()


extra_calc_keywords = dict(
    pseudopotentials={pseudopotentials},
    pseudo_dir='{pseudo_dir}',
    label=label
)

ts_geom = read('./{ts_xyz}')
ts_geom.set_constraint(FixAtoms([
    x.index for x in ts_geom if x.position[2] < ts_geom.cell[2, 2] / 2.
]))

ts_geom.calc = EspressoBalsamSocketIO(
    workflow='QE_Socket',
    job_kwargs=balsam_exe_settings,
    **calc_keywords
)

ts_geom.calc.set(**extra_calc_keywords)

opt = IRC(ts_geom, trajectory=trajdir, dx=0.1, eta=1e-4, gamma=1e-3)
opt.run(fmax=0.1, steps=1000, direction='reverse')

ts_geom.calc.close()

png_write_dir_r = os.path.join(prefix, prefix + '_' + rxn + '_irc_r.png')
write(png_write_dir_r, read(trajdir))

#####

end = datetime.datetime.now()
with open(label + '_irc_r_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()
