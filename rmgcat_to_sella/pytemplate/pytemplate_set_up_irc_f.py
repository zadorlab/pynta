#!/usr/bin/env python3
import os
import datetime

from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO

from ase.io import read, write
from ase.constraints import FixAtoms
from sella import IRC


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


EspressoBalsamSocketIO.exe = executable
extra_calc_keywords = dict(
        pseudopotentials={pseudopotentials},
        pseudo_dir='{pseudo_dir}',
        label=label
        )

TS_geom = read('./{TS_xyz}')
TS_geom.set_constraint(FixAtoms([
    x.index for x in TS_geom if x.position[2] < TS_geom.cell[2, 2] / 2.
]))

TS_geom.calc = EspressoBalsamSocketIO(
        workflow='QE_Socket',
        job_kwargs=balsam_exe_settings,
        **calc_keywords
        )

TS_geom.calc.set(**extra_calc_keywords)

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
