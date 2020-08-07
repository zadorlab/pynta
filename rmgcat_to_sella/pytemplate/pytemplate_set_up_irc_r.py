#!/usr/bin/env python3
import os
import shutil

from ase.io import read, write

from sella import IRC

from ase.constraints import FixAtoms

import datetime

rxn = '{rxn}'
prefix = '{prefix}'
trajdir = os.path.join(prefix, prefix + '_' + rxn + '_irc_r.traj')
# jobdir = os.path.join()
label = os.path.join(prefix, prefix)

start = datetime.datetime.now()
with open(label + '_irc_r_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()

# unixsocket = '_'.join([rxn, prefix])
# unixsocket = '{prefix}/{prefix}'.format(prefix=prefix)

from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO
EspressoBalsamSocketIO.exe = executable
extra_calc_keywords = dict(
        pseudopotentials={pseudopotentials},
        pseudo_dir='{pseudo_dir}',
        label=label
        )

TS_geom = read('./{TS_xyz}')
TS_geom.set_constraint(FixAtoms(
    [atom.index for atom in TS_geom if atom.position[2] < TS_geom.cell[2, 2] / 2.]))

TS_geom.calc = EspressoBalsamSocketIO(
        workflow='QE_Socket',
        job_kwargs=balsam_exe_settings,
        **calc_keywords
        )

TS_geom.calc.set(**extra_calc_keywords)

opt = IRC(TS_geom, trajectory=trajdir, dx=0.1, eta=1e-4, gamma=0.4)
opt.run(fmax=0.1, steps=1000, direction='reverse')
TS_geom.calc.close()

pngWriteDir_f = os.path.join(prefix, prefix + '_' + rxn + '_irc_r.png')
write(pngWriteDir_f, read(trajdir))

#####

end = datetime.datetime.now()
with open(label + '_irc_r_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()
