#!/usr/bin/env python3
#SBATCH -J {prefix}_{rxn}_irc_f
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p week-long-cpu
#SBATCH -t 1-00:00:00
#SBATCH -e %x.err
#SBATCH -o %x.out

import os
import shutil

from ase.io import read, write

from sella import IRC

from ase.calculators.espresso import Espresso
from ase.calculators.socketio import SocketIOCalculator
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

# unixsocket = '_'.join([rxn, prefix])
# unixsocket = '{prefix}/{prefix}'.format(prefix=prefix)
unixsocket = '{{prefix}}'.format(prefix=prefix)
socketpath = f'/tmp/ipi_{{unixsocket}}'
if os.path.exists(socketpath):
    os.remove(socketpath)


espresso = Espresso(command='/home/ehermes/local/bin/mpirun -np 48 /home/ehermes/local/bin/pw.x -inp PREFIX.pwi --ipi {{unixsocket}}:UNIX > irc_f.pwo'
                            .format(unixsocket=unixsocket),
                    label=label,
                    pseudopotentials=dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                                          H='H.pbe-kjpaw_psl.1.0.0.UPF',
                                          O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                                          C='C.pbe-n-kjpaw_psl.1.0.0.UPF',
                                          ),
                    pseudo_dir='/home/mgierad/espresso/pseudo',
                    kpts=(3, 3, 1),
                    occupations='smearing',
                    smearing='marzari-vanderbilt',
                    degauss=0.01,  # Rydberg
                    ecutwfc=40,  # Rydberg
                    nosym=True,  # Allow symmetry breaking during optimization
                    conv_thr=1e-11,
                    mixing_mode='local-TF',
                    )

TS_geom = read('./{TS_xyz}')
TS_geom.set_constraint(FixAtoms(
    [atom.index for atom in TS_geom if atom.position[2] < TS_geom.cell[2, 2] / 2.]))

with SocketIOCalculator(espresso, unixsocket=unixsocket) as calc:
    TS_geom.calc = calc
    opt = IRC(TS_geom, trajectory=trajdir, dx=0.1, eta=1e-4, gamma=0.4)
    opt.run(fmax=0.1, steps=1000, direction='forward')

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
