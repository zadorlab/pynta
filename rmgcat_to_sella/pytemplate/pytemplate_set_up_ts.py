#!/usr/bin/env python3
#SBATCH -J {prefix}_{rxn}_ts
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p week-long-cpu
#SBATCH -t 1-00:00:00
#SBATCH -e %x.err
#SBATCH -o %x.out

import os
import shutil

from ase.io import read, write

from sella import Sella

from ase.calculators.espresso import Espresso
from ase.calculators.socketio import SocketIOCalculator
from ase.constraints import FixAtoms

import datetime

rxn = '{rxn}'
prefix = '{prefix}'
trajdir = os.path.join(prefix, prefix + '_' + rxn + '.traj')
# jobdir = os.path.join()
label = os.path.join(prefix, prefix)

start = datetime.datetime.now()
with open(label + '_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()

# unixsocket = '_'.join([rxn, prefix])
# unixsocket = '{prefix}/{prefix}'.format(prefix=prefix)
unixsocket = '{{prefix}}'.format(prefix=prefix)
socketpath = f'/tmp/ipi_{{unixsocket}}'
if os.path.exists(socketpath):
    os.remove(socketpath)


espresso = Espresso(command='/home/ehermes/local/bin/mpirun -np 48 /home/ehermes/local/bin/pw.x -inp PREFIX.pwi --ipi {{unixsocket}}:UNIX > PREFIX.pwo'
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

TS_est = read('{TS}')
TS_est.set_constraint(FixAtoms([atom.index for atom in TS_est if atom.position[2] < TS_est.cell[2, 2] / 2.]))

with SocketIOCalculator(espresso, unixsocket=unixsocket) as calc:
    TS_est.calc = calc
    opt = Sella(TS_est, order=1, delta0=1e-2, trajectory = trajdir)
    opt.run(fmax=0.01)

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
