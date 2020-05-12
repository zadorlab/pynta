#!/usr/bin/env python3
#SBATCH -J {adsorbate}_{prefix}_relax
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p day-long-cpu
#SBATCH -t 1-00:00:00
#SBATCH -e %x.err
#SBATCH -o %x.out

import os
import shutil
import getpass

from ase.io import read, write
from ase.constraints import FixAtoms

from sella import Sella

from ase.calculators.espresso import Espresso
from ase.calculators.socketio import SocketIOCalculator

import datetime

adsorbate = '{adsorbate}'
prefix = '{prefix}'

user = getpass.getuser()

unixsocket = '{{adsorbate}}_{{prefix}}'.format(adsorbate=adsorbate, prefix=prefix)
socketpath = f'/tmp/ipi_{{unixsocket}}'
if os.path.exists(socketpath):
    os.remove(socketpath)

jobdir = os.path.join(adsorbate, prefix)
outdir = os.path.join(jobdir, prefix)

if os.path.exists(jobdir):
    shutil.rmtree(jobdir)
os.mkdir(jobdir)

label = os.path.join(jobdir, prefix)

start = datetime.datetime.now()

with open(outdir + '_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()

atoms = read(jobdir + '.xyz')
atoms.set_constraint(FixAtoms([atom.index for atom in atoms if atom.position[2] < atoms.cell[2, 2] / 2.]))

espresso = Espresso(command='/home/ehermes/local/bin/mpirun -np 48 /home/ehermes/local/bin/pw.x -inp PREFIX.pwi --ipi {{unixsocket}}:UNIX > PREFIX.pwo'
                            .format(unixsocket=unixsocket),
                    label=outdir,
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

with SocketIOCalculator(espresso, unixsocket=unixsocket) as calc:
    atoms.calc = calc
    opt = Sella(atoms, order=0, delta0=1e-2, trajectory=jobdir + '.traj')
    opt.run(fmax=0.01)

end = datetime.datetime.now()

with open(outdir + '_time.log', 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()

