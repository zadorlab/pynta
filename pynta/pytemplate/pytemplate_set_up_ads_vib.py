#!/usr/bin/env python3
import os
import datetime

from ase.io import read
from ase.constraints import FixAtoms
from ase.vibrations import Vibrations

geom = '{geom}'
nimages = {nimages}
n = {n}
socket_calculator = '{socket_calculator}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
outdir = os.path.dirname(geom)
vib_files_loc = os.path.join(outdir, 'vib')
timelog_file = os.path.join(outdir, '00_time.log')

atoms = read(geom)

# freeze half botom of the slab
atoms.set_constraint(FixAtoms([
    atom.index for atom in atoms if atom.position[2] < atoms.cell[2, 2] / 2.
]))
# freeze all surface atoms
# atoms.set_constraint(FixAtoms(
#     [atom.index for atom in atoms if atom.symbol == 'Cu']))

# get indices of adsorbed species only
# indices_ads = [atom.index for atom in atoms if atom.symbol != 'Cu']

# get indices of 2 first layers of the slab
indices = [atom.index for atom in atoms if atom.position[2]
           > atoms.cell[2, 2] / 2.]

extra_calc_keywords = dict(
    pseudopotentials={pseudopotentials},
    pseudo_dir='{pseudo_dir}',
    label=geom
)

balsamcalc_module = __import__('pynta.balsamcalc', fromlist=[
    socket_calculator])
sock_calc = getattr(balsamcalc_module, socket_calculator)

atoms.calc = sock_calc(
    workflow='QE_Socket',
    job_kwargs=balsam_exe_settings,
    **calc_keywords
)

atoms.calc.set(**extra_calc_keywords)

start = datetime.datetime.now()
with open(timelog_file, 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()

# start vibrations calculations
vib = Vibrations(atoms, indices=indices, name=vib_files_loc)
vib.run()
vib.summary()

# write the first vibration mode to vib.0.traj file (default) - imaginary freq
vib.write_mode(n=n, nimages=nimages)

end = datetime.datetime.now()

with open(timelog_file, 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()
