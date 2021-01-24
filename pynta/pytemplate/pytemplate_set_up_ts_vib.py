#!/usr/bin/env python3

import os

import datetime

from ase.constraints import FixAtoms
from ase.vibrations import Vibrations
from ase.io import read

geom = '{geom}'
prefix = geom[:2]
socket_calculator = '{socket_calculator}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
nimages = {nimages}
n = {n}
vib_files_loc = os.path.join(os.getcwd(), prefix, 'vib')
geom_prefix = os.path.join(prefix, geom[:-10])

start = datetime.datetime.now()

with open(geom_prefix + '_time.log', 'w+') as f:
    f.write(str(start))
    f.write("\n")
    f.close()

atoms = read(os.path.join(prefix, geom))

# freeze all surface atoms
# atoms.set_constraint(FixAtoms(
#     [atom.index for atom in atoms if atom.symbol == 'Cu']))
# freeze half botom of the slab
atoms.set_constraint(FixAtoms([
    atom.index for atom in atoms if atom.position[2] < atoms.cell[2, 2] / 2.
]))
# vibrate only adsorbed species
# indices = [atom.index for atom in atoms if atom.symbol != 'Cu']
# vibrate adsorbates and 2 first layesr of the slab
indices = [atom.index for atom in atoms if atom.position[2]
           > atoms.cell[2, 2] / 2.]

# update balsam_exe_settings with info about a new num_nodes
# balsam_exe_settings['num_nodes'] = {n_kpts}

extra_calc_keywords = dict(
    pseudopotentials={pseudopotentials},
    pseudo_dir='{pseudo_dir}',
    label=geom
)

# kpts={repeats},
# jobs_args='-nk {n_kpts}',

balsamcalc_module = __import__('pynta.balsamcalc', fromlist=[
    socket_calculator])
sock_calc = getattr(balsamcalc_module, socket_calculator)

atoms.calc = sock_calc(
    workflow='QE_Socket',
    job_kwargs=balsam_exe_settings,
    **calc_keywords
)

atoms.calc.set(**extra_calc_keywords)

# start vibrations calculations
vib = Vibrations(atoms, indices=indices, name=vib_files_loc)
vib.run()
vib.summary()
# vib.clean()

# write the first vibration mode to vib.0.traj file (default) - imaginary freq
vib.write_mode(n=n, nimages=nimages)

end = datetime.datetime.now()

with open(os.path.join(prefix, geom[:-10] + '_time.log'), 'a+') as f:
    f.write(str(end))
    f.write("\n")
    f.write(str(end - start))
    f.write("\n")
    f.close()
