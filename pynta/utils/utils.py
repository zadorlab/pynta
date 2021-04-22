#!/usr/bin/env python

# by default, each atom is bonded to the surface through the atom
# at index 0. This can be overwrite by using this dict
import os
from ase.io import Trajectory, write
edge_cases_bonded_dict = dict(CO3=1, CH3OH=1, CH3O2=2, HCO=1)

# by default, the first topology structure is used to generate
# adsorbates. This can be modified using this dict
edge_cases_topology_dict = dict(COOH=1, HCOOH=1, CH3O2=1, HCOOCH3=17)


def convert_traj_to_xyz(path_to_traj):

    traj = Trajectory(path_to_traj)

    with open('traj.xyz', 'a+') as main:
        for i, step in enumerate(traj):
            i = str(i).zfill(2)
            f_name = os.path.join(os.path.dirname(
                path_to_traj), 'step_{}.xyz'.format(i))
            write(f_name, step)
            with open(f_name, 'r') as inner:
                inner_file = inner.read()
            main.write(inner_file)
            os.remove(f_name)
