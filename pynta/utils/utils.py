#!/usr/bin/env python
import os
from ase.io import (Trajectory, write, read)


# by default, each atom is bonded to the surface through the atom
# at index 0. This can be overwrite by using this dict
edge_cases_bonded_dict = dict(CO3=1, CH3OH=1, CH3O2=2, HCO=1)

# by default, the first topology structure is used to generate
# adsorbates. This can be modified using this dict
edge_cases_topology_dict = dict(COOH=1, HCOOH=1, CH3O2=1, HCOOCH3=17)


def convert_traj_to_xyz(path_to_traj: str) -> None:
    ''' Convert all geometries in path_to_traj file to xyz format and combine
    everything into one file

    Parameters
    ----------
    path_to_traj : str
        a path to *.traj file

    '''

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


def traj_to_xyz(path_to_traj: str) -> None:
    ''' Convert last geometry of path_to_traj to xyz file

    Parameters
    ----------
    path_to_traj : str
        a path to *.traj file

    '''
    dir_name = os.path.dirname(path_to_traj)
    f_name = os.path.join(dir_name, 'last_geom.xyz')
    write(f_name, read(path_to_traj))
