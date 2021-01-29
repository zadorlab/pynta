# from .adjacency_to_3d import adjacency_to_3d
# from .compare_structures import find_all_unique
# from .relax_3d import create_relax_jobs
import os

"""Catalysis Generator."""

from collections import MutableMapping
import numpy as np
import ase

from pynta.__version__ import __version__


radicals = np.ones(92)
radicals[[6, 7, 8, 9, 15, 16]] = [4, 3, 2, 1, 3, 2]


class Defaults(MutableMapping, dict):
    '''No frills default dictionary class.'''

    def __init__(self):
        self.update({
            'radii': ase.data.covalent_radii.copy(),
            'radicals': radicals,
            'orthogonal': False
        })

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)

    def __getitem__(self, key):
        return dict.__getitem__(self, key)


defaults = Defaults()


class Licence():
    def __init__(self):
        self.path = os.path.dirname(__file__)
        self.licence_info = os.path.join(self.path, 'license', 'banner.txt')

    def show_banner(self) -> None:
        ''' Print a copyright banner everytime Pynta is imported'''
        path_to_excatkit = os.path.join(self.path, 'excatkit')
        with open(self.licence_info, 'r') as infile:
            banner = infile.read()
            with open(self.licence_info, 'w') as outfile:
                outfile.write(banner.format(
                    path_to_excatkit=path_to_excatkit))
        with open(self.licence_info, 'r') as infile:
            banner = infile.read()
        print(banner)


Licence().show_banner()
