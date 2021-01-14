# from .adjacency_to_3d import adjacency_to_3d
# from .compare_structures import find_all_unique
# from .relax_3d import create_relax_jobs
import os

"""Catalysis Generator."""

from collections import MutableMapping
import numpy as np
import ase

radicals = np.ones(92)
radicals[[6, 7, 8, 9, 15, 16]] = [4, 3, 2, 1, 3, 2]


class Defaults(MutableMapping, dict):
    """No frills default dictionary class."""

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

# __all__ = ['defaults', 'symmetry', 'adsorption', 'surface', 'molecules']


class Licence():
    def __init__(self):
        path = os.path.dirname(__file__)
        self.licence_info = os.path.join(path, 'license', 'license_info.py')

    def show_banner(self):
        with open(self.licence_info) as f:
            banner = f.read()
        print(banner)


Licence().show_banner()
