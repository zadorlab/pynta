#!/usr/bin/env python3
import numpy as np
import spglib


class Symmetry():
    ''' Wrapper for the ``spglib`` package. '''

    def __init__(self, atoms, tol=1e-5, ang_tol=-1):
        ''' Atoms object interface with spglib `symmetry finder
            <https://atztogo.github.io/spglib/python-spglib.html#python-spglib>`:

        Parameters
        ----------
        atoms : Atoms object
            Atomic structure to return the symmetry operations for.
        tol : float
            Tolerance for floating point precision errors.

        '''
        self.lattice = atoms.cell
        self.positions = atoms.get_scaled_positions()
        self.numbers = atoms.get_atomic_numbers()
        self.magmoms = atoms.get_initial_magnetic_moments()
        self.modified_numbers = get_modified_spin_symbols(
            self.numbers, self.magmoms)
        self.tol = tol

        cell = (self.lattice, self.positions, self.modified_numbers)
        self.data = spglib.get_symmetry_dataset(
            cell, symprec=tol, angle_tolerance=ang_tol)

    def get_symmetry_operations(
            self,
            affine=True):
        ''' Return the symmetry operations for a given atomic structure.

        Parameters
        ----------
        affine : bool
            Whether to return the affine matrix operations.

        Returns
        -------
        rotations : ndarray (N, 3, 3)
            Rotation matices of the symmetry operations.
        translations ndarray (N, 3)
            Translation vector components of the symmetry operations.
        affine_matrices ndarray (N, 4, 4)
            Affine matrix operations, combinations of the rotation and
            translation with ones along the diagonal.
        '''
        rotations = self.data['rotations'][1:]
        translations = self.data['translations'][1:]

        if affine:
            affine_matrices = np.zeros((rotations.shape[0], 4, 4))
            affine_matrices[:, :3, :3] = rotations
            affine_matrices[:, -1, :3] = translations
            affine_matrices[:, -1, -1] = 1
            return affine_matrices

        return rotations, translations

    def get_pointgroup(
            self,
            check_laue=False):
        ''' Return the point group operations of a systems.

        Parameters
        ----------
        check_laue: bool
            Return if the pointgroup is a laue symmetry.

        Returns
        -------
        pointgroup: str
            The pointgroup symmetry of the atomic structure.
        is_laue: bool
            Whether the pointgroup is a laue symmetry.

        '''
        pointgroup = self.data['pointgroup']

        if check_laue:
            laue = ['-1', '2/m', 'mmm', '4/m', '4/mmm',
                    '-3', '-3m', '6/m', '6/mmm', 'm-3', 'm-3m']
            is_laue = pointgroup in laue

            return pointgroup, is_laue

        return pointgroup

    def get_lattice_name(self):
        ''' Return the lattice name of an atoms object based
        on its `spacegroup number <https://en.wikipedia.org/wiki/List_of_space_groups/>`:

        Returns
        -------
        lattice: str
            The name of the structures lattice.

        '''
        space_group_number = self.data['number']

        if space_group_number in [146, 148, 155, 160, 161, 166, 167]:
            return 'rhombohedral'

        lattices = {
            'triclinic': 2,
            'monoclinic': 15,
            'orthorhombic': 74,
            'tetragonal': 142,
            'hexagonal': 194,
            'cubic': 230}

        for lattice, max_number in lattices.items():
            if space_group_number <= max_number:
                return lattice


def get_modified_spin_symbols(numbers, magmoms):
    '''Return a representation of atomic symbols which is
    unique to the magnetic moment as well.

    This is effectivly creating a single integer which contains the
    atomic number and the magnetic moment multiplied by 10.

    Parameters
    ----------
    numbers: ndarray(N,)
        Atomic numbers to be joined with the magnetic moments.
    magmoms: ndarray(N,)
        Magnetic moments to be joined to the atomic numbers.

    Returns
    -------
    spin_mod_symbols: ndarray(N,)
        The spin modified symbols representation for each atom.

    '''
    spin_mod_symbols = numbers.copy()
    magmoms = magmoms * 10
    magmoms = magmoms.astype(int)

    sign = np.sign(magmoms)
    spin_mod_symbols *= 1000
    spin_mod_symbols += np.abs(magmoms)
    ind = np.where(sign)
    spin_mod_symbols[ind] *= sign[ind]

    return spin_mod_symbols
