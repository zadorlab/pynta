#!/usr/bin/env python3
from pynta.excatkit.gratoms import Gratoms
from pynta import defaults

from ase.data import chemical_symbols
from networkx import dfs_successors, is_connected, from_numpy_matrix
from itertools import combinations
import numpy as np
import re


class Molecule():
    def molecule(
            self,
            species,
            bond_index=None,
            vacuum=0):
        '''Return list of enumerated gas-phase molecule structures based
        on species and topology.

        Parameters
        ----------
        species : str
            The chemical symbols to construct a molecule from.
        bond_index : int
            Construct the molecule as though it were adsorbed to a surface
            parallel to the z-axis. Will bond by the atom index given.
        vacuum : float
            Angstroms of vacuum to pad the molecules with.

        Returns
        -------
        images : list of Gratoms objects
            3D structures of the requested chemical species and topologies.

        '''
        molecule_graphs = self.get_topologies(species)
        images = []
        for atoms in molecule_graphs:
            atoms = Molecule.get_3D_positions(atoms, bond_index)
            atoms.center(vacuum)
            images += [atoms]

        return images

    @staticmethod
    def get_atomic_numbers(
            formula,
            return_count=False):
        ''' Return the atomic numbers associated with a chemical formula.

        Parameters
        ----------
        formula : string
            A chemical formula to parse into atomic numbers.
        return_count : bool
            Return the count of each element in the formula.

        Returns
        -------
        numbers : ndarray (n,)
            Element numbers in associated species.
        counts : ndarray (n,)
            Count of each element in a species.

        '''
        parse = re.findall('[A-Z][a-z]?|[0-9]+', formula)

        values = {}
        for i, e in enumerate(parse):
            if e.isdigit():
                values[parse[i - 1]] += int(e) - 1
            else:
                if e not in values:
                    values[e] = 1
                else:
                    values[e] += 1

        numbers = np.array([
            chemical_symbols.index(k) for k in values.keys()])
        srt = np.argsort(numbers)
        numbers = numbers[srt]

        if return_count:
            counts = np.array([v for v in values.values()])[srt]

            return numbers, counts

        return numbers

    @staticmethod
    def hydrogenate(
            atoms,
            bins,
            copy=True):
        ''' Add hydrogens to a gratoms object via provided bins

        '''
        h_index = len(atoms)

        edges = []
        for i, j in enumerate(bins):
            for _ in range(j):
                edges += [(i, h_index)]
                h_index += 1

        if copy:
            atoms = atoms.copy()
        atoms += Gratoms('H{}'.format(sum(bins)))
        atoms.graph.add_edges_from(edges)

        return atoms

    def bin_hydrogen(
            self,
            hydrogens=1,
            bins=1):
        ''' Recursive function for determining distributions of
            hydrogens across bins.

        '''
        if bins == 1:
            yield [hydrogens]

        elif hydrogens == 0:
            yield [0] * bins

        else:
            for i in range(hydrogens + 1):
                for j in self.bin_hydrogen(hydrogens - i, 1):
                    for k in self.bin_hydrogen(i, bins - 1):
                        yield j + k

    def get_topologies(
            self,
            symbols,
            saturate=False):
        ''' Return the possible topologies of a given chemical species.

        Parameters
        ----------
        symbols : str
            Atomic symbols to construct the topologies from.
        saturate : bool
            Saturate the molecule with hydrogen based on the
            default.radicals set.

        Returns
        -------
        molecules : list (N,)
            Gratoms objects with unique connectivity matrix attached.
            No 3D positions will be provided for these structures.

        '''
        num, cnt = Molecule.get_atomic_numbers(symbols, True)
        # print(num, cnt)
        mcnt = cnt[num != 1]
        mnum = num[num != 1]

        if cnt[num == 1]:
            hcnt = cnt[num == 1][0]
        else:
            hcnt = 0

        elements = np.repeat(mnum, mcnt)
        max_degree = defaults.get('radicals')[elements]
        n = mcnt.sum()

        hmax = int(max_degree.sum() - (n - 1) * 2)
        if hcnt > hmax:
            hcnt = hmax

        if saturate:
            hcnt = hmax
        if n == 1:
            atoms = Gratoms(elements, cell=[1, 1, 1])
            hatoms = Molecule.hydrogenate(atoms, np.array([hcnt]))
            return [hatoms]
        elif n == 0:
            hatoms = Gratoms('H{}'.format(hcnt))
            if hcnt == 2:
                hatoms.graph.add_edge(0, 1, bonds=1)
            return [hatoms]

        ln = np.arange(n).sum()
        il = np.tril_indices(n, -1)

        backbones, molecules = [], []
        combos = combinations(np.arange(ln), n - 1)
        for c in combos:
            # Construct the connectivity matrix
            ltm = np.zeros(ln)
            ltm[np.atleast_2d(c)] = 1

            connectivity = np.zeros((n, n))
            connectivity[il] = ltm
            connectivity = np.maximum(connectivity, connectivity.T)

            degree = connectivity.sum(axis=0)

            # Not fully connected (subgraph)
            if np.any(degree == 0) or not \
                    is_connected(from_numpy_matrix(connectivity)):
                continue

            # Overbonded atoms.
            remaining_bonds = (max_degree - degree).astype(int)
            if np.any(remaining_bonds < 0):
                continue

            atoms = Gratoms(
                numbers=elements,
                edges=connectivity,
                cell=[1, 1, 1])

            isomorph = False
            for G0 in backbones:
                if atoms.is_isomorph(G0):
                    isomorph = True
                    break

            if not isomorph:
                backbones += [atoms]

                # The backbone is saturated, do not enumerate
                if hcnt == hmax:
                    hatoms = Molecule.hydrogenate(atoms, remaining_bonds)
                    molecules += [hatoms]
                    continue

                # Enumerate hydrogens across backbone
                for bins in self.bin_hydrogen(hcnt, n):
                    if not np.all(bins <= remaining_bonds):
                        continue

                    hatoms = Molecule.hydrogenate(atoms, bins)

                    isomorph = False
                    for G0 in molecules:
                        if hatoms.is_isomorph(G0):
                            isomorph = True
                            break

                    if not isomorph:
                        molecules += [hatoms]

        return molecules

    @staticmethod
    def get_basis_vectors(coordinates):
        ''' Return a set of basis vectors for a given array of
            3D coordinates.

        Parameters
        ----------
        coordinates : array_like (3, 3) | (2, 3)
            Cartesian coordinates to determine the basis of. If
            only 2 positions are given 3rd is chosen as the positive
            y-axis.

        Returns
        -------
        basis_vectors : ndarray (3, 3)
            Automatically generated basis vectors from the given
            positions.

        '''
        if len(coordinates) == 3:
            c0, c1, c2 = coordinates
        else:
            c0, c1 = coordinates
            c2 = np.array([0, 1, 0])

        basis1 = c0 - c1
        basis2 = np.cross(basis1, c0 - c2)
        basis3 = np.cross(basis1, basis2)

        basis_vectors = np.vstack([basis1, basis2, basis3])
        basis_vectors /= np.linalg.norm(
            basis_vectors, axis=1, keepdims=True)

        return basis_vectors

    @staticmethod
    def branch_molecule(
            atoms,
            branch,
            basis=None,
            adsorption=None):
        ''' Return the positions of a Gratoms object for a segment of its
        attached graph. This function is mean to be iterated over by a
        depth first search form NetworkX.

        Parameters
        ----------
        atoms : Gratoms object
            Gratoms object to be iterated over. Will have its positions
            altered in-place.
        branch : tuple (1, [N,])
            A single entry from the output of nx.bfs_successors. The
            first entry is the root node and the following list is the
            nodes branching from the root.
        basis : ndarray (3, 3)
            The basis vectors to use for this step of the branching.
        adsorption : bool
            If True, will construct the molecule as though there is a
            surface normal to the negative z-axis. Must be None for
            all but the first index in the depth first search.

        Returns
        -------
        positions : ndarray (N, 3)
            Estimated positions for the branch nodes.

        '''
        root, nodes = branch
        root_position = atoms[root].position

        radii = defaults.get('radii')
        atomic_numbers = atoms.numbers[[root] + nodes]
        atomic_radii = radii[atomic_numbers]
        dist = (atomic_radii[0] + atomic_radii[1:])[:, None]

        angles = np.array([109.47, 109.47, 109.47, 0])
        dihedral = np.array([0, 120, -120, 0])

        if adsorption:
            # Move adsorption structures away from surface
            angles += 15

        # Tetrahedral bond arrangement by default
        n = len(nodes)
        if n == 1:
            # Linear bond arrangement
            angles[0] = 180
        elif n == 2:
            # Trigonal-planer bond arrangement
            angles[:2] = 120
            dihedral[1] = 180

        # Place the atoms of this segment of the branch
        if basis is None:
            basis = Molecule.get_basis_vectors(
                [root_position, [0, 0, -1]])
        basis = np.repeat(basis[None, :, :], len(dist), axis=0)

        ang = np.deg2rad(angles)[:len(dist), None]
        tor = np.deg2rad(dihedral)[:len(dist), None]

        basis[:, 1] *= -np.sin(tor)
        basis[:, 2] *= np.cos(tor)

        vectors = basis[:, 1] + basis[:, 2]
        vectors *= dist * np.sin(ang)
        basis[:, 0] *= dist * np.cos(ang)

        positions = vectors + root_position - basis[:, 0]

        return positions

    @staticmethod
    def get_3D_positions(
            atoms,
            bond_index=None):
        ''' Return an estimation of the 3D structure of a Gratoms object
        based on its graph.

        .. warning:: This function operates on the atoms object in-place.

        Parameters
        ----------
        atoms : Gratoms object
            Structure with connectivity matrix to provide a 3D structure.
        bond_index : int
            Index of the atoms to consider as the origin of a surface
            bonding site.

        Returns
        -------
        atoms : Gratoms object
            Structure with updated 3D positions.

        '''
        branches = dfs_successors(atoms.graph, bond_index)

        complete = []
        for i, branch in enumerate(branches.items()):
            root, nodes = branch

            if len(nodes) == 0:
                continue

            c0 = atoms[root].position
            if i == 0:
                basis = Molecule.get_basis_vectors([c0, [0, 0, -1]])
            else:
                bond_index = None
                for j, base_root in enumerate(complete):
                    if root in branches[base_root]:
                        c1 = atoms[base_root].position
                        # Flip the basis for every alternate step down
                        # the chain.
                        basis = Molecule.get_basis_vectors([c0, c1])
                        if (i - j) % 2 != 0:
                            basis[2] *= -1
                        break
            complete.insert(0, root)

            positions = Molecule.branch_molecule(
                atoms, branch, basis, bond_index)
            atoms.positions[nodes] = positions

        return atoms
