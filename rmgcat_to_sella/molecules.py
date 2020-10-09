from rmgcat_to_sella.gratoms import Gratoms
from ase.data import chemical_symbols
from networkx import dfs_successors
import numpy as np
import re


class Molecule():
    def molecule(self, species, bond_index=None, vacuum=0):
        """Return list of enumerated gas-phase molecule structures based
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
        """
        molecule_graphs = self.get_topologies(species)

        images = []
        for atoms in molecule_graphs:
            atoms = self.get_3D_positions(atoms, bond_index)
            atoms.center(vacuum)
            images += [atoms]

        return images

    def get_atomic_numbers(self, formula, return_count=False):
        """Return the atomic numbers associated with a chemical formula.

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
        """
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

    def hydrogenate(self, atoms, bins, copy=True):
        """Add hydrogens to a gratoms object via provided bins"""
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

    def get_topologies(self, symbols, saturate=False):
        """Return the possible topologies of a given chemical species.

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
        """
        num, cnt = self.get_atomic_numbers(symbols, True)
        mcnt = cnt[num != 1]
        mnum = num[num != 1]

        if cnt[num == 1]:
            hcnt = cnt[num == 1][0]
        else:
            hcnt = 0

        elements = np.repeat(mnum, mcnt)
        max_degree = catkit.gen.defaults.get('radicals')[elements]
        n = mcnt.sum()

        hmax = int(max_degree.sum() - (n - 1) * 2)
        if hcnt > hmax:
            hcnt = hmax

        if saturate:
            hcnt = hmax

        if n == 1:
            atoms = Gratoms(elements, cell=[1, 1, 1])
            hatoms = hydrogenate(atoms, np.array([hcnt]))
            return [hatoms]
        elif n == 0:
            hatoms = catkit.Gratoms('H{}'.format(hcnt))
            if hcnt == 2:
                hatoms.graph.add_edge(0, 1, bonds=1)
            return [hatoms]

        ln = np.arange(n).sum()
        il = np.tril_indices(n, -1)

        backbones, molecules = [], []
        combos = itertools.combinations(np.arange(ln), n - 1)
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
                    nx.is_connected(nx.from_numpy_matrix(connectivity)):
                continue

            # Overbonded atoms.
            remaining_bonds = (max_degree - degree).astype(int)
            if np.any(remaining_bonds < 0):
                continue

            atoms = catkit.Gratoms(
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
                    hatoms = hydrogenate(atoms, remaining_bonds)
                    molecules += [hatoms]
                    continue

                # Enumerate hydrogens across backbone
                for bins in bin_hydrogen(hcnt, n):
                    if not np.all(bins <= remaining_bonds):
                        continue

                    hatoms = hydrogenate(atoms, bins)

                    isomorph = False
                    for G0 in molecules:
                        if hatoms.is_isomorph(G0):
                            isomorph = True
                            break

                    if not isomorph:
                        molecules += [hatoms]

        return molecules

    def get_3D_positions(self, atoms, bond_index=None):
        """Return an estimation of the 3D structure of a Gratoms object
        based on its graph.

        WARNING: This function operates on the atoms object in-place.

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
        """
        branches = dfs_successors(atoms.graph, bond_index)

        complete = []
        for i, branch in enumerate(branches.items()):
            root, nodes = branch

            if len(nodes) == 0:
                continue

            c0 = atoms[root].position
            if i == 0:
                basis = catkit.gen.utils.get_basis_vectors([c0, [0, 0, -1]])
            else:
                bond_index = None
                for j, base_root in enumerate(complete):
                    if root in branches[base_root]:
                        c1 = atoms[base_root].position
                        # Flip the basis for every alternate step down the chain.
                        basis = catkit.gen.utils.get_basis_vectors([c0, c1])
                        if (i - j) % 2 != 0:
                            basis[2] *= -1
                        break
            complete.insert(0, root)

            positions = _branch_molecule(atoms, branch, basis, bond_index)
            atoms.positions[nodes] = positions

        return atoms
