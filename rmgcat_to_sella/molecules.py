def molecule(species, bond_index=None, vacuum=0):
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
    molecule_graphs = catkit.gen.molecules.get_topologies(species)

    images = []
    for atoms in molecule_graphs:
        atoms = catkit.gen.molecules.get_3D_positions(atoms, bond_index)
        atoms.center(vacuum)
        images += [atoms]

    return images


def get_3D_positions(atoms, bond_index=None):
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
    branches = nx.dfs_successors(atoms.graph, bond_index)

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
