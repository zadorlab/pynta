from pynta.excatkit.gratoms import Gratoms
from molecule.molecule import Molecule
from ase.io import read, write
from ase.data import covalent_radii
import numpy as np

def molecule_to_gratoms(mol):
    symbols = []
    edges = []
    surf_indexes = []
    surf_index_atom_map = {}
    atom_map = {}
    c = 0
    tags = []
    for i,atm in enumerate(mol.atoms):
        if not atm.is_surface_site():
            symbols.append(atm.element.symbol)
            tags.append(c)
            atom_map[c] = i
            if atm.is_bonded_to_surface():
                surf_indexes.append(c)
                surf_index_atom_map[c] = i
            c += 1

    for edge in mol.get_all_edges():
        if not edge.atom1.is_surface_site() and not edge.atom2.is_surface_site():
            ind1 = mol.atoms.index(edge.atom1)
            ind2 = mol.atoms.index(edge.atom2)
            edges.append((ind1,ind2))

    gra = Gratoms(symbols=symbols,edges=edges)
    gra.set_tags(tags)

    return gra,surf_indexes,atom_map,surf_index_atom_map

def get_edges(slab_path, find_surface=False):
    ''' Get adsorption edges

    Parameters
    ___________
    find_surface : bool
        specify whether to include surface or not
        default = False

    Returns
    ________
    edges : list[tuple]
        adsobrtion edges
    surface : numpy.ndarray
        an array with tags describing:
        top surface atoms (1)
        bottom surface (-1)
        and bulk atoms (0)
        Atoms with tags '1' are considered as the possible binding spots

    '''
    # read slab as an Atom object
    slab_atom = read(slab_path)

    # If the Atoms object is periodic, we need to check connectivity
    # across the unit cell boundary as well.
    tvecs = np.array([[0., 0., 0.]])
    if np.any(slab_atom.pbc):
        cell = slab_atom.cell
        # We are looking for atoms that are at most 2.5 times the
        # greatest covalent radius of all atoms in the system, so we
        # limit the number of periodic images we consider in each
        # direction to those that are within this cutoff radius of any
        # point on the face of the unit cell
        cutoff = max(covalent_radii[slab_atom.numbers]) * 2.5

        # If you want to understand the code below, read the `find_mic`
        # code in ASE, located at ase/geometry/geometry.py
        latt_len = np.sqrt((cell**2).sum(1))
        V = slab_atom.get_volume()
        padding = slab_atom.pbc * np.array(np.ceil(
            cutoff * np.prod(latt_len) / (V * latt_len)),
            dtype=int
        )
        offsets = np.mgrid[-padding[0]:padding[0] + 1,
                           -padding[1]:padding[1] + 1,
                           -padding[2]:padding[2] + 1].T
        tvecs = np.dot(offsets, cell).reshape(-1, 3)

    edges = []

    # pairvecs = np.ndarray(0)
    pairvecs = np.zeros_like(slab_atom.positions)

    # if find_surface:
    #     pairvecs = np.zeros_like(slab_atom.positions)
    # else:
    #     pairvecs = np.ndarray(0)

    for atomi in slab_atom:
        for atomj in slab_atom:
            i = atomi.index
            j = atomj.index
            # Like above, consider only bonds where j >= i. Note, we do
            # need to consider bonds where j == i because of periodic
            # boundary conditions.
            if j < i:
                continue
            # 1.25 times the sum of the covalent radii was chosen based
            # on trial and error. Too small, and you miss neighbors.
            # Too big, and you start including next-nearest-neighbors.
            cutoff = 1.25 * (
                covalent_radii[atomi.number] + covalent_radii[atomj.number]
            )
            # xij is the direct displacement vector in the central unit
            # cell.
            xij = atomj.position - atomi.position
            nbonds = 0
            # Loop over all neighboring unit cells
            for tvec in tvecs:
                # ...including the central unit cell. If i == j, then
                # explicitly skip the central unit cell.
                if i == j and np.all(tvec == [0., 0., 0.]):
                    continue
                # Count up the number of times i bonds to j via pbc
                if np.linalg.norm(xij + tvec) < cutoff:
                    if find_surface:
                        dx = xij + tvec
                        dx /= np.linalg.norm(dx)
                        pairvecs[i] -= dx
                        pairvecs[j] += dx
                    nbonds += 1
            # CatKit uses NetworkX "MultiGraph"s for periodic systems.
            # I'm not entirely sure how to interpret the edges for a
            # multigraph, but this is how CatKit wants multiple edges
            # between the same two atoms to be specified.
            for k in range(nbonds):
                edges.append((i, j, k))
    if not find_surface:
        return edges

    surface = np.zeros(len(slab_atom), dtype=int)

    for i, pairvec in enumerate(pairvecs):
        # Big pairvec means highly asymmetrical atoms, which
        # implies that it is near or at a surface
        if np.linalg.norm(pairvec) > 1.:
            # Use sign to determine whether it is at the top or
            # the bottom of the slab
            surface[i] = int(np.sign(pairvec[2]))
    return edges, surface

def get_grslab(slab_path):
    ''' Convert surface slab Atoms object into Gratoms object

    Returns
    -------
    grslab : Gratoms
        Gratoms representation of the surface slab - ready to place
        adsorbates

    '''
    slabedges, tags = get_edges(slab_path, True)
    slab_atom = read(slab_path)
    grslab = Gratoms(numbers=slab_atom.numbers,
                     positions=slab_atom.positions,
                     cell=slab_atom.cell,
                     pbc=slab_atom.pbc,
                     edges=slabedges)
    grslab.arrays['surface_atoms'] = tags

    return grslab
