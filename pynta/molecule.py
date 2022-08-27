from pynta.excatkit.gratoms import Gratoms
from molecule.molecule import Molecule
from ase.io import read, write
from ase.data import covalent_radii
from acat.adsorption_sites import SlabAdsorptionSites
from acat.adsorbate_coverage import SlabAdsorbateCoverage
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

def get_labeled_bonds(mol):
    labeled_atoms = mol.get_all_labeled_atoms()

    bonds = dict()
    for label,atm in labeled_atoms.items():
        bonds[label] = []
        bds = mol.get_bonds(atm)
        for bd in bds.values():
            if bd.atom1 is atm:
                if bd.atom2.label:
                    bonds[label].append(bd.atom2.label)
            else:
                if bd.atom1.label:
                    bonds[label].append(bd.atom1.label)

    return bonds

def get_broken_formed_bonds(reactant,product):
    reactant_labeled_bonds = get_labeled_bonds(reactant)
    product_labeled_bonds = get_labeled_bonds(product)
    broken_bonds = set()
    formed_bonds = set()
    for label,rlabels in reactant_labeled_bonds.items():
        if label in product_labeled_bonds.keys():
            plabels = product_labeled_bonds[label]
            missing_from_p = [label for label in rlabels if label not in plabels]
            missing_from_r = [label for label in plabels if label not in rlabels]
            for label2 in missing_from_p:
                broken_bonds.add(frozenset((label,label2)))
            for label2 in missing_from_r:
                formed_bonds.add(frozenset((label,label2)))

    return broken_bonds,formed_bonds

def get_template_mol_map(template,mols):
    """
    template is the combined labeled Molecule object
    find dictionary mapping indices of the template to the molecule objects of interest
    output is a list with each item coresponding to a Molecule object in mols the dictionary
    maps indices of the template to indices of the corresponding mol object
    """
    tempmol_mol_map = []
    tempmols = [x for x in template.split() if not x.is_surface_site()]
    ordered_tempmols = []
    for i,mol in enumerate(mols):
        for j,tempmol in enumerate(tempmols):
            if tempmol.is_isomorphic(mol,save_order=True):
                mapv = tempmol.find_isomorphism(mol,save_order=True)
                indmap = {tempmol.atoms.index(key): mol.atoms.index(val) for key,val in mapv[0].items()}
                tempmol_mol_map.append(indmap)
                ordered_tempmols.append(tempmol)
                tempmols.pop(j)
                break
        else:
            print(mol.to_adjacency_list())
            for tempmol in tempmols:
                print(tempmol.to_adjacency_list())
            raise ValueError("mapping could not be found")

    temp_tempmol_map = []
    for i,tempmol in enumerate(ordered_tempmols):
        grptempmol = tempmol.to_group()
        grptempmol.multiplicity = [template.multiplicity]
        mapv = template.find_subgraph_isomorphisms(grptempmol,save_order=True)
        maps = get_nonintersectingkeys_maps(mapv)
        c = 0
        m = maps[c]
        indmap = {template.atoms.index(key): grptempmol.atoms.index(val) for key,val in m.items()}

        while set(list(indmap.keys())) in [set(list(v.keys())) for v in temp_tempmol_map]:
            c += 1
            m = maps[c]
            indmap = {template.atoms.index(key): grptempmol.atoms.index(val) for key,val in m.items()}
        temp_tempmol_map.append(indmap)

    temp_mol_map = []
    for i in range(len(tempmol_mol_map)):
        d = {tind: tempmol_mol_map[i][tmolind] for tind,tmolind in temp_tempmol_map[i].items()}
        temp_mol_map.append(d)
    return temp_mol_map

def get_nonintersectingkeys_maps(maps):
    mvals_init = [frozenset(list(m.keys())) for m in maps]
    mvals_unique = list(set(mvals_init))
    valid_map_inds = []
    for mval in mvals_unique:
        for inds in valid_map_inds:
            if mval.intersection(inds) != frozenset():
                break
        else:
            valid_map_inds.append(mval)

    inds = [mvals_init.index(mval) for mval in valid_map_inds]

    return [maps[ind] for ind in inds]

def ads_size(mol):
    return len(mol.atoms) - len(mol.get_surface_sites())

def get_mol_index(ind,template_mol_map):
    """
    ind is the index in the template
    first index corresponds to the molecule object
    second index corresponds to the index in the molecule object
    """
    for i,d in enumerate(template_mol_map):
        if ind in d.keys():
            return (i,d[ind])
    else:
        return None #surface_site

def get_ase_index(ind,template_mol_map,molecule_to_gratom_maps,nslab,ads_sizes):
    n = nslab
    for i,d in enumerate(template_mol_map):
        if ind in d.keys() and d[ind] in molecule_to_gratom_maps[i].keys():
            return molecule_to_gratom_maps[i][d[ind]]+n
        n += ads_sizes[i]
    else:
        return None #surface site

def get_bond_lengths_sites(mol,ads,atom_map,surf_atom_map,nslab,facet="fcc111",metal="Cu",cas=None):
    """
    gets bond lengths and site information indexed to the Molecule object atoms
    bond lengths is a matrix with the distances between all atoms that share bonds
    site information is the atom index mapped to the site type
    """
    rev_atom_map = {value: key for key,value in atom_map.items()} #map from mol to ads
    rev_surf_atom_map = { value: key for key,value in surf_atom_map.items()}

    if cas is None:
        cas = SlabAdsorptionSites(ads,facet,allow_6fold=False,composition_effect=False,
                            label_sites=True,
                            surrogate_metal=metal)
    adcov = SlabAdsorbateCoverage(ads,adsorption_sites=cas)
    occ = adcov.get_sites(occupied_only=True)
    surface_dict = [{"atom_index":x["bonding_index"]-nslab, "site":x["site"],
                       "bond_length":x["bond_length"], "position":x["position"]} for x in occ]
    ad = ads[nslab:]
    bondlengths = np.zeros((len(mol.atoms),len(mol.atoms)))
    sites = dict()
    sitelengths = dict()

    for bond in mol.get_all_edges():
        if bond.atom1.is_surface_site():
            ind1 = mol.atoms.index(bond.atom1)
            ind2 = mol.atoms.index(bond.atom2)
            aind = rev_surf_atom_map[ind2]
            surfd = [s for s in surface_dict if s["atom_index"] == aind][0]
            bondlengths[ind1,ind2] = surfd["bond_length"]
            bondlengths[ind2,ind1] = surfd["bond_length"]
            sites[ind1] = surfd["site"]
            sitelengths[ind1] = surfd["bond_length"]
        elif bond.atom2.is_surface_site():
            ind1 = mol.atoms.index(bond.atom1)
            ind2 = mol.atoms.index(bond.atom2)
            aind = rev_surf_atom_map[ind1]
            surfd = [s for s in surface_dict if s["atom_index"] == aind][0]
            bondlengths[ind1,ind2] = surfd["bond_length"]
            bondlengths[ind2,ind1] = surfd["bond_length"]
            sites[ind2] = surfd["site"]
            sitelengths[ind2] = surfd["bond_length"]
        else: #not a surface bond
            ind1 = mol.atoms.index(bond.atom1)
            ind2 = mol.atoms.index(bond.atom2)
            aind1 = rev_atom_map[ind1]
            aind2 = rev_atom_map[ind2]
            d = ad.get_distance(aind1,aind2)
            bondlengths[ind1,ind2] = d
            bondlengths[ind2,ind1] = d

    return bondlengths,sites,sitelengths
