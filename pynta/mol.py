from molecule.molecule import Molecule,Atom,Bond,GroupBond,ATOMTYPES
from ase.io import read, write
from ase.data import covalent_radii
from ase import Atoms
from ase.geometry import get_distances
from acat.adsorption_sites import SlabAdsorptionSites
from acat.adsorbate_coverage import SlabAdsorbateCoverage
from acat.settings import site_heights
from acat.utilities import get_mic
from acat.utilities import (custom_warning,
                         is_list_or_tuple,
                         get_close_atoms,
                         get_rodrigues_rotation_matrix,
                         get_angle_between,
                         get_rejection_between)
from pynta.utils import get_unique_sym_struct_index_clusters, get_unique_sym, get_unique_sym_structs, get_unique_sym_struct_indices, get_occupied_sites, sites_match
from pynta.calculator import run_harmonically_forced
from rdkit import Chem
from copy import deepcopy
import numpy as np
import random
import itertools
import logging

def get_desorbed_with_map(mol):
    molcopy = mol.copy(deep=True)
    init_map = {i:a for i,a in enumerate(molcopy.atoms)}

    for bd in molcopy.get_all_edges():
        if bd.atom1.is_surface_site():
            bd.atom2.radical_electrons += round(bd.order)
            molcopy.remove_bond(bd)
            molcopy.remove_atom(bd.atom1)
        elif bd.atom2.is_surface_site():
            bd.atom1.radical_electrons += round(bd.order)
            molcopy.remove_bond(bd)
            molcopy.remove_atom(bd.atom2)

    molcopy.sort_atoms()
    out_map = {i:molcopy.atoms.index(a) for i,a in init_map.items() if a in molcopy.atoms}
    return molcopy,out_map

def get_conformer(desorbed):
    try:
        rdmol,rdmap = desorbed.to_rdkit_mol(remove_h=False,return_mapping=True)
    except Exception as e:
        syms = [a.symbol for a in desorbed.atoms]
        indmap = {i:i for i in range(len(desorbed.atoms))}
        if len(desorbed.atoms) == 1:
            atoms = Atoms(syms[0],positions=[(0,0,0)])
            return atoms,indmap
        elif len(desorbed.atoms) == 2:
            atoms = Atoms(syms[0]+syms[1],positions=[(0,0,0),(1.3,0,0)])
            return atoms,indmap
        else:
            raise e

    indmap = {i:rdmap[a] for i,a in enumerate(desorbed.atoms)}
    Chem.AllChem.EmbedMultipleConfs(rdmol,numConfs=1,randomSeed=1)
    conf = rdmol.GetConformer()
    pos = conf.GetPositions()
    syms = [a.GetSymbol() for a in rdmol.GetAtoms()]
    atoms = Atoms(symbols=syms,positions=pos)
    return atoms,indmap

def get_adsorbate(mol):
    desorbed,mol_to_desorbed_map = get_desorbed_with_map(mol)
    atoms,desorbed_to_atoms_map = get_conformer(desorbed)
    mol_to_atoms_map = {key:desorbed_to_atoms_map[val] for key,val in mol_to_desorbed_map.items()}
    return atoms,mol_to_atoms_map


site_bond_length_dict = {
        ("ontop",None,None): 1.826370311,
        ("bridge",None,None): 1.806089179,
        ("fcc",None,None): 1.372599861,
        ("hcp",None,None): 1.397379832,

        ("ontop","C",None): 2.056904761,
        ("bridge","C",None): 1.920118777,
        ("fcc","C",None): 1.795024649,
        ("hcp","C",None): 1.764312871,
        ("ontop","O",None): 1.872603673,
        ("bridge","O",None): 1.69205958,
        ("fcc","O",None): 1.408365497,
        ("hcp","O",None): 1.510464567,
        ("ontop","H",None): 1.5496025,
        ("fcc","H",None): 1.0321708,
        ("hcp","H",None): 1.0321708,
        ("fcc","N",None): 1.2548385,
        ("hcp","N",None): 1.28257109,

        ("ontop","C","Cu"): 2.056904761,
        ("bridge","C","Cu"): 1.920118777,
        ("fcc","C","Cu"): 1.795024649,
        ("hcp","C","Cu"): 1.764312871,
        ("ontop","O","Cu"): 1.872603673,
        ("bridge","O","Cu"): 1.69205958,
        ("fcc","O","Cu"): 1.408365497,
        ("hcp","O","Cu"): 1.510464567,
        ("ontop","H","Cu"): 1.5496025,
        ("fcc","H","Cu"): 1.0321708,
        ("hcp","H","Cu"): 1.0321708,
        ("fcc","N","Cu"): 1.2548385,
        ("hcp","N","Cu"): 1.28257109,
}

def get_site_bond_length(sitetype,atomtype=None,metal=None):
    if "fold" in sitetype:
        sitetype = "fcc"
    if sitetype == "longbridge" or sitetype == "shortbridge":
        sitetype = "bridge"
    if (sitetype,atomtype,metal) in site_bond_length_dict.keys():
        return site_bond_length_dict[(sitetype,atomtype,metal)]
    elif (sitetype,atomtype,None) in site_bond_length_dict.keys():
        return site_bond_length_dict[(sitetype,atomtype,None)]
    else:
        return site_bond_length_dict[(sitetype,None,None)]

def add_adsorbate_to_site(atoms, adsorbate, surf_ind, site, height=None,
                          orientation=None, tilt_angle=0.):
    """The base function for adding one adsorbate to a site.
    Site must include information of 'normal' and 'position'.
    Useful for adding adsorbate to multiple sites or adding
    multidentate adsorbates.

    Parameters
    ----------
    atoms : ase.Atoms object
        Accept any ase.Atoms object. No need to be built-in.

    adsorbate : str or ase.Atom object or ase.Atoms object
        The adsorbate species to be added onto the surface.

    site : dict
        The site that the adsorbate should be added to.
        Must contain information of the position and the
        normal vector of the site.

    height : float, default None
        The height of the added adsorbate from the surface.
        Use the default settings if not specified.

    orientation : list or numpy.array, default None
        The vector that the multidentate adsorbate is aligned to.

    tilt_angle: float, default None
        Tilt the adsorbate with an angle (in degrees) relative to
        the surface normal.

    """
    if height is None:
        height = site_heights[site['site']]

    # Make the correct position
    normal = np.array(site['normal'])
    if np.isnan(np.sum(normal)):
        normal = np.array([0., 0., 1.])
    pos = np.array(site['position']) + normal * height

    # Convert the adsorbate to an Atoms object
    if isinstance(adsorbate, Atoms):
        ads = adsorbate
    elif isinstance(adsorbate, Atom):
        ads = Atoms([adsorbate])

    # Or assume it is a string representing a molecule
    else:
        ads = adsorbate_molecule(adsorbate)
        if not ads:
            warnings.warn('Nothing is added.')
            return

    bondpos = ads[surf_ind].position
    ads.translate(-bondpos)
    z = -1. if adsorbate in ['CH','NH','OH','SH'] else 1.
    ads.rotate(np.asarray([0., 0., z]) - bondpos, normal)
    if tilt_angle > 0.:
        pvec = np.cross(np.random.rand(3) - ads[0].position, normal)
        ads.rotate(tilt_angle, pvec, center=ads[0].position)

#     if adsorbate not in adsorbate_list:
#         # Always sort the indices the same order as the input symbol.
#         # This is a naive sorting which might cause H in wrong order.
#         # Please sort your own adsorbate atoms by reindexing as has
#         # been done in the adsorbate_molecule function in acat.settings.
#         symout = list(Formula(adsorbate))
#         symin = list(ads.symbols)
#         newids = []
#         for elt in symout:
#             idx = symin.index(elt)
#             newids.append(idx)
#             symin[idx] = None
#         ads = ads[newids]
    if orientation is not None:
        orientation = np.asarray(orientation)
        oripos = next((a.position for a in ads[1:] if
                       a.symbol != 'H'), ads[1].position)

        v1 = get_rejection_between(oripos - bondpos, normal)
        v2 = get_rejection_between(orientation, normal)
        theta = get_angle_between(v1, v2)

        # Flip the sign of the angle if the result is not the closest
        rm_p = get_rodrigues_rotation_matrix(axis=normal, angle=theta)
        rm_n = get_rodrigues_rotation_matrix(axis=normal, angle=-theta)
        npos_p, npos_n = rm_p @ oripos, rm_n @ oripos
        nbpos_p = npos_p + pos - bondpos
        nbpos_n = npos_n + pos - bondpos
        d_p = np.linalg.norm(nbpos_p - pos - orientation)
        d_n = np.linalg.norm(nbpos_n - pos - orientation)
        if d_p <= d_n:
            for a in ads:
                a.position = rm_p @ a.position
        else:
            for a in ads:
                a.position = rm_n @ a.position

    ads.translate(pos - bondpos)
    # Randomly offsetting atoms to avoid highly symmetric structures
    for atom in ads:
        x_trans = random.choice([-0.05, 0.05])
        y_trans = random.choice([-0.05, 0.05])
        atom.position[0] += x_trans
        atom.position[1] += y_trans

    atoms += ads
    if ads.get_chemical_formula() == 'H2':
        shift = (atoms.positions[-2] - atoms.positions[-1]) / 2
        atoms.positions[-2:,:] += shift

def place_adsorbate_covdep(ads,slab,atom_surf_inds,sites,metal):
    if len(atom_surf_inds) == 1:
        geo = slab.copy()
        h = get_site_bond_length(sites[0]["site"],ads.get_chemical_symbols()[atom_surf_inds[0]],metal)
        add_adsorbate_to_site(geo, ads, atom_surf_inds[0], sites[0], height=h)
        return geo,h,None
    elif len(atom_surf_inds) == 2:
        geo = slab.copy()
        h1 = get_site_bond_length(sites[0]["site"],ads.get_chemical_symbols()[atom_surf_inds[0]],metal)
        h2 = get_site_bond_length(sites[1]["site"],ads.get_chemical_symbols()[atom_surf_inds[1]],metal)
        ori = get_mic(sites[0]['position'], sites[1]['position'], geo.cell)
        add_adsorbate_to_site(geo, deepcopy(ads), atom_surf_inds[0], sites[0], height=h1, orientation=ori)
        if np.isnan(geo.positions).any(): #if nans just ignore orientation and let it optimize
            geo = slab.copy()
            add_adsorbate_to_site(geo, deepcopy(ads), atom_surf_inds[0], sites[0], height=h1, orientation=None)
        return geo,h1,h2
    else:
        raise ValueError
    
def get_unique_sites(site_list, cell, unique_composition=False,
                         unique_subsurf=False,
                         return_signatures=False,
                         return_site_indices=False,
                         about=None):
        """Function mostly copied from the ACAT software 
        Get all symmetry-inequivalent adsorption sites (one
        site for each type).

        Parameters
        ----------
        unique_composition : bool, default False
            Take site composition into consideration when
            checking uniqueness.

        unique_subsurf : bool, default False
            Take subsurface element into consideration when
            checking uniqueness.

        return_signatures : bool, default False
            Whether to return the unique signatures of the
            sites instead.

        return_site_indices: bool, default False
            Whether to return the indices of each unique
            site (in the site list).

        about: numpy.array, default None
            If specified, returns unique sites closest to
            this reference position.

        """ 
        sl = site_list[:]
        key_list = ['site', 'morphology']
        if unique_composition:
            key_list.append('composition')
            if unique_subsurf:
                key_list.append('subsurf_element')
        else:
            if unique_subsurf:
                raise ValueError('to include the subsurface element, ' +
                                 'unique_composition also need to be set to True')

        seen_tuple = []
        uni_sites = []
        if about is not None:
            sl = sorted(sl, key=lambda x: get_mic(x['position'],
                        about, cell, return_squared_distance=True))
        for i, s in enumerate(sl):
            sig = tuple(s[k] for k in key_list)
            if sig not in seen_tuple:
                seen_tuple.append(sig)
                if return_site_indices:
                    s = i
                uni_sites.append(s)

        return uni_sites

def generate_unique_placements(slab,sites):
    nslab = len(slab)
    middle = sum(slab.cell)/2.0

    unique_single_sites = get_unique_sites(sites,slab.cell,about=middle)

    unique_site_pairs = dict() # (site1site,site1morph,),(site2site,site2morph),xydist,zdist
    for unique_site in unique_single_sites:
        uni_site_fingerprint = (unique_site["site"],unique_site["morphology"])
        for site in sites:
            site_fingerprint = (site["site"],site["morphology"])
            bd,d = get_distances([unique_site["position"]], [site["position"]], cell=slab.cell, pbc=(True,True,False))
            xydist = np.linalg.norm(bd[0][0][:1])
            zdist = bd[0][0][2]

            fingerprint = (uni_site_fingerprint,site_fingerprint,round(xydist,3),round(zdist,3))

            if fingerprint in unique_site_pairs.keys():
                current_sites = unique_site_pairs[fingerprint]
                current_dist = np.linalg.norm(sum([s["position"][:1] for s in current_sites])/2-middle[:1])
                possible_dist = np.linalg.norm((unique_site["position"][:1]+site["position"][:1])/2-middle[:1])
                if possible_dist < current_dist:
                    unique_site_pairs[fingerprint] = [unique_site,site]
            else:
                unique_site_pairs[fingerprint] = [unique_site,site]

    unique_site_pairs_lists = list(unique_site_pairs.values())
    unique_site_lists = [[unique_site] for unique_site in unique_single_sites]

    single_site_bond_params_lists = []
    for unique_site_list in unique_site_lists:
        pos = deepcopy(unique_site_list[0]["position"])
        single_site_bond_params_lists.append([{"site_pos": pos,"ind": None, "k": 100.0, "deq": 0.0}])

    double_site_bond_params_lists = []
    for unique_site_pair_list in unique_site_pairs_lists:
        bond_params_list = []
        for site in unique_site_pair_list:
            pos = deepcopy(site["position"])
            bond_params_list.append({"site_pos": pos,"ind": None, "k": 100.0, "deq": 0.0})
        double_site_bond_params_lists.append(bond_params_list)

    return unique_site_lists,unique_site_pairs_lists,single_site_bond_params_lists,double_site_bond_params_lists

def generate_unique_site_additions(geo,sites,slab,nslab,site_bond_params_list=[],sites_list=[]):
    nads = len(geo) - nslab
    #label sites with unique noble gas atoms
    he = Atoms('He',positions=[[0, 0, 0],])
    ne = Atoms('Ne',positions=[[0, 0, 0],])
    ar = Atoms('Ar',positions=[[0, 0, 0],])
    kr = Atoms('Kr',positions=[[0, 0, 0],])
    xe = Atoms('Xe',positions=[[0, 0, 0],])
    rn = Atoms('Rn',positions=[[0, 0, 0],])
    site_tags = [he,ne,ar,kr,xe,rn]
    tag = site_tags[nads]
    occ = get_occupied_sites(geo,sites,nslab)
    unocc = [site for site in sites if not any(sites_match(site,osite,slab) for osite in occ)]

    site_bond_params_lists = [deepcopy(site_bond_params_list) for i in range(len(unocc))]
    sites_lists = [deepcopy(sites_list) for i in range(len(unocc))]

    geoms = []
    for i,site in enumerate(unocc):
        geom = geo.copy()
        add_adsorbate_to_site(geom,adsorbate=tag,surf_ind=0,site=site,height=1.5)
        pos = site["position"].tolist()
        params = {"site_pos": pos,"ind": None, "k": 100.0, "deq": 0.0} #just the site position, will shift up later
        site_bond_params_lists[i].append(params)
        sites_lists[i].append(site)
        geoms.append(geom)

    indclusters = get_unique_sym_struct_index_clusters(geoms)

    inds = []
    for cluster in indclusters: #choose symmetric geometries closer to center of cell
        min_dist = np.inf
        indout = None
        for ind in cluster:
            d = get_adsorbate_dist_from_center(geoms[ind],nslab)
            if d < min_dist:
                indout = ind
                min_dist = d
        inds.append(indout)

    return [geoms[ind] for ind in inds], [site_bond_params_lists[ind] for ind in inds], [sites_lists[ind] for ind in inds]

def get_adsorbate_dist_from_center(atoms,nslab):
    cell_center = sum(atoms.cell)/3.5
    adatoms = atoms[nslab:]
    adcenter = sum([a.position for a in adatoms])/len(adatoms)
    return np.linalg.norm(adcenter - cell_center)


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

def get_labeled_bonds(mol):
    """
    generate a list of all Bonds between labeled atoms in the Molecule object
    """
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
    """
    compares the reactant and product structures assuming bonds only form and break between labeled atoms
    returns frozensets of pairs of labels associated with atoms whose bond breaks or is formed during the reaction
    """
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
    for tempmol in tempmols:
        if tempmol.multiplicity == -187: #handle surface molecules
            tempmol.multiplicity = tempmol.get_radical_count() + 1

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
    """
    process the subgraph isomorphisms found to generate a list of valid maps
    """
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
    """
    get number of atoms in the adsorbate part of the Molecule
    """
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
    """
    get the index associated with the ase.Atoms object
    """
    n = nslab
    for i,d in enumerate(template_mol_map):
        if ind in d.keys() and d[ind] in molecule_to_gratom_maps[i].keys():
            return molecule_to_gratom_maps[i][d[ind]]+n
        n += ads_sizes[i]
    else:
        return None #surface site

def get_bond_lengths_sites(mol,ads,atom_map,surf_atom_map,nslab,sites,site_adjacency,facet="fcc111",metal="Cu",):
    """
    gets bond lengths and site information indexed to the Molecule object atoms
    bond lengths is a matrix with the distances between all atoms that share bonds
    site information is the atom index mapped to the site type
    """
    rev_atom_map = {value: key for key,value in atom_map.items()} #map from mol to ads
    rev_surf_atom_map = { value: key for key,value in surf_atom_map.items()}

    occ = get_occupied_sites(ads,sites,nslab)
    surface_dict = [{"atom_index":x["bonding_index"]-nslab, "site":x["site"],
                       "bond_length":x["bond_length"], "position":x["position"]} for x in occ]

    if len(occ) < len(surf_atom_map): #number of sites on geometry disagrees with surf_atom_map
#         print("occupational analysis in get_bond_lengths_sites failed")
#         print(mol.to_adjacency_list())
#         print("expected {} sites".format(len(surf_atom_map)))
#         print("found {} sites".format(len(occ)))
        return None,None,None

    if mol.contains_surface_site():
        ad = ads[nslab:]
    else:
        ad = ads

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

def get_name(mol):
    try:
        return mol.to_smiles()
    except:
        return mol.to_adjacency_list().replace("\n"," ")[:-1].replace(' ','')


def remove_slab(mol,remove_slab_bonds=False,update_atomtypes=True):
    m = mol.copy(deep=True)
    for site in m.get_surface_sites():
        for a in site.bonds.keys():
            if not a.is_surface_site():
                break
        else:
            m.remove_atom(site)
            
    if remove_slab_bonds:
        bonds_to_remove = []
        for bd in m.get_all_edges():
            if bd.atom1.is_surface_site() and bd.atom2.is_surface_site():
                m.remove_bond(bd)
    
    if update_atomtypes:           
        m.update_atomtypes()
    m.update_connectivity_values()
    
    return m

def pluck_subgraph(mol,atom):
    subgraph_atoms = []
    new_subgraph_atoms = [atom]
    while new_subgraph_atoms != []:
        temp = []
        subgraph_atoms.extend(new_subgraph_atoms)
        for a in new_subgraph_atoms:
            for a2 in a.bonds.keys():
                if a2 not in subgraph_atoms:
                    temp.append(a2)
        new_subgraph_atoms = temp
    
    inds = [mol.atoms.index(a) for a in subgraph_atoms]
    struct = Molecule(atoms=subgraph_atoms)
    struct.update_atomtypes()
    struct.update_connectivity_values()
    return struct,inds
    
def generate_without_site_info(m):
    mol = m.copy(deep=True)
    for a in mol.atoms:
        if a.is_surface_site():
            a.site = ""
            a.morphology = ""
    return mol

class FindingPathError(Exception):
    pass

def reduce_graph_to_pairs(admol):
    adatoms = [a for a in admol.atoms if a.is_bonded_to_surface() and not a.is_surface_site()]
    surface_save = set()
    for i,a1 in enumerate(adatoms):
        for a2 in adatoms[i+1:]:
            paths = find_shortest_paths(a1,a2)
            if paths is None: #separation between pair is so large the generated site graphs don't connect
                raise FindingPathError(admol.to_adjacency_list())
            for path in paths:
                surface_save = surface_save | set([a for a in path if a.is_surface_site()])
    
    atoms_to_remove = []
    for i,a in enumerate(admol.atoms):
        if a.is_surface_site() and not (a in surface_save):
            atoms_to_remove.append(a)
    
    for a in atoms_to_remove:
        admol.remove_atom(a)

    admol.update_atomtypes()
    admol.update_connectivity_values()
    
    return admol

def find_shortest_paths(start, end, path=None):
    paths = [[start]]
    outpaths = []
    while paths != []:
        newpaths = []
        for path in paths:
            for node in path[-1].edges.keys():
                if node in path:
                    continue
                elif node is end:
                    outpaths.append(path[:]+[node])
                elif outpaths == []:
                    newpaths.append(path[:]+[node])
        if outpaths:
            return outpaths
        else:
            paths = newpaths
    
    return None

def find_shortest_paths_sites(start, end, path=None):
    paths = [[start]]
    outpaths = []
    while paths != []:
        newpaths = []
        for path in paths:
            for node in path[-1].edges.keys():
                if node in path:
                    continue
                elif node is end:
                    outpaths.append(path[:]+[node])
                elif outpaths == [] and node.is_surface_site():
                    newpaths.append(path[:]+[node])
        if outpaths:
            return outpaths
        else:
            paths = newpaths
    
    return None
    
def find_adsorbate_atoms_surface_sites(atom,mol):
    """
    one atom of the associated adsorbate on the surface Molecule object mol
    returns all of the surface atoms bonded to adsorbate
    """
    atsout = [atom]
    ats = [atom]
    surf_sites = []
    while ats != []:
        new_ats = []
        for a in ats:
            for a2 in atom.edges.keys():
                if a2 in atsout:
                    continue
                elif a2.is_surface_site() and a2 not in surf_sites:
                    surf_sites.append(a2)
                elif not a2.is_surface_site():
                    new_ats.append(a2)
                    atsout.append(a2)
        ats = new_ats

    return atsout,surf_sites

def split_adsorbed_structures(admol,clear_site_info=True,adsorption_info=False,atoms_to_skip=None,
                              atom_mapping=False, split_sites_with_multiple_adsorbates=False):
    """_summary_

    Args:
        admol (_type_): a Molecule object resolving a slab
        clear_site_info (bool, optional): clears site identify information. Defaults to True.
        adsorption_info (bool, optional): also returns the map from atom to index for admol for each surface site structures are split off from. Defaults to False.
        atoms_to_skip (_type_, optional): list of atoms in admol to not include in the split structures. Defaults to None.
        split_sites_with_multiple_adsorbates (bool, optional): If a site has multiple bonds to adsorbates splits that site so that each adsorption bond is to a different site
        
    Returns:
        _type_: _description_
    """
    m = admol.copy(deep=True)
    
    if atoms_to_skip is None:
        atoms_to_remove = []
        skip_atoms = []
    else:
        skip_atoms = [m.atoms[admol.atoms.index(a)] for a in atoms_to_skip]
        atoms_to_remove = skip_atoms[:]
    
    for i,at in enumerate(m.atoms):
        if at.is_surface_site():
            bdict = m.get_bonds(at)
            if len(bdict) > 0:
                for a,b in bdict.items():
                    if not a.is_surface_site() and (a not in skip_atoms):
                        break
                else:
                    atoms_to_remove.append(at)
            else:
                atoms_to_remove.append(at)

    reaction_bonds = False
    bd_to_remove = []
    for bd in m.get_all_edges():
        if bd.is_reaction_bond():
            reaction_bonds = True
        if bd.atom1.is_surface_site() and bd.atom2.is_surface_site():
            bd_to_remove.append(bd)
    
    split_site_atom_mapping = dict()
    new_site_atoms = []
    if split_sites_with_multiple_adsorbates:
        ad_ind_to_original_site_atom_bond_map = dict()
        adsorbate_new_site_inds = []
        for i,a in enumerate(m.atoms):
            if a.is_surface_site() and len([k for k in a.bonds.keys() if not k.is_surface_site()]) > 1:
                ind = max(m.atoms.index(a2) for a2 in a.bonds.keys())
                adsorbate_new_site_inds.append(ind)
                bd = a.bonds[m.atoms[ind]]
                ad_ind_to_original_site_atom_bond_map[ind] = (a,bd)
                bd_to_remove.append(bd)

        for ind in adsorbate_new_site_inds:
            orig_site,orig_bond = ad_ind_to_original_site_atom_bond_map[ind]
            a = m.atoms[ind]
            newat = Atom(element="X", lone_pairs=0, site=a.site, morphology=a.morphology)
            m.add_atom(newat)
            new_site_atoms.append(newat)
            bd = Bond(a,newat,order=orig_bond.order)
            m.add_bond(bd)
            split_site_atom_mapping[orig_site] = [m.atoms.index(orig_site),m.atoms.index(newat)]
        
    if adsorption_info:
        adsorbed_atom_dict = {at: i for i,at in enumerate(m.atoms) if at.is_surface_site() and at not in atoms_to_remove and at not in new_site_atoms}
        for k,v in split_site_atom_mapping.items():
            adsorbed_atom_dict[k] = v

    if atom_mapping:
        mapping = {at: i for i,at in enumerate(m.atoms) if at not in atoms_to_remove and at not in new_site_atoms}
        for k,v in split_site_atom_mapping.items():
            mapping[k] = v
    
    for bd in bd_to_remove:
        m.remove_bond(bd)
        
    for at in atoms_to_remove:
        m.remove_atom(at) 
    
    if clear_site_info:
        for at in m.atoms:
            if at.is_surface_site():
                at.site = ""
                at.morphology = ""

    split_structs = m.split()
    for newmol in split_structs:
        if reaction_bonds:
            newmol.update_multiplicity()
            newmol.identify_ring_membership()
            newmol.update_connectivity_values()
        else:
            newmol.update(sort_atoms=False)
            newmol.update_connectivity_values()

    if atom_mapping and not adsorption_info:
        return split_structs,mapping 
    elif atom_mapping and adsorption_info:
        return split_structs,adsorption_info,mapping 
    elif not atom_mapping and adsorption_info:
        return split_structs,adsorbed_atom_dict
    else:
        return split_structs

def get_full_mol(admol,slab_mol,site_map):
    """
    Generates the full Molecule adsorbate/s and slab representation from a unresolved adsorbate/s and fully resolved slab representation
    and an atom mapping between the adsorbates/s sites and slab representation sites
    admol is a Molecule object without the slab resolved with all adsorbates/TS bound to sites where bound to the surface
    slab_mol is a Molecule object representing the full periodic slab
    site_map is a dictionary mapping atoms in admol to slab_mol
    """
    index_map = {admol.atoms.index(k):slab_mol.atoms.index(v) for k,v in site_map.items()}
    slabm = slab_mol.copy(deep=True)
    admolm = admol.copy(deep=True)
    smap = {admolm.atoms[k]:slabm.atoms[v] for k,v in index_map.items()}
    mmol = slabm.merge(admolm)
    for adsite,slabsite in smap.items():
        if adsite.label:
            slabsite.label = adsite.label
        adatom = list(adsite.bonds.keys())[0]
        adbd = mmol.get_bond(adsite,adatom)
        order = adbd.order
        mmol.add_bond(Bond(adatom,slabsite,order=order))
        mmol.remove_bond(adbd)

    for adsite in smap.keys():
        mmol.remove_atom(adsite)

    mmol.multiplicity = mmol.get_radical_count() + 1
    return mmol

def get_labeled_full_TS_mol(template,mol):
    """
    takes in a template with only the labels that participate in reactions
    and an unlabeled slab resolved Molecule object 
    """
    m = mol.copy(deep=True)
    m_single = m.to_single_bonds(raise_atomtype_exception=False)
    m_single.multiplicity = 1
    
    tempmol_mol_map = []
    tempmols = [x for x in template.split() if not x.is_surface_site()]
    tempgrps = []
    for tempmol in tempmols:
        if tempmol.multiplicity == -187: #handle surface molecules
            tempmol.multiplicity = tempmol.get_radical_count() + 1
        tgp = tempmol.to_single_bonds(raise_atomtype_exception=False).to_group()
        for i in range(len(tempmol.atoms)):
            tgp.atoms[i].label = tempmol.atoms[i].label
        for a in tgp.atoms:
            if a.is_surface_site():
                a.atomtype = [ATOMTYPES["X"]]
        tgp.multiplicity = [1]
        tempgrps.append(tgp)

    map_list = []
    for i,tempgrp in enumerate(tempgrps):
        mapvs = m_single.find_subgraph_isomorphisms(tempgrp,save_order=True)
        if len(mapvs) == 0 and len(tempgrp.get_surface_sites()) == 2: #the molecule may bind two atoms to one site
            sites = tempgrp.get_surface_sites()
            for at,bd in sites[1].bonds.items():
                if not tempgrp.has_bond(at,sites[0]):
                    b = GroupBond(sites[0],at,order=bd.order)
                    tempgrp.add_bond(b)
            tempgrp.remove_atom(sites[1])
            mapvs = m_single.find_subgraph_isomorphisms(tempgrp,save_order=True)
        if len(mapvs) == 0:
            logging.error(m_single.to_adjacency_list())
            logging.error(tempgrp.to_adjacency_list())
            raise IndexError
        map_list.append([{m.atoms[m_single.atoms.index(k)]: v for k,v in x.items()} for x in mapvs]) #remap back to proper not single bonded Molecule

    fullmaps = []
    for pd in itertools.product(*map_list):
        key_atoms = sum([list(x.keys()) for x in pd],[])
        value_atoms = sum([list(x.values()) for x in pd],[])
        if len(set(key_atoms)) == len(key_atoms) and len(set(value_atoms)) == len(value_atoms):
            fullmaps.append({k:v for x in pd for k,v in x.items()})

    unique_fullmaps = []
    for mapm in fullmaps:
        for map2 in unique_fullmaps:
            if mapm == map2:
                break
        else:
            unique_fullmaps.append(mapm)
    
    labeled_mols = []
    for mapm in unique_fullmaps:
        labeled_mol = m.copy(deep=True)
        for k,v in mapm.items():
            labeled_mol.atoms[m.atoms.index(k)].label = v.label
        labeled_mols.append(labeled_mol)

    return labeled_mols