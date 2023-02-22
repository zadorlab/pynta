from molecule.molecule import Molecule
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
from pynta.utils import get_unique_sym_struct_index_clusters, get_unique_sym, get_unique_sym_structs, get_unique_sym_struct_indices
from pynta.calculator import run_harmonically_forced_xtb
from rdkit import Chem
from copy import deepcopy
import numpy as np

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

def generate_adsorbate_guesses(mol,ads,slab,repeats,cas,mol_to_atoms_map,metal,
                               single_site_bond_params_lists,single_sites_lists,double_site_bond_params_lists,double_sites_lists,
                               Eharmtol,Eharmfiltertol,Ntsmin):
    full_slab = slab * repeats
    mol_surf_inds = [mol.atoms.index(a) for a in mol.get_adatoms()]
    atom_surf_inds = [mol_to_atoms_map[i] for i in mol_surf_inds]
    if len(atom_surf_inds) == 1:
        site_bond_params_lists = deepcopy(single_site_bond_params_lists)
        sites_lists = single_sites_lists
        for site_bond_params_list in site_bond_params_lists:
            site_bond_params_list[0]["ind"] = atom_surf_inds[0]+len(full_slab)

            #add up pulling potential
            for ind in range(len(ads)):
                if ind in atom_surf_inds:
                    continue
                pos = deepcopy(site_bond_params_list[0]["site_pos"])
                pos[2] += 8.5
                site_bond_params_list.append({"site_pos": pos,"ind": ind+len(full_slab), "k": 0.1, "deq": 0.0})

    elif len(atom_surf_inds) == 2:
        site_bond_params_lists = deepcopy(double_site_bond_params_lists)
        sites_lists = double_sites_lists
        for site_bond_params_list in site_bond_params_lists:
            site_bond_params_list[0]["ind"] = atom_surf_inds[0]+len(full_slab)
            site_bond_params_list[1]["ind"] = atom_surf_inds[1]+len(full_slab)
    else:
        raise ValueError("Only monodentate and bidentate guesses currently allowed. The infrastructure can support tridenate and higher, but the filtering process may be very expensive above bidentate.")


    mol_fixed_bond_pairs = [[mol.atoms.index(bd.atom1),mol.atoms.index(bd.atom2)] for bd in mol.get_all_edges() if (not bd.atom1.is_surface_site()) and (not bd.atom2.is_surface_site())]
    atom_fixed_bond_pairs = [[mol_to_atoms_map[pair[0]]+len(full_slab),mol_to_atoms_map[pair[1]]+len(full_slab)]for pair in mol_fixed_bond_pairs]
    constraint_list = [{"type": "fix_bond", "indices": pair} for pair in atom_fixed_bond_pairs]+["freeze slab"]

    geos = []
    for i,sites_list in enumerate(sites_lists):
        geo,h1,h2 = place_adsorbate(ads,full_slab,atom_surf_inds,sites_list,metal)
        if h1:
            site_bond_params_lists[i][0]["site_pos"][2] += h1
        if h2:
            site_bond_params_lists[i][1]["site_pos"][2] += h2
        geos.append(geo)

    print("initial geometries")
    print(len(geos))
    geos_out = []
    Eharms = []
    site_bond_params_lists_out = []
    for i,geo in enumerate(geos):
        #freeze bonds for messier first opt
        geo_out,Eharm,Fharm = run_harmonically_forced_xtb(geo,[],site_bond_params_lists[i],len(full_slab),
                                molecule_to_atom_maps=mol_to_atoms_map,ase_to_mol_num=None,
                                method="GFN1-xTB",constraints=constraint_list)
        if geo_out:
            geo_out.calc = None
            geos_out.append(geo_out)
            Eharms.append(Eharm)
            site_bond_params_lists_out.append(site_bond_params_lists[i])

    print("optimized geometries")
    print(len(geos_out))
    inds = get_unique_sym_struct_indices(geos_out)

    print("after symmetry")
    print(len(inds))

    geos_out = [geos_out[ind] for ind in inds]
    Eharms = [Eharms[ind] for ind in inds]

    if len(atom_surf_inds) == 1: #should be small, don't bother filtering
        xyzsout = geos_out
        site_bond_params_lists_final = [site_bond_params_lists_out[ind] for ind in inds]
        return xyzsout
    else:
        Einds = np.argsort(np.array(Eharms))
        Emin = np.min(np.array(Eharms))
        xyzsout = []
        site_bond_params_lists_final = []
        for Eind in Einds:
            if Eharms[Eind]/Emin < Eharmtol: #include all TSs with energy close to Emin
                xyzsout.append(geos_out[Eind])
                site_bond_params_lists_final.append(site_bond_params_lists_out[Eind])
            elif Eharms[Eind]/Emin > Eharmfiltertol: #if the energy is much larger than Emin skip it
                continue
            elif len(xyzsout) < Ntsmin: #if the energy isn't similar, but isn't much larger include the smallest until Ntsmin is reached
                xyzsout.append(geos_out[Eind])
                site_bond_params_lists_final.append(site_bond_params_lists_out[Eind])

    return xyzsout

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
    normal = site['normal']
    if np.isnan(np.sum(normal)):
        normal = np.array([0., 0., 1.])
    pos = site['position'] + normal * height

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
    atoms += ads
    if ads.get_chemical_formula() == 'H2':
        shift = (atoms.positions[-2] - atoms.positions[-1]) / 2
        atoms.positions[-2:,:] += shift

def place_adsorbate(ads,slab,atom_surf_inds,sites,metal):
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

def get_unique_sites(cas, unique_composition=False,
                         unique_subsurf=False,
                         return_signatures=False,
                         return_site_indices=False,
                         about=None, site_list=None):
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

        if site_list is None:
            sl = cas.site_list
        else:
            sl = site_list
        key_list = ['site', 'morphology']
        if unique_composition:
            if not cas.composition_effect:
                raise ValueError('the site list does not include '
                                 + 'information of composition')
            key_list.append('composition')
            if unique_subsurf:
                key_list.append('subsurf_element')
        else:
            if unique_subsurf:
                raise ValueError('to include the subsurface element, ' +
                                 'unique_composition also need to be set to True')
        if return_signatures:
            sklist = sorted([[s[k] for k in key_list] for s in sl])
            return sorted(list(sklist for sklist, _ in groupby(sklist)))
        else:
            seen_tuple = []
            uni_sites = []
            if about is not None:
                sl = sorted(sl, key=lambda x: get_mic(x['position'],
                            about, cas.cell, return_squared_distance=True))
            for i, s in enumerate(sl):
                sig = tuple(s[k] for k in key_list)
                if sig not in seen_tuple:
                    seen_tuple.append(sig)
                    if return_site_indices:
                        s = i
                    uni_sites.append(s)

            return uni_sites

def generate_unique_placements(slab,cas):
    nslab = len(slab)
    middle = sum(slab.cell)/2.0

    unique_single_sites = get_unique_sites(cas,about=middle)

    unique_site_pairs = dict() # (site1site,site1morph,),(site2site,site2morph),xydist,zdist
    for unique_site in unique_single_sites:
        uni_site_fingerprint = (unique_site["site"],unique_site["morphology"])
        for site in cas.get_sites():
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

def generate_unique_site_additions(geo,cas,nslab,site_bond_params_list=[],sites_list=[]):
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
    adcov = SlabAdsorbateCoverage(geo,adsorption_sites=cas)
    sites = adcov.get_sites()
    unocc = [site for site in sites if site["occupied"] == False]

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

    if len(occ) < len(surf_atom_map): #number of sites on geometry disagrees with surf_atom_map
        print("occupational analysis in get_bond_lengths_sites failed")
        print(mol.to_adjacency_list())
        print("expected {} sites".format(len(surf_atom_map)))
        print("found {} sites".format(len(occ)))
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
        return mol.to_adjacency_list().replace("\n"," ")[:-1]
