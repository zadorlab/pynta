from pynta.excatkit.adsorption import Builder
from pynta.excatkit.gratoms import Gratoms
from molecule.molecule import Molecule
from pynta.io import IO
import yaml
from ase.io import read, write
import numpy as np
from ase.data import covalent_radii
from acat.adsorption_sites import SlabAdsorptionSites
from acat.adsorbate_coverage import SlabAdsorbateCoverage
from acat.settings import site_heights
import os
from pynta.symmetry import get_unique_sym, get_unique_sym_structs
from ase.visualize import view
from ase.atoms import Atoms
from ase.geometry import get_distances
import itertools
from pynta.calculator import HarmonicallyForcedXTB
from pynta.molecule import *

def get_unique_optimized_adsorbates(rxn,adsorbates_path):
    adsorbates = dict()

    for name in rxn["reactant_names"]+rxn["product_names"]:
        prefixes = os.listdir(os.path.join(adsorbates_path,name))
        geoms = []
        for prefix in prefixes:
            path = os.path.join(adsorbates_path,name,prefix,prefix+".xyz")
            if os.path.exists(path):
                geoms.append(path)
        xyzs = get_unique_sym(geoms)
        adsorbates[name] = [read(xyz) for xyz in xyzs]


    return adsorbates

def determine_TS_construction(reactant_names,reactant_mols,product_names,product_mols):
    forward = None
    rnum_surf_sites = [len(mol.get_surface_sites()) for i,mol in enumerate(reactant_mols)]
    pnum_surf_sites = [len(mol.get_surface_sites()) for i,mol in enumerate(product_mols)]

    rnummultidentate = len([i for i in rnum_surf_sites if i>1])
    pnummultidentate = len([i for i in pnum_surf_sites if i>1])

    rnumgas = len([i for i in rnum_surf_sites if i==0])
    pnumgas = len([i for i in pnum_surf_sites if i==0])

    #avoid directions starting with multiple multidentate species
    if rnummultidentate > 2 and pnummultidentate > 2:
        raise ValueError("Cannot handle multiple bidentates on both sides of reaction")
    elif rnummultidentate > 2:
        forward = False
    elif pnummultidentate > 2:
        forward = True

    #prefer directions with fewer gas phase species
    if forward is None:
        if rnumgas > pnumgas:
            forward = False
        elif pnumgas > rnumgas:
            forward = True

    #prefer directions with fewer species
    if forward is None:
        if len(reactant_mols) > len(product_mols):
            forward = False
        elif len(reactant_mols) < len(product_mols):
            forward = True

    if forward is None:
        forward = True

    ordered_reacting_species = []
    if forward:
        inds_to_order = list(range(len(reactant_names)))
        while len(ordered_reacting_species) < len(reactant_names):
            bidentateinds = [ind for ind in inds_to_order if rnum_surf_sites[ind] > 1]
            if len(bidentateinds) > 0:
                indout = bidentateinds[0] #should only be one if any
                ordered_reacting_species.append(reactant_names[indout])
                inds_to_order.pop(indout)
                continue

            monodentateinds = [ind for ind in inds_to_order if rnum_surf_sites[ind] == 1]
            if len(monodentateinds) > 0:
                #biggest monodentate first
                sorted_inds = sorted(monodentateinds,key=lambda ind: len(reactant_mols[ind].atoms))
                indout = sorted_inds[-1]
                ordered_reacting_species.append(reactant_names[indout])
                inds_to_order.pop(indout)
                continue

            #gas phase biggest first
            sorted_inds = sorted(inds_to_order,key=lambda ind: len(reactant_mols[ind].atoms))
            indout = sorted_inds[-1]
            ordered_reacting_species.append(reactant_names[indout])
            inds_to_order.pop(indout)
            continue
    else:
        inds_to_order = list(range(len(product_names)))
        while len(ordered_reacting_species) < len(product_names):
            bidentateinds = [ind for ind in inds_to_order if pnum_surf_sites[ind] > 1]
            if len(bidentateinds) > 0:
                indout = bidentateinds[0] #should only be one if any
                ordered_reacting_species.append(product_names[indout])
                inds_to_order.pop(indout)
                continue

            monodentateinds = [ind for ind in inds_to_order if pnum_surf_sites[ind] == 1]
            if len(monodentateinds) > 0:
                #biggest monodentate first
                sorted_inds = sorted(monodentateinds,key=lambda ind: len(product_mols[ind].atoms))
                indout = sorted_inds[-1]
                ordered_reacting_species.append(product_names[indout])
                inds_to_order.pop(indout)
                continue

            #gas phase biggest first
            sorted_inds = sorted(inds_to_order,key=lambda ind: len(product_mols[ind].atoms))
            indout = sorted_inds[-1]
            ordered_reacting_species.append(product_names[indout])
            inds_to_order.pop(indout)
            continue

    return forward,ordered_reacting_species

def add_adsorbate_to_site(atoms, adsorbate, site, height=None,
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
        warnings.warn('The normal vector is NaN, use [0., 0., 1.] instead.')
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

    bondpos = ads[0].position
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

def get_unique_TS_structs(adsorbates,species_names,cas,nslab,num_surf_sites,mol_dict,
                          gratom_to_molecule_atom_maps,gratom_to_molecule_surface_atom_maps,
                          facet,metal,gas_height=5.0):
    tsstructs = []
    ordered_adsorbates = [adsorbates[name] for name in species_names]
    for adss in itertools.product(*ordered_adsorbates):
        adslab = adss[0]
        if len(adss) == 1:
            tsstructs.append(adslab)
        else:
            if num_surf_sites[1] == 1:
                name = species_names[1]
                ad = adss[1][nslab:]
                bondlengths1,sites1,site_lengths1 = get_bond_lengths_sites(mol_dict[name],adss[1],gratom_to_molecule_atom_maps[name],
                                                                        gratom_to_molecule_surface_atom_maps[name],nslab,
                                                                        facet=facet,metal=metal,cas=cas)
                sitetype1 = list(sites1.values())[0]
                height1 = list(site_lengths1.values())[0]
                adcov = SlabAdsorbateCoverage(adslab,adsorption_sites=cas)
                for site in adcov.get_sites():
                    adslab2 = adslab.copy()
                    if site['occupied'] == False and site["site"] == sitetype1:
                        add_adsorbate_to_site(adslab2,adsorbate=ad,site=site,height=height1)
                        if len(adss) == 2:
                            tsstructs.append(adslab2)
                        else:
                            if num_surf_sites[2] == 1:
                                name = species_names[2]
                                ad2 = adss[2][nslab:]
                                bondlengths2,sites2,site_lengths2 = get_bond_lengths_sites(mol_dict[name],adss[2],gratom_to_molecule_atom_maps[name],
                                                                        gratom_to_molecule_surface_atom_maps[name],nslab,
                                                                        facet=facet,metal=metal,cas=cas)
                                sitetype2 = list(sites2.values())[0]
                                height2 = list(site_lengths2.values())[0]
                                adcov2 = SlabAdsorbateCoverage(adslab2,adsorption_sites=cas)
                                for site2 in adcov2.get_sites():
                                    adslab3 = adslab2.copy()
                                    if site2['occupied'] == False and site["site"] == sitetype2:
                                        add_adsorbate_to_site(adslab3,adsorbate=ad2,site=site2,height=height2)
                                        if len(adss) == 3:
                                            tsstructs.append(adslab3)
                                        else:
                                            raise ValueError("Cannot handle more than three reactants")
                            else:
                                adcovl2 = SlabAdsorbateCoverage(adslab2,adsorption_sites=cas)
                                sites = adcovl2.get_sites()
                                c = 0
                                site = sites[c]
                                while site["occupied"] == True:
                                    c += 1
                                    site = sites[c]
                                add_adsorbate_to_site(adslab,adsorbate=adss[2],site=site,height=gas_height)

            elif num_surf_sites[1] == 0:
                adcovl1 = SlabAdsorbateCoverage(adslab,adsorption_sites=cas)
                sites = adcovl1.get_sites()
                c = 0
                site = sites[c]
                while site["occupied"] == True:
                    c += 1
                    site = sites[c]

                add_adsorbate_to_site(adslab,adsorbate=adss[1],site=site,height=gas_height)
                if len(adss) == 2:
                    tsstructs.append(adslab)
                else:
                    c = 0
                    site2 = sites[c]
                    while site2["occupied"] == True and site2 != site:
                        c += 1
                        site2 = sites[c]
                    add_adsorbate_to_site(adslab,adsorbate=adss[2],site=site)
                    if len(adss) == 3:
                        tsstructs.append(adslab)
                    else:
                        raise ValueError("Cannot handle more than three reactants")

    unique_tsstructs = get_unique_sym_structs(tsstructs)
    return unique_tsstructs

def generate_constraints_harmonic_parameters(tsstructs,adsorbates,slab,forward_template,
                                             reverse_template,template_name,template_reversed,
                                            ordered_names,reverse_names,mol_dict,gratom_to_molecule_atom_maps,
                                            gratom_to_molecule_surface_atom_maps,nslab,facet,metal,cas):
    constraint_lists = []
    atom_bond_potential_lists = []
    site_bond_potential_lists = []
    site_bond_dict_list = []
    site_fixed_bond_dict_list = []

    mols = [mol_dict[name] for name in ordered_names]
    rev_mols = [mol_dict[name] for name in reverse_names]
    ads_sizes = [ads_size(mol) for mol in mols]
    mols_rev = [mol_dict[name] for name in reverse_names]
    template_mol_map = get_template_mol_map(forward_template,mols)
    reverse_template_mol_map = get_template_mol_map(reverse_template,rev_mols)
    molecule_to_gratom_maps = [{value:key for key,value in gratom_to_molecule_atom_maps[name].items()} for name in ordered_names]
    broken_bonds,formed_bonds = get_broken_formed_bonds(forward_template,reverse_template)

    mols_info = []
    for name in ordered_names:
        bdlength_list = []
        sites_list = []
        site_lengths_list = []
        for ads in adsorbates[name]:
            bondlengths,sites,site_lengths = get_bond_lengths_sites(mol_dict[name],ads,
                                                                gratom_to_molecule_atom_maps[name],
                                                                gratom_to_molecule_surface_atom_maps[name],
                                                                nslab,facet=facet,metal=metal,cas=cas)
            bdlength_list.append(bondlengths)
            sites_list.append(sites)
            site_lengths_list.append(site_lengths)

        ave_bdlength = sum(bdlength_list)/len(bdlength_list)
        site_dict = dict()
        for i,sites in enumerate(sites_list):
            for key,site in sites.items():
                if (key,site) in site_dict.keys():
                    site_dict[(key,site)].append(site_lengths_list[i][key])
                else:
                    site_dict[(key,site)] = [site_lengths_list[i][key]]

        site_dict = {key: np.mean(value) for key,value in site_dict.items()}
        mols_info.append({"bondlengths":ave_bdlength,"site_dict":site_dict})

    rev_mols_info = []
    for name in reverse_names:
        bdlength_list = []
        sites_list = []
        site_lengths_list = []
        for ads in adsorbates[name]:
            bondlengths,sites,site_lengths = get_bond_lengths_sites(mol_dict[name],ads,
                                                                gratom_to_molecule_atom_maps[name],
                                                                gratom_to_molecule_surface_atom_maps[name],
                                                                nslab,facet=facet,metal=metal,cas=cas)
            bdlength_list.append(bondlengths)
            sites_list.append(sites)
            site_lengths_list.append(site_lengths)

        ave_bdlength = sum(bdlength_list)/len(bdlength_list)
        site_dict = dict()
        for i,sites in enumerate(sites_list):
            for key,site in sites.items():
                if (key,site) in site_dict.keys():
                    site_dict[(key,site)].append(site_lengths_list[i][key])
                else:
                    site_dict[(key,site)] = [site_lengths_list[i][key]]

        site_dict = {key: np.mean(value) for key,value in site_dict.items()}
        rev_mols_info.append({"bondlengths":ave_bdlength,"site_dict":site_dict})

    for tsstruct in tsstructs:
        atom_bond_potentials = []
        site_bond_potentials = []
        fixed_bond_pairs = []
        site_bond_dict = dict()
        site_fixed_bond_dict = dict()

        adcov = SlabAdsorbateCoverage(tsstruct,adsorption_sites=cas)
        sites = adcov.get_sites()
        occ = [site for site in sites if site["occupied"]]
        occ_atom_inds = [site["bonding_index"] for site in occ]
        occ_bd_lengths = {site["bonding_index"]:site['bond_length'] for site in occ}
        occ_site_types = {site["bonding_index"]:site['site'] for site in occ}
        occ_site_pos = {site["bonding_index"]:site['position'] for site in occ}
        unocc = [site for site in sites if not site["occupied"]]

        for edge in forward_template.get_all_edges():
            ind1 = forward_template.atoms.index(edge.atom1)
            ind2 = forward_template.atoms.index(edge.atom2)
            label1 = edge.atom1.label
            label2 = edge.atom2.label
            ase_ind1 = get_ase_index(ind1,template_mol_map,molecule_to_gratom_maps,nslab,ads_sizes)
            ase_ind2 = get_ase_index(ind2,template_mol_map,molecule_to_gratom_maps,nslab,ads_sizes)
            labels = set([label1,label2])
            if labels in broken_bonds: #bonds that break
                if ase_ind1 and ase_ind2:
                    dwell = tsstruct.get_distance(ase_ind1,ase_ind2,mic=True)
                    sitetype = None
                    pos = None
                elif ase_ind1 is None:
                    assert edge.atom1.is_surface_site()
                    dwell = occ_bd_lengths[ase_ind2]
                    sitetype = occ_site_types[ase_ind2]
                    pos = occ_site_pos[ase_ind2]
                    ind = ase_ind2
                else:
                    assert edge.atom2.is_surface_site()
                    dwell = occ_bd_lengths[ase_ind1]
                    sitetype = occ_site_types[ase_ind1]
                    pos = occ_site_pos[ase_ind1]
                    ind = ase_ind1
                deq,k = estimate_deq_k(labels,dwell,forward_template,reverse_template,template_name,template_reversed,sitetype=sitetype)
                if pos is not None:
                    d = {"ind":ind,"site_pos":pos.tolist(),"k":k,"deq":deq}
                    site_bond_potentials.append(d)
                else:
                    d = {"ind1":ase_ind1,"ind2":ase_ind2,"k":k,"deq":deq}
                    atom_bond_potentials.append(d)
            else: #bond that doesn't break
                if ase_ind1 and ase_ind2:
                    fixed_bond_pairs.append([ase_ind1,ase_ind2])
                elif ase_ind1 is None: #surface bonds that don't break
                    dwell = occ_bd_lengths[ase_ind2]
                    sitetype = occ_site_types[ase_ind2]
                    pos = occ_site_pos[ase_ind2]
                    ind = ase_ind2
                    deq,k = estimate_deq_k_fixed_surf_bond(labels,dwell,forward_template,reverse_template,template_name,template_reversed,sitetype=sitetype)
                    d = {"ind":ind,"site_pos":pos.tolist(),"k":k,"deq":deq}
                    site_bond_potentials.append(d)

                    site_fixed_bond_dict[ase_ind2] = dict()
                    ind_mol,ind1_mol = get_mol_index(ind1,template_mol_map)
                    site_dict = mols_info[ind_mol]["site_dict"]
                    for (mol_ind,sitetype) in site_dict.keys():
                        if mol_ind == ind1_mol:
                            dwell = site_dict[(mol_ind,sitetype)]
                            deq,k = estimate_deq_k_fixed_surf_bond(labels,dwell,forward_template,reverse_template,template_name,template_reversed,sitetype=sitetype)
                            site_fixed_bond_dict[ase_ind2][sitetype] = {"deq":deq,"k":k}
                else:
                    dwell = occ_bd_lengths[ase_ind1]
                    sitetype = occ_site_types[ase_ind1]
                    pos = occ_site_pos[ase_ind1]
                    ind = ase_ind1
                    deq,k = estimate_deq_k_fixed_surf_bond(labels,dwell,forward_template,reverse_template,template_name,template_reversed,sitetype=sitetype)
                    d = {"ind":ind,"site_pos":pos.tolist(),"k":k,"deq":deq}
                    site_bond_potentials.append(d)

                    site_fixed_bond_dict[ase_ind1] = dict()
                    ind_mol,ind2_mol = get_mol_index(ind2,template_mol_map)
                    site_dict = mols_info[ind_mol]["site_dict"]
                    for (mol_ind,sitetype) in site_dict.keys():
                        if mol_ind == ind2_mol:
                            dwell = site_dict[(mol_ind,sitetype)]
                            deq,k = estimate_deq_k_fixed_surf_bond(labels,dwell,forward_template,reverse_template,template_name,template_reversed,sitetype=sitetype)
                            site_fixed_bond_dict[ase_ind1][sitetype] = {"deq":deq,"k":k}

        for labels in formed_bonds: #bonds that form
            atmsf = []
            for label in labels:
                atmsf.extend(forward_template.get_labeled_atoms(label))
            indsf = [forward_template.atoms.index(atm) for atm in atmsf]
            ind1f = indsf[0]
            ind2f = indsf[1]
            ase_ind1 = get_ase_index(ind1f,template_mol_map,molecule_to_gratom_maps,nslab,ads_sizes)
            ase_ind2 = get_ase_index(ind2f,template_mol_map,molecule_to_gratom_maps,nslab,ads_sizes)
            atmsr = []
            for label in labels:
                atmsr.extend(reverse_template.get_labeled_atoms(label))
            indsr = [reverse_template.atoms.index(atm) for atm in atmsr]
            ind1r = indsr[0]
            ind2r = indsr[1]
            ind1_mol_r,ind1r_mol = get_mol_index(ind1r,reverse_template_mol_map)
            ind2_mol_r,ind2r_mol = get_mol_index(ind2r,reverse_template_mol_map)
            assert ind1_mol_r == ind2_mol_r

            if ase_ind1 and ase_ind2: #bond between atoms
                dwell = rev_mols_info[ind1_mol_r]["bondlengths"][ind1r_mol,ind2r_mol]
                sitetype = None
                deq,k = estimate_deq_k(labels,dwell,forward_template,reverse_template,template_name,template_reversed,sitetype=sitetype)
                d = {"ind1":ase_ind1,"ind2":ase_ind2,"k":k,"deq":deq}
                atom_bond_potentials.append(d)
            elif ase_ind1 is None: #surface bond to ase_ind2
                site_bond_dict[ase_ind2] = dict()
                site_dict = rev_mols_info[ind1_mol_r]["site_dict"]
                for (mol_ind,sitetype) in site_dict.keys():
                    if mol_ind == ind1r_mol:
                        dwell = site_dict[(mol_ind,sitetype)]
                        deq,k = estimate_deq_k(labels,dwell,forward_template,reverse_template,template_name,template_reversed,sitetype=sitetype)
                        site_bond_dict[ase_ind2][sitetype] = {"deq":deq,"k":k}

            else:
                site_bond_dict[ase_ind1] = dict()
                site_dict = rev_mols_info[ind2_mol_r]["site_dict"]
                for (mol_ind,sitetype) in site_dict.keys():
                    if mol_ind == ind2r_mol:
                        dwell = site_dict[(mol_ind,sitetype)]
                        deq,k = estimate_deq_k(labels,dwell,forward_template,reverse_template,template_name,template_reversed,sitetype=sitetype)
                        site_bond_dict[ase_ind1][sitetype] = {"deq":deq,"k":k}

        if fixed_bond_pairs:
            constraint_lists.append([{"type": "FixBondLengths","pairs":fixed_bond_pairs},"freeze slab"])
        else:
            constraint_lists.append(["freeze slab"])
        atom_bond_potential_lists.append(atom_bond_potentials)
        site_bond_potential_lists.append(site_bond_potentials)
        site_bond_dict_list.append(site_bond_dict)
        site_fixed_bond_dict_list.append(site_fixed_bond_dict)

    return constraint_lists,atom_bond_potential_lists,site_bond_potential_lists,site_bond_dict_list,site_fixed_bond_dict_list

def estimate_deq_k(labels,dwell,forward_template,reverse_template,template_name,template_reversed,sitetype=None):
    return dwell*1.4,100.0

def estimate_deq_k_fixed_surf_bond(labels,dwell,forward_template,reverse_template,template_name,template_reversed,sitetype=None):
    return dwell,100.0

def get_surface_forming_bond_pairings(tsstructs,atom_bond_potential_lists,
                                      site_bond_potential_lists,constraints_list,site_bond_dict_list,
                                      site_fixed_bond_dict_list,cas):

    if len(site_bond_dict_list[0]) == 0: #no site bonds formed
        return tsstructs,site_bond_potential_lists,None

    adslablen = len(tsstructs[0])

    #label sites with unique noble gas atoms
    he = Atoms('He',positions=[[0, 0, 0],])
    ne = Atoms('Ne',positions=[[0, 0, 0],])
    ar = Atoms('Ar',positions=[[0, 0, 0],])
    kr = Atoms('Kr',positions=[[0, 0, 0],])
    xe = Atoms('Xe',positions=[[0, 0, 0],])
    rn = Atoms('Rn',positions=[[0, 0, 0],])
    site_tags = [he,ne,ar,kr,xe,rn]

    atm_inds = list(site_bond_dict_list[0].keys())
    fixed_atm_inds = list(site_fixed_bond_dict_list[0].keys())
    site_tags = {atm_ind : site_tags[i] for i,atm_ind in enumerate(atm_inds)}

    new_atom_bond_potential_lists = []
    new_site_bond_potential_lists = []
    new_constraints_list = []
    site_bond_potential_check_lists = []
    out_tsstructs = []

    for i,tsstruct in enumerate(tsstructs):
        atom_bond_potential_list = atom_bond_potential_lists[i]
        constraints = constraints_list[i]
        adcov = SlabAdsorbateCoverage(tsstruct,adsorption_sites=cas)
        unocc = [site for site in adcov.get_sites() if site["occupied"] == False]
        site_bond_dict = site_bond_dict_list[i]
        site_fixed_bond_dict = site_fixed_bond_dict_list[i]
        site_bond_potential_dict = sites_to_site_bond_potentials(unocc,site_bond_dict,atm_inds)
        site_fixed_bond_potential_dict = sites_to_site_bond_potentials(unocc,site_fixed_bond_dict,fixed_atm_inds)
        geoms = []
        site_parameters = []

        pd = [unocc for i in range(len(atm_inds))]
        for sites in itertools.product(*pd):
            if len(set([repr(x) for x in sites])) != len(sites) or any([sites[j]["site"] not in site_bond_dict[atm_ind].keys() for j,atm_ind in enumerate(atm_inds)]):
                continue
            else:
                param_list = []
                geo = tsstruct.copy()
                for j,atm_ind in enumerate(atm_inds):
                    tag = site_tags[atm_ind]
                    bond_dict = site_bond_dict[atm_ind]
                    site = sites[j]
                    add_adsorbate_to_site(geo,adsorbate=tag,site=site,height=1.5)
                    params = bond_dict[site["site"]].copy()
                    params["site_pos"] = site["position"].tolist()
                    params["ind"] = atm_ind
                    param_list.append(params)

                site_parameters.append(param_list)
                geoms.append(geo)

        unique_geoms = get_unique_sym_structs(geoms)
        unique_site_parameters = []

        for geom in unique_geoms:
            ind = geoms.index(geom)
            out_tsstructs.append(geom[:adslablen])
            new_atom_bond_potential_lists.append(atom_bond_potential_list)
            new_constraints_list.append(constraints)

            params = site_parameters[ind]

            site_bond_potential = site_bond_potential_lists[i].copy()
            site_bond_potential.extend(params)

            new_site_bond_potential_lists.append(site_bond_potential)

            site_bond_potential_dict_temp = site_bond_potential_dict.copy()
            for atm_ind in atm_inds:
                for param in params:
                    if param["ind"] == atm_ind:
                        k = 0
                        while any([site_bond_potential_dict_temp[atm_ind][k]["site_pos"][v] != param["site_pos"][v] for v in range(3)]):
                            k += 1
                        del site_bond_potential_dict_temp[atm_ind][k]

            site_fixed_bond_potential_dict_temp = site_fixed_bond_potential_dict.copy()

            site_bond_potential_dict_sum = {**site_bond_potential_dict_temp,**site_fixed_bond_potential_dict_temp}
            site_bond_potential_check_lists.append(site_bond_potential_dict_sum)




    return out_tsstructs,new_atom_bond_potential_lists,new_site_bond_potential_lists,new_constraints_list,site_bond_potential_check_lists

def sites_to_site_bond_potentials(sites,site_bond_dicts,atm_inds):
    site_bond_potentials = dict()
    for atm_ind in atm_inds:
        params_list = []
        for site in sites:
            bond_dict = site_bond_dicts[atm_ind]
            if site["site"] in bond_dict.keys():
                params = bond_dict[site["site"]].copy()
                params["site_pos"] = site["position"].tolist()
                params["ind"] = atm_ind
                params_list.append(params)
        site_bond_potentials[atm_ind] = params_list

    return site_bond_potentials
