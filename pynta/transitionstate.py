from molecule.molecule import Molecule, Bond
from molecule.molecule.pathfinder import find_shortest_path
import yaml
from ase.io import read, write
import numpy as np
from ase.data import covalent_radii
from acat.adsorption_sites import SlabAdsorptionSites
from acat.adsorbate_coverage import SlabAdsorbateCoverage
from acat.settings import adsorbate_molecule
from acat.utilities import (custom_warning,
                         is_list_or_tuple,
                         get_close_atoms,
                         get_rodrigues_rotation_matrix,
                         get_angle_between,
                         get_rejection_between)
import os
from pynta.utils import get_unique_sym, get_unique_sym_structs
from ase.visualize import view
from ase.atoms import Atoms, Atom
from ase.geometry import get_distances
import itertools
from pynta.calculator import HarmonicallyForcedXTB
from pynta.mol import *
from copy import deepcopy
from pysidt import *
import pynta.models

def get_unique_optimized_adsorbates(rxn,adsorbates_path,mol_dict,cas,gratom_to_molecule_surface_atom_maps,nslab):
    """
    load the adsorbates associated with the reaction and find the unique optimized
    adsorbate structures for each species
    returns a dictionary mapping each adsorbate name to a list of ase.Atoms objects
    """
    adsorbates = dict()

    for name in rxn["reactant_names"]+rxn["product_names"]:
        prefixes = os.listdir(os.path.join(adsorbates_path,name))
        ase_to_mol_surface_atom_map = gratom_to_molecule_surface_atom_maps[name]
        geoms = []
        for prefix in prefixes:
            path = os.path.join(adsorbates_path,name,prefix,prefix+".xyz")
            if os.path.exists(path):
                geoms.append(path)
        xyzs = get_unique_sym(geoms)
        adsorbates[name] = []
        for xyz in xyzs:
            geo = read(xyz)
            adcov = SlabAdsorbateCoverage(geo,adsorption_sites=cas)
            sites = adcov.get_sites()
            occ = [site for site in sites if site["occupied"]]
            required_surface_inds = set([ind+nslab for ind in ase_to_mol_surface_atom_map.keys()])
            found_surface_inds = set([site["bonding_index"] for site in occ])
            if len(occ) >= len(mol_dict[name].get_adatoms()) and required_surface_inds.issubset(found_surface_inds):
                adsorbates[name].append(geo)

    return adsorbates

def determine_TS_construction(reactant_names,reactant_mols,product_names,product_mols):
    """
    determine the direction to generate the TS guess in and the order in which the
    reactants will be placed on the surface
    returns forward a boolean indicating if the reaction should be estimated in the forward (True)
    or reverse (False) direction and ordered_reacting_species which is a list of species in the order
    they should be placed on the surface
    """
    forward = None
    rnum_surf_sites = [len(mol.get_surface_sites()) for i,mol in enumerate(reactant_mols)]
    pnum_surf_sites = [len(mol.get_surface_sites()) for i,mol in enumerate(product_mols)]

    rnummultidentate = len([i for i in rnum_surf_sites if i>1])
    pnummultidentate = len([i for i in pnum_surf_sites if i>1])

    rnumgas = len([i for i in rnum_surf_sites if i==0])
    pnumgas = len([i for i in pnum_surf_sites if i==0])

    #avoid directions starting with multiple multidentate species
    if rnummultidentate >= 2 and pnummultidentate >= 2:
        raise ValueError("Cannot handle multiple bidentates on both sides of reaction")
    elif rnummultidentate >= 2:
        forward = False
    elif pnummultidentate >= 2:
        forward = True

    #prefer directions with multidentate species, this reduces the time and cost of finding TS guesses
    if forward is None:
        if rnummultidentate == 1 and pnummultidentate == 0:
            forward = True
        elif pnummultidentate == 1 and rnummultidentate == 0:
            forward = False

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
                inds_to_order.remove(indout)
                continue

            monodentateinds = [ind for ind in inds_to_order if rnum_surf_sites[ind] == 1]
            if len(monodentateinds) > 0:
                #biggest monodentate first
                sorted_inds = sorted(monodentateinds,key=lambda ind: len(reactant_mols[ind].atoms))
                indout = sorted_inds[-1]
                ordered_reacting_species.append(reactant_names[indout])
                inds_to_order.remove(indout)
                continue

            #gas phase biggest first
            sorted_inds = sorted(inds_to_order,key=lambda ind: len(reactant_mols[ind].atoms))
            indout = sorted_inds[-1]
            ordered_reacting_species.append(reactant_names[indout])
            inds_to_order.remove(indout)
            continue
    else:
        inds_to_order = list(range(len(product_names)))
        while len(ordered_reacting_species) < len(product_names):
            bidentateinds = [ind for ind in inds_to_order if pnum_surf_sites[ind] > 1]
            if len(bidentateinds) > 0:
                indout = bidentateinds[0] #should only be one if any
                ordered_reacting_species.append(product_names[indout])
                inds_to_order.remove(indout)
                continue

            monodentateinds = [ind for ind in inds_to_order if pnum_surf_sites[ind] == 1]

            if len(monodentateinds) > 0:
                #biggest monodentate first
                sorted_inds = sorted(monodentateinds,key=lambda ind: len(product_mols[ind].atoms))
                indout = sorted_inds[-1]
                ordered_reacting_species.append(product_names[indout])
                inds_to_order.remove(indout)
                continue

            #gas phase biggest first
            sorted_inds = sorted(inds_to_order,key=lambda ind: len(product_mols[ind].atoms))
            indout = sorted_inds[-1]
            ordered_reacting_species.append(product_names[indout])
            inds_to_order.remove(indout)
            continue

    return forward,ordered_reacting_species


def get_unique_TS_structs(adsorbates,species_names,slab,cas,nslab,num_surf_sites,mol_dict,
                          gratom_to_molecule_atom_maps,gratom_to_molecule_surface_atom_maps,
                          facet,metal,gas_height=5.0):
    """
    Generate unique initial structures for TS guess generation
    """
    tsstructs = []
    slab_sites = cas.get_sites()
    site_adjacency = cas.get_neighbor_site_list()
    ordered_adsorbates = [adsorbates[name] for name in species_names]
    for adss in itertools.product(*ordered_adsorbates):
        if num_surf_sites[0] > 0:
            adslab = adss[0].copy()
        else:
            adslab = slab.copy()
            site = cas.get_sites()[0]
            add_adsorbate_to_site(adslab,adsorbate=adss[0],surf_ind=0,site=site,height=gas_height)
        if len(adss) == 1:
            tsstructs.append(adslab)
        else:
            if num_surf_sites[1] == 1:
                name = species_names[1]
                ad = adss[1][nslab:]
                bondlengths1,sites1,site_lengths1 = get_bond_lengths_sites(mol_dict[name],adss[1],gratom_to_molecule_atom_maps[name],
                                                                        gratom_to_molecule_surface_atom_maps[name],nslab,
                                                                        slab_sites,site_adjacency,facet=facet,metal=metal)
                sitetype1 = list(sites1.values())[0]
                height1 = list(site_lengths1.values())[0]
                surf_ind1 = list(gratom_to_molecule_surface_atom_maps[name].keys())[0]
                adcov = SlabAdsorbateCoverage(adslab,adsorption_sites=cas)
                for site in adcov.get_sites():
                    adslab2 = adslab.copy()
                    if site['occupied'] == False and site["site"] == sitetype1:
                        add_adsorbate_to_site(adslab2,adsorbate=ad,surf_ind=surf_ind1,site=site,height=height1)
                        if len(adss) == 2:
                            tsstructs.append(adslab2)
                        else:
                            if num_surf_sites[2] == 1:
                                name2 = species_names[2]
                                ad2 = adss[2][nslab:]
                                bondlengths2,sites2,site_lengths2 = get_bond_lengths_sites(mol_dict[name2],adss[2],gratom_to_molecule_atom_maps[name2],
                                                                        gratom_to_molecule_surface_atom_maps[name2],nslab,
                                                                        slab_sites,site_adjacency,facet=facet,metal=metal)
                                sitetype2 = list(sites2.values())[0]
                                height2 = list(site_lengths2.values())[0]
                                surf_ind2 = list(gratom_to_molecule_surface_atom_maps[name2].keys())[0]
                                adcov2 = SlabAdsorbateCoverage(adslab2,adsorption_sites=cas)
                                for site2 in adcov2.get_sites():
                                    adslab3 = adslab2.copy()
                                    if site2['occupied'] == False and site2["site"] == sitetype2:
                                        add_adsorbate_to_site(adslab3,adsorbate=ad2,surf_ind=surf_ind2,site=site2,height=height2)
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
                                add_adsorbate_to_site(adslab,adsorbate=adss[2],surf_ind=0,site=site,height=gas_height)

            elif num_surf_sites[1] == 0:
                adcovl1 = SlabAdsorbateCoverage(adslab,adsorption_sites=cas)
                sites = adcovl1.get_sites()
                c = 0
                site = sites[c]
                while site["occupied"] == True:
                    c += 1
                    site = sites[c]

                add_adsorbate_to_site(adslab,adsorbate=adss[1],surf_ind=0,site=site,height=gas_height)
                if len(adss) == 2:
                    tsstructs.append(adslab)
                else:
                    c = 0
                    site2 = sites[c]
                    while site2["occupied"] == True and site2 != site:
                        c += 1
                        site2 = sites[c]
                    add_adsorbate_to_site(adslab,adsorbate=adss[2],surf_ind=0,site=site2,height=gas_height)
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
    """
    Generate constraints and harmonic parameters for the harmonically forced optimization
    """
    constraint_lists = []
    atom_bond_potential_lists = []
    site_bond_potential_lists = []
    site_bond_dict_list = []
    out_structs = []
    slab_sites = cas.get_sites()
    site_adjacency = cas.get_neighbor_site_list()
    
    nodes_file = os.path.join(os.path.split(pynta.models.__file__)[0],"TS_tree_nodes.json")
    nodes = read_nodes(nodes_file)
    sidt = SubgraphIsomorphicDecisionTree(nodes=nodes)
    
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
                                                                nslab,slab_sites,site_adjacency,facet=facet,metal=metal)
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
                                                                nslab,slab_sites,site_adjacency,facet=facet,metal=metal)
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
        tsstruct_valid = True
        atom_bond_potentials = []
        site_bond_potentials = []
        fixed_bond_pairs = []
        site_bond_dict = dict()

        adcov = SlabAdsorbateCoverage(tsstruct,adsorption_sites=cas)
        sites = adcov.get_sites()
        occ = [site for site in sites if site["occupied"]]
        if len(occ) < len(forward_template.get_adatoms()): #if we don't identify as many occupied sites as we should skip this structure
            print("occupied not equal to forward_template")
            print("occ={0}".format(len(occ)))
            print("forward_template_occs={0}".format(len(forward_template.get_adatoms())))
            continue
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
            if ((ase_ind1 is None) and (ase_ind2 is None)):
                print(ind1)
                print(ind2)
                print(nslab)
                print(ads_sizes)
                raise ValueError
            labels = set([label1,label2])
            if labels in broken_bonds: #bonds that break
                if ase_ind1 and ase_ind2:
                    dwell = tsstruct.get_distance(ase_ind1,ase_ind2,mic=True)
                    sitetype = None
                    pos = None
                elif ase_ind1 is None:
                    assert edge.atom1.is_surface_site()
                    if ase_ind2 not in occ_bd_lengths.keys():
                        print("acat occupational detection has failed for a TS guess structure...may be skipping important TS guesses")
                        tsstruct_valid = False
                        break
                    dwell = occ_bd_lengths[ase_ind2]
                    sitetype = occ_site_types[ase_ind2]
                    pos = deepcopy(occ_site_pos[ase_ind2])
                    pos[2] += dwell
                    ind = ase_ind2
                else:
                    assert edge.atom2.is_surface_site()
                    if ase_ind1 not in occ_bd_lengths.keys():
                        print("acat occupational detection has failed for a TS guess structure...may be skipping important TS guesses")
                        tsstruct_valid = False
                        break
                    dwell = occ_bd_lengths[ase_ind1]
                    sitetype = occ_site_types[ase_ind1]
                    pos = deepcopy(occ_site_pos[ase_ind1])
                    pos[2] += dwell
                    ind = ase_ind1
                deq,k = estimate_deq_k(sidt,labels,dwell,forward_template,reverse_template,template_name,template_reversed,broken_bonds,formed_bonds,sitetype=sitetype)
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
                    pos = deepcopy(occ_site_pos[ase_ind2])
                    pos[2] += dwell #shift position up to atom position
                    ind = ase_ind2
                    deq,k = estimate_deq_k_fixed_surf_bond(labels,dwell,forward_template,reverse_template,template_name,template_reversed,broken_bonds,formed_bonds,sitetype=sitetype)
                    d = {"ind":ind,"site_pos":pos.tolist(),"k":k,"deq":deq}
                    site_bond_potentials.append(d)
                else:
                    dwell = occ_bd_lengths[ase_ind1]
                    sitetype = occ_site_types[ase_ind1]
                    pos = deepcopy(occ_site_pos[ase_ind1])
                    pos[2] += dwell #shift position up to atom position
                    ind = ase_ind1
                    deq,k = estimate_deq_k_fixed_surf_bond(labels,dwell,forward_template,reverse_template,template_name,template_reversed,broken_bonds,formed_bonds,sitetype=sitetype)
                    d = {"ind":ind,"site_pos":pos.tolist(),"k":k,"deq":deq}
                    site_bond_potentials.append(d)

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
                deq,k = estimate_deq_k(sidt,labels,dwell,forward_template,reverse_template,template_name,template_reversed,broken_bonds,formed_bonds,sitetype=sitetype)
                d = {"ind1":ase_ind1,"ind2":ase_ind2,"k":k,"deq":deq}
                atom_bond_potentials.append(d)
            elif ase_ind1 is None: #surface bond to ase_ind2
                site_dict = rev_mols_info[ind1_mol_r]["site_dict"]
                reusing_site = False

                #check if site reused
                if len(forward_template.atoms[ind1f].bonds) > 0:
                    reusing_site = True
                    atm = forward_template.atoms[ind1f]
                    bd = list(forward_template.atoms[ind1f].bonds.values())[0]
                    if bd.atom1 == atm:
                        bonded_atom = bd.atom2
                    else:
                        bonded_atom = bd.atom1
                    template_bonded_index = forward_template.atoms.index(bonded_atom)
                    ase_bonded_index = get_ase_index(template_bonded_index,template_mol_map,molecule_to_gratom_maps,nslab,ads_sizes)
                    reused_site_type = occ_site_types[ase_bonded_index]
                    reused_site_pos = occ_site_pos[ase_bonded_index]

                if reusing_site and (reused_site_type not in [x[1] for x in site_dict.keys()]): #the site one of the reactants is attached to is not stable in the products
                    print("{reused_site} not in list of sites {site_list}".format(reused_site=reused_site_type,site_list=[x[1] for x in site_dict.keys()]))
                    tsstruct_valid = False
                    break

                site_bond_dict[ase_ind2] = dict()
                for (mol_ind,sitetype) in site_dict.keys():
                    if mol_ind == ind1r_mol:
                        dwell = site_dict[(mol_ind,sitetype)]
                        deq,k = estimate_deq_k(sidt,labels,dwell,forward_template,reverse_template,template_name,template_reversed,broken_bonds,formed_bonds,sitetype=sitetype)
                        if reusing_site: #attaching to a site that was bonded to another atom
                            if reused_site_type == sitetype:
                                pos = deepcopy(reused_site_pos)
                                pos[2] += dwell
                                d = {"ind":ase_ind2,"site_pos":pos.tolist(),"k":k,"deq":deq}
                                site_bond_potentials.append(d)
                                if ase_ind2 in site_bond_dict.keys():
                                    del site_bond_dict[ase_ind2]
                        else: #attaching to a new empty site
                            site_bond_dict[ase_ind2][sitetype] = {"deq":deq,"k":k,"dwell":dwell}

            else:
                site_dict = rev_mols_info[ind2_mol_r]["site_dict"]
                reusing_site = False

                #check if site reused
                if len(forward_template.atoms[ind2f].bonds) > 0:
                    reusing_site = True
                    atm = forward_template.atoms[ind2f]
                    bd = list(forward_template.atoms[ind2f].bonds.values())[0]
                    if bd.atom1 == atm:
                        bonded_atom = bd.atom2
                    else:
                        bonded_atom = bd.atom1
                    template_bonded_index = forward_template.atoms.index(bonded_atom)
                    ase_bonded_index = get_ase_index(template_bonded_index,template_mol_map,molecule_to_gratom_maps,nslab,ads_sizes)
                    reused_site_type = occ_site_types[ase_bonded_index]
                    reused_site_pos = occ_site_pos[ase_bonded_index]

                if reusing_site and (reused_site_type not in [x[1] for x in site_dict.keys()]): #the site one of the reactants is attached to is not stable in the products
                    print("{reused_site} not in list of sites {site_list}".format(reused_site=reused_site_type,site_list=[x[1] for x in site_dict.keys()]))
                    tsstruct_valid = False
                    break

                site_bond_dict[ase_ind1] = dict()
                for (mol_ind,sitetype) in site_dict.keys():
                    if mol_ind == ind2r_mol:
                        dwell = site_dict[(mol_ind,sitetype)]
                        deq,k = estimate_deq_k(sidt,labels,dwell,forward_template,reverse_template,template_name,template_reversed,broken_bonds,formed_bonds,sitetype=sitetype)
                        if reusing_site: #attaching to a site that was bonded to another atom
                            if reused_site_type == sitetype:
                                pos = deepcopy(reused_site_pos)
                                pos[2] += dwell
                                d = {"ind":ase_ind1,"site_pos":pos.tolist(),"k":k,"deq":deq}
                                site_bond_potentials.append(d)
                                if ase_ind1 in site_bond_dict.keys():
                                    del site_bond_dict[ase_ind1]
                        else:
                            site_bond_dict[ase_ind1][sitetype] = {"deq":deq,"k":k,"dwell":dwell}

        if tsstruct_valid:
            if fixed_bond_pairs:
                constraint_list = [{"type": "FixBondLength", "a1": pair[0], "a2": pair[1]} for pair in fixed_bond_pairs]+["freeze slab"]
                constraint_lists.append(constraint_list)
            else:
                constraint_lists.append(["freeze slab"])
            atom_bond_potential_lists.append(atom_bond_potentials)
            site_bond_potential_lists.append(site_bond_potentials)
            site_bond_dict_list.append(site_bond_dict)
            out_structs.append(tsstruct)

    return out_structs,constraint_lists,atom_bond_potential_lists,site_bond_potential_lists,site_bond_dict_list

def estimate_deq_k(sidt,labels,dwell,forward_template,reverse_template,template_name,template_reversed,
    broken_bonds,formed_bonds,sitetype=None):
    """
    Estimate the equilibrium bond length and force constant for broken/forming bonds
    0--(1--2)--3
    """
    interactions = broken_bonds | formed_bonds
    formed = labels in formed_bonds
    mol = deepcopy(forward_template)
    all_labels = set()
    for inter in interactions:
        label_list = list(inter)
        all_labels.add(label_list[0])
        all_labels.add(label_list[1])
        a1 = mol.get_labeled_atoms(label_list[0])[0]
        a2 = mol.get_labeled_atoms(label_list[1])[0]
        if mol.has_bond(a1,a2):
            bd = mol.get_bond(a1,a2)
            bd.set_order_str('R')
        else:
            mol.add_bond(Bond(a1,a2,order='R'))
    
    for label in all_labels:
        a = mol.get_labeled_atoms(label)[0]
        if label in labels:
            a.label = '*'
        else:
            a.label = ''
    
    v = sidt.evaluate(mol)
    
    if any([a.is_surface_site() for a in mol.get_labeled_atoms('*')]):
        return dwell*(v-1.0),100.0
    else:
        return dwell*v,100.0

def estimate_deq_k_fixed_surf_bond(labels,dwell,forward_template,reverse_template,template_name,template_reversed,
    broken_bonds,formed_bonds,sitetype=None):
    """
    Estimate the equilibrium bond length and force constant for surface bonds
    """
    return 0.0,100.0

def get_surface_forming_bond_pairings(tsstructs,atom_bond_potential_lists,
                                      site_bond_potential_lists,constraints_list,site_bond_dict_list,
                                      cas):
    """
    Identify unique site targets for the harmonically forced optimization
    """
    if len(site_bond_dict_list[0]) == 0: #no site bonds formed
        return tsstructs,atom_bond_potential_lists,site_bond_potential_lists,constraints_list

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
    site_tags = {atm_ind : site_tags[i] for i,atm_ind in enumerate(atm_inds)}

    new_atom_bond_potential_lists = []
    new_site_bond_potential_lists = []
    new_constraints_list = []
    out_tsstructs = []

    for i,tsstruct in enumerate(tsstructs):
        atom_bond_potential_list = atom_bond_potential_lists[i]
        constraints = constraints_list[i]
        adcov = SlabAdsorbateCoverage(tsstruct,adsorption_sites=cas)
        unocc = [site for site in adcov.get_sites() if site["occupied"] == False]
        site_bond_dict = site_bond_dict_list[i]
        site_bond_potential_dict = sites_to_site_bond_potentials(unocc,site_bond_dict,atm_inds)
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
                    add_adsorbate_to_site(geo,adsorbate=tag,surf_ind=0,site=site,height=1.5)
                    params = bond_dict[site["site"]].copy()
                    params["site_pos"] = site["position"].tolist()
                    params["site_pos"][2] += params["dwell"]
                    del params["dwell"]
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

    return out_tsstructs,new_atom_bond_potential_lists,new_site_bond_potential_lists,new_constraints_list

def sites_to_site_bond_potentials(sites,site_bond_dicts,atm_inds):
    """
    Generate a list of possible site bond formation harmonic parameters for each
    atom that forms bonds with a site
    """
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
