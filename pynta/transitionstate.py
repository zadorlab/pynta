from molecule.molecule import Molecule, Bond
import yaml
from ase.io import read
import numpy as np
import os
from pynta.utils import get_unique_sym, get_unique_sym_structs, sites_match
from pynta.geometricanalysis import *
import itertools
from pynta.mol import *
from copy import deepcopy
from pysidt.sidt import *
import pynta.models
import logging

def get_unique_optimized_adsorbates(rxn,adsorbates_path,mol_dict,gratom_to_molecule_surface_atom_maps,sites,nslab):
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
            occ = get_occupied_sites(geo,sites,nslab)
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

    rnum_vdW_species = len([i for i,mol in enumerate(reactant_mols) if any(bd.is_van_der_waals() for bd in mol.get_all_edges())])
    pnum_vdW_species = len([i for i,mol in enumerate(product_mols) if any(bd.is_van_der_waals() for bd in mol.get_all_edges())])
    
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

    #prefer directions with fewer species with vdW bonds
    if forward is None:
        if rnum_vdW_species > pnum_vdW_species:
            forward = False
        elif pnum_vdW_species > rnum_vdW_species:
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


def get_unique_TS_structs(adsorbates,species_names,slab,slab_sites,site_adjacency,nslab,num_surf_sites,mol_dict,
                          gratom_to_molecule_atom_maps,gratom_to_molecule_surface_atom_maps,
                          facet,metal,gas_height=5.0):
    """
    Generate unique initial structures for TS guess generation
    """
    tsstructs = []
    tsmols = []
    ordered_adsorbates = [adsorbates[name] for name in species_names]
    for adss in itertools.product(*ordered_adsorbates):
        if num_surf_sites[0] > 0:
            adslab = adss[0].copy()
            try:
                adslabmol,neighbor_sites,ninds = generate_adsorbate_2D(adslab, slab_sites, site_adjacency, nslab, max_dist=np.inf, cut_off_num=None, allowed_structure_site_structures=None,
                          keep_binding_vdW_bonds=True, keep_vdW_surface_bonds=any(bd.is_van_der_waals() for bd in mol_dict[species_names[0]].get_all_edges()))
            except (SiteOccupationException,FailedFixBondsException,TooManyElectronsException) as e:
                logging.warning("a configuration for {} was not convertable to 2D and has been skipped".format(species_names[0]))
                continue
        else:
            adslab = slab.copy()
            adslabmol,neighbor_sites,ninds = generate_adsorbate_2D(adslab, slab_sites, site_adjacency, nslab, max_dist=np.inf, cut_off_num=None, allowed_structure_site_structures=None,
                          keep_binding_vdW_bonds=True, keep_vdW_surface_bonds=False)
            site = slab_sites[0]
            add_adsorbate_to_site(adslab,adsorbate=adslab,surf_ind=0,site=site,height=gas_height)
            adslabmol = adslabmol.merge(mol_dict[species_names[0]])
            adslabmol.update_multiplicity()
            adslabmol.update_atomtypes()
            adslabmol.update_connectivity_values()
            adslabmol.identify_ring_membership()
        if len(adss) == 1:
            tsstructs.append(adslab)
            tsmols.append(adslabmol)
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
                occ = get_occupied_sites(adslab,slab_sites,nslab)
                for site in slab_sites:
                    adslab2 = adslab.copy()
                    adslabmol2 = adslabmol.copy(deep=True)
                    if not any(sites_match(site,osite,slab) for osite in occ) and site["site"] == sitetype1:
                        add_adsorbate_to_site(adslab2,adsorbate=ad,surf_ind=surf_ind1,site=site,height=height1)
                        adslabmol_ind2 = [i for i,s in enumerate(neighbor_sites) if sites_match(s,site,slab)][0]
                        adslabmol2 = add_coadsorbate_2D(adslabmol2,site,mol_dict[species_names[1]],slab,neighbor_sites,[adslabmol_ind2])
                        if len(adss) == 2:
                            tsstructs.append(adslab2)
                            tsmols.append(adslabmol2)
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
                                occ2 = get_occupied_sites(adslab2,slab_sites,nslab)
                                for site2 in slab_sites:
                                    adslab3 = adslab2.copy()
                                    adslabmol3 = adslabmol2.copy(deep=True)
                                    if not any(sites_match(site2,osite,slab) for osite in occ2) and site2["site"] == sitetype2:
                                        add_adsorbate_to_site(adslab3,adsorbate=ad2,surf_ind=surf_ind2,site=site2,height=height2)
                                        adslabmol_ind3 = [i for i,s in enumerate(neighbor_sites) if sites_match(s,site2,slab)][0]
                                        adslabmol3 = add_coadsorbate_2D(adslabmol3,site2,mol_dict[species_names[2]],slab,neighbor_sites,[adslabmol_ind3])
                                        if len(adss) == 3:
                                            tsstructs.append(adslab3)
                                            tsmols.append(adslabmol3)
                                        else:
                                            raise ValueError("Cannot handle more than three reactants")
                            else:
                                occl2 = get_occupied_sites(adslab2,slab_sites,nslab)
                                c = 0
                                sitel2 = slab_sites[c]
                                while any(sites_match(sitel2,osite,slab) for osite in occl2): #while occupied
                                    c += 1
                                    sitel2 = slab_sites[c]
                                add_adsorbate_to_site(adslab2,adsorbate=adss[2],surf_ind=0,site=sitel2,height=gas_height)
                                adslabmol3 = adslabmol2.merge(mol_dict[species_names[2]])
                                adslabmol3.update_multiplicity()
                                adslabmol3.update_atomtypes()
                                adslabmol3.update_connectivity_values()
                                adslabmol3.identify_ring_membership()
                                if len(adss) == 3:
                                    tsstructs.append(adslab2)
                                    tsmols.append(adslabmol3)
                                else:
                                    raise ValueError("Cannot handle more than three reactants")
            elif num_surf_sites[1] == 0:
                occl1 = get_occupied_sites(adslab,slab_sites,nslab)
                c = 0
                sitel1 = slab_sites[c]
                while any(sites_match(sitel1,osite,slab) for osite in occl1): #while occupied
                    c += 1
                    sitel1 = slab_sites[c]

                add_adsorbate_to_site(adslab,adsorbate=adss[1],surf_ind=0,site=sitel1,height=gas_height)
                adslabmol2 = adslabmol.merge(mol_dict[species_names[1]])
                adslabmol2.update_multiplicity()
                adslabmol2.update_atomtypes()
                adslabmol2.update_connectivity_values()
                adslabmol2.identify_ring_membership()
                if len(adss) == 2:
                    tsstructs.append(adslab)
                    tsmols.append(adslabmol2)
                else:
                    c = 0
                    site2 = slab_sites[c]
                    while any(sites_match(site2,osite,slab) for osite in occl1) and site2 != sitel1:
                        c += 1
                        site2 = slab_sites[c]
                    add_adsorbate_to_site(adslab,adsorbate=adss[2],surf_ind=0,site=site2,height=gas_height)
                    adslabmol3 = adslabmol2.merge(mol_dict[species_names[2]])
                    adslabmol3.update_multiplicity()
                    adslabmol3.update_atomtypes()
                    adslabmol3.update_connectivity_values()
                    adslabmol3.identify_ring_membership()
                    if len(adss) == 3:
                        tsstructs.append(adslab)
                        tsmols.append(adslabmol3)
                    else:
                        raise ValueError("Cannot handle more than three reactants")

    unique_tsstructs = []
    unique_tsmols = []
    for i,m1 in enumerate(tsmols):
        for j,m2 in enumerate(unique_tsmols):
            if m1.is_isomorphic(m2,save_order=True):
                break 
        else:
            unique_tsmols.append(m1)
            unique_tsstructs.append(tsstructs[i])
    
    return unique_tsstructs,unique_tsmols,neighbor_sites,ninds

def get_unique_TS_templates_site_pairings(tsstructs,tsmols,forward_template,reverse_template,nsites,slab,neighbor_sites,ninds,slab_sites,nslab):
    broken_bonds,formed_bonds = get_broken_formed_bonds(forward_template,reverse_template)
    template = forward_template.copy(deep=True)
    for i,a in enumerate(template.atoms): #unlabel atoms that do not participate in reactions
        if a.label != "":
            for bd in list(broken_bonds)+list(formed_bonds):
                label1,label2 = list(bd)
                if label1 == a.label or label2 == a.label:
                    break
            else:
                a.label = ""
    
    empty_site_labels = [a.label for a in template.atoms if a.is_surface_site() and all(a2.is_surface_site() for a2 in a.bonds.keys())]
    
    unique_tsstructs = []
    unique_tsmols = []
    for i,tsmol in enumerate(tsmols):
        try:
            labeled_tsmols = get_labeled_full_TS_mol(template,tsmol)
        except IndexError:
            logging.warning("could not label a TS structure, usually means a reactant/product structure is bad")
            continue
        for labeled_tsmol in labeled_tsmols:
            unique_tsstructs.append(tsstructs[i].copy())
            unique_tsmols.append(labeled_tsmol)
    
    out_tsstructs = []
    out_tsmols = [] 
    out_target_sites = [] 
    label_site_mappings = []
    for i,tsstruct in enumerate(unique_tsstructs):
        occ = get_occupied_sites(tsstruct,neighbor_sites,nslab)
        unocc = [site for site in neighbor_sites if not any(sites_match(site,osite,slab) for osite in occ)]
        pd = [unocc for i in range(nsites)]
        for sites in itertools.product(*pd):
            unique_sites = []
            for site in sites:
                for site2 in unique_sites:
                    if site == site2:
                        break
                else:
                    unique_sites.append(site)
            if len(unique_sites) < len(sites): #no duplicate sites
                continue
            label_site_mapping = dict()
            tsmol = unique_tsmols[i].copy(deep=True)
            for j,s in enumerate(sites):
                ind = [i for i,x in enumerate(neighbor_sites) if sites_match(s,x,slab)][0]
                label = empty_site_labels[j]
                tsmol.atoms[ind].label = label
                label_site_mapping[label] = s
            
            for bd in list(broken_bonds)+list(formed_bonds): #create reaction bonds in tsmol
                label1,label2 = list(bd)
                if label1 == "" or label2 == "":
                    continue 
                else:
                    a1 = tsmol.get_labeled_atoms(label1)[0]
                    a2 = tsmol.get_labeled_atoms(label2)[0]
                    if tsmol.has_bond(a1,a2):
                        bdts = tsmol.get_bond(a1,a2)
                        bdts.set_order_str("R")
                    else:
                        tsmol.add_bond(Bond(a1,a2,order="R"))
            
            for bd in tsmol.get_all_edges(): #Our TS representation isn't entirely unique, but the representation fix_bond_orders generates is consistent and what we train the SIDT on
                if bd.is_reaction_bond() or (bd.atom1.is_surface_site() and bd.atom2.is_surface_site()):
                    continue 
                elif bd.atom1.is_surface_site() or bd.atom2.is_surface_site(): #surface to adsorbate bond
                    bd.set_order_str("vdW")
                else:
                    bd.set_order_str("S")
                
            fix_bond_orders(tsmol,allow_failure=True,keep_binding_vdW_bonds=True,keep_vdW_surface_bonds=True)
            tsmol.update_multiplicity()
            tsmol.update_atomtypes()
            tsmol.update_connectivity_values()
            tsmol.identify_ring_membership()
            out_tsstructs.append(tsstruct.copy())
            out_tsmols.append(tsmol)
            out_target_sites.append(sites)
            label_site_mappings.append(label_site_mapping)
    
    unique_out_tsstructs = []
    unique_out_tsmols = [] 
    unique_out_target_sites = [] 
    unique_label_site_mappings = []
    for i,m in enumerate(out_tsmols):
        for j,m2 in enumerate(unique_out_tsmols):
            if m.is_isomorphic(m2,save_order=True):
                break 
        else:
            unique_out_tsstructs.append(out_tsstructs[i])
            unique_out_tsmols.append(out_tsmols[i])
            unique_out_target_sites.append(out_target_sites[i])
            unique_label_site_mappings.append(label_site_mappings[i])
    
    return unique_out_tsstructs,unique_out_tsmols,unique_out_target_sites,unique_label_site_mappings
                
        
    
def generate_constraints_harmonic_parameters(tsstructs,tsmols,label_site_mappings,adsorbates,slab,forward_template,
                                             reverse_template,template_name,template_reversed,
                                            ordered_names,reverse_names,mol_dict,gratom_to_molecule_atom_maps,
                                            gratom_to_molecule_surface_atom_maps,nslab,facet,metal,slab_sites,site_adjacency):
    """
    Generate constraints and harmonic parameters for the harmonically forced optimization
    """
    constraint_lists = []
    atom_bond_potential_lists = []
    site_bond_potential_lists = []
    out_structs = []
    
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

    for i,tsstruct in enumerate(tsstructs):
        tsstruct_valid = True
        atom_bond_potentials = []
        site_bond_potentials = []
        fixed_bond_pairs = []
        tsmol = tsmols[i]
        label_site_mapping = label_site_mappings[i]
        
        occ = get_occupied_sites(tsstruct,slab_sites,nslab)
        occ_atom_inds = [site["bonding_index"] for site in occ]
        
        #if we don't identify all the template sites in the occupied site list skip this structure
        if not ({get_ase_index(forward_template.atoms.index(a),template_mol_map,molecule_to_gratom_maps,nslab,ads_sizes) for a in forward_template.get_adatoms()} <= set(occ_atom_inds)): 
            print("occupied sites and forward_template sites do not match")
            print("occ={0}".format(len(occ)))
            print("forward_template_occs={0}".format(len(forward_template.get_adatoms())))
            continue
        
        occ_bd_lengths = {site["bonding_index"]:site['bond_length'] for site in occ}
        occ_site_types = {site["bonding_index"]:site['site'] for site in occ}
        occ_site_pos = {site["bonding_index"]:site['position'] for site in occ}

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
                deq,k = estimate_deq_k(sidt,labels,dwell,tsmol)
                if pos is not None:
                    d = {"ind":ind,"site_pos":pos,"k":k,"deq":deq}
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
                    deq,k = estimate_deq_k_fixed_surf_bond(labels,dwell,tsmol)
                    d = {"ind":ind,"site_pos":pos,"k":k,"deq":deq}
                    site_bond_potentials.append(d)
                else:
                    dwell = occ_bd_lengths[ase_ind1]
                    sitetype = occ_site_types[ase_ind1]
                    pos = deepcopy(occ_site_pos[ase_ind1])
                    pos[2] += dwell #shift position up to atom position
                    ind = ase_ind1
                    deq,k = estimate_deq_k_fixed_surf_bond(labels,dwell,tsmol)
                    d = {"ind":ind,"site_pos":pos,"k":k,"deq":deq}
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
                deq,k = estimate_deq_k(sidt,labels,dwell,tsmol)
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
                
                site = label_site_mapping[[x for x in labels if tsmol.get_labeled_atoms(x)[0].is_surface_site()][0]]
                for (mol_ind,sitetype) in site_dict.keys():
                    if mol_ind == ind1r_mol and site["site"] == sitetype:
                        dwell = site_dict[(mol_ind,sitetype)]
                        deq,k = estimate_deq_k(sidt,labels,dwell,tsmol)
                        if reusing_site: #attaching to a site that was bonded to another atom
                            if reused_site_type == sitetype:
                                pos = deepcopy(reused_site_pos)
                                pos[2] += dwell
                                d = {"ind":ase_ind2,"site_pos":pos,"k":k,"deq":deq}
                                site_bond_potentials.append(d)
                        else: #attaching to a new empty site
                            pos = deepcopy(site["position"])
                            pos[2] += dwell
                            d = {"ind":ase_ind2,"site_pos":pos,"k":k,"deq":deq}
                            site_bond_potentials.append(d)
                        break
                else:
                    tsstruct_valid = False 
                    break

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
                
                site = label_site_mapping[[x for x in labels if tsmol.get_labeled_atoms(x)[0].is_surface_site()][0]]
                for (mol_ind,sitetype) in site_dict.keys():
                    if mol_ind == ind2r_mol and site["site"] == sitetype:
                        dwell = site_dict[(mol_ind,sitetype)]
                        deq,k = estimate_deq_k(sidt,labels,dwell,tsmol)
                        if reusing_site: #attaching to a site that was bonded to another atom
                            if reused_site_type == sitetype:
                                pos = deepcopy(reused_site_pos)
                                pos[2] += dwell
                                d = {"ind":ase_ind1,"site_pos":pos,"k":k,"deq":deq}
                                site_bond_potentials.append(d)
                        else:
                            pos = deepcopy(site["position"])
                            pos[2] += dwell
                            d = {"ind":ase_ind1,"site_pos":pos,"k":k,"deq":deq}
                            site_bond_potentials.append(d)
                        break
                else:
                    tsstruct_valid = False 
                    break

        if tsstruct_valid:
            if fixed_bond_pairs:
                constraint_list = [{"type": "FixBondLength", "a1": pair[0], "a2": pair[1]} for pair in fixed_bond_pairs]+["freeze slab"]
                constraint_lists.append(constraint_list)
            else:
                constraint_lists.append(["freeze slab"])
            atom_bond_potential_lists.append(atom_bond_potentials)
            site_bond_potential_lists.append(site_bond_potentials)
            out_structs.append(tsstruct)

    return out_structs,constraint_lists,atom_bond_potential_lists,site_bond_potential_lists

def estimate_deq_k(sidt,labels,dwell,tsmol):
    """
    Estimate the equilibrium bond length and force constant for broken/forming bonds
    0--(1--2)--3
    """
    mol = tsmol.copy(deep=True)
    
    surface_bond = False
    for a in mol.atoms:
        if a.label == "":
            continue 
        elif a.label in labels:
            a.label = "*"
            if a.is_surface_site():
                surface_bond = True
        else:
            a.label = ""
    
    logv,tr = sidt.evaluate(mol,trace=True)
    
    if tr == "Root":
        logging.warning("Bond factor estimation used Root group in SIDT")
    
    v = np.exp(logv)
    
    if surface_bond:
        v = max(v,1.0) #don't allow surface bond stretches to be smaller than 1.0
        return dwell*(v-1.0),100.0
    else:
        v = max(v,1.1) #require covalent bond stretches to be at least 1.1
        return dwell*v,100.0

def estimate_deq_k_fixed_surf_bond(labels,dwell,tsmol):
    """
    Estimate the equilibrium bond length and force constant for surface bonds
    """
    return 0.0,100.0

def get_surface_forming_bond_pairings(tsstructs,slab,atom_bond_potential_lists,
                                      site_bond_potential_lists,constraints_list,site_bond_dict_list,
                                      slab_sites):
    """
    Identify unique site targets for the harmonically forced optimization
    """
    if len(site_bond_dict_list[0]) == 0: #no site bonds formed
        return tsstructs,atom_bond_potential_lists,site_bond_potential_lists,constraints_list

    nslab = len(slab)
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
        occ = get_occupied_sites(tsstruct,slab_sites,nslab)
        unocc = [site for site in slab_sites if not any(sites_match(site,osite,slab) for osite in occ)]
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
                    params["site_pos"] = deepcopy(site["position"])
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
                params["site_pos"] = deepcopy(site["position"])
                params["ind"] = atm_ind
                params_list.append(params)
        site_bond_potentials[atm_ind] = params_list

    return site_bond_potentials
