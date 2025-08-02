from molecule.molecule import Molecule,Atom,Bond,Group,GroupAtom,GroupBond
from molecule.exceptions import AtomTypeError
from ase.io import read, write
from ase.geometry import get_distances
from ase.visualize import view
from ase.neighborlist import natural_cutoffs
from acat.adsorption_sites import SlabAdsorptionSites
from pynta.utils import get_unique_sym, get_occupied_sites, sites_match, SiteOccupationException
from pynta.mol import *
from pynta.geometricanalysis import *
from pysidt import *
from pysidt.extensions import split_mols
from pysidt.sidt import *
from fireworks import LaunchPad, Workflow
from fireworks.queue.queue_launcher import rapidfire as rapidfirequeue
from fireworks.features.multi_launcher import launch_multiprocess
from fireworks.utilities.fw_serializers import load_object_from_file
from fireworks.core.rocket_launcher import rapidfire
from fireworks.core.fworker import FWorker
from copy import deepcopy
import numpy as np
import json
import shutil
import os
import itertools 
import logging

def get_unstable_pairs(pairsdir,adsorbate_dir,sites,site_adjacency,nslab,max_dist=3.0,show=False,metal=None,facet=None,
                       infopath_dict=None,imag_freq_path_dict=None):
    out_pairs = []
    allowed_structure_site_structure_map = generate_allowed_structure_site_structures(adsorbate_dir,sites,site_adjacency,nslab,max_dist=max_dist,cut_off_num=None)
    if show:
        config_show = []
    for pair in os.listdir(pairsdir):
        if not "_" in pair or pair[0] == ".":
            continue
        adname = pair.split("_")[0]
        coadname = pair.split("_")[1]
        info_coad_path = os.path.join(adsorbate_dir,coadname,"info.json")
        with open(info_coad_path,'r') as f:
            info_coad = json.load(f)
        coad_struct = Molecule().from_adjacency_list(info_coad["adjlist"])
        cut_off_num = len([at for at in coad_struct.atoms if not at.is_surface_site()])
        if infopath_dict:
            info_path = infopath_dict[adname]
        else:
            info_path = None
        if imag_freq_path_dict:
            imag_freq_path = imag_freq_path_dict[adname]
        else:
            imag_freq_path = None
        for num in os.listdir(os.path.join(pairsdir,pair)):
            p = os.path.join(pairsdir,pair,num) 
            if num.isdigit() and os.path.isdir(p):
                init = read(os.path.join(p,"init.xyz"))
                with open(os.path.join(p,"info.json"),'r') as f:
                    initinfo = json.load(f)
                try:
                    m = Molecule().from_adjacency_list(initinfo["adjlist"],check_consistency=False)
                    #m.update(sort_atoms=False)
                except AtomTypeError:
                    m = Molecule().from_adjacency_list(initinfo["adjlist"],raise_atomtype_exception=False,
                            raise_charge_exception=False,check_consistency=False)
                    for a in m.atoms:
                        if a.charge != 0:
                            a.charge = 0
                            v = 0
                            for a2,edge in a.edges.items():
                                v += edge.order*2
                            if not a.is_hydrogen():
                                a.lone_pairs = int(np.round((8-v)/2))
#                     try:
#                         m.update(sort_atoms=False)
#                     except Exception as e:
#                         print(m.to_adjacency_list())
#                         raise e
                is_ts = any(bd.get_order_str() == "R" for bd in m.get_all_edges())
                g = reduce_graph_to_pairs(m)
                
                #g.update(sort_atoms=False)
                outpath = os.path.join(p,"out.xyz")
                if not os.path.exists(outpath):
                    if show:
                        pass
                        
                    out_pairs.append(g.to_group())
                else:
                    final = read(outpath)
                    try:
                        gout = extract_pair_graph(final,sites,site_adjacency,nslab,max_dist=max_dist,
                                                  allowed_structure_site_structures=allowed_structure_site_structure_map,is_ts=is_ts,
                                                  metal=metal,facet=facet,info_path=info_path,imag_freq_path=imag_freq_path)
                        if len(gout.atoms) != len(g.atoms):
                            out_pairs.append(g.to_group())
                            if show:
                                config_show.append(init)
                                config_show.append(final)
                            continue
                        #gout.update(sort_atoms=False)
                    except (FindingPathError, SiteOccupationException, TooManyElectronsException):
                        continue
                    except FailedFixBondsException: #Pynta is unable to understand the final structure in a resonance sense...definitely not isomorphic
                        logging.error("Failed to fix bonds for structure: {}".format(p))
                        out_pairs.append(g.to_group())
                        if show:
                            config_show.append(init)
                            config_show.append(final)
                        continue
                    
                    if not gout.is_isomorphic(g,save_order=True,strict=False):
                        out_pairs.append(g.to_group())
                        if show:
                            config_show.append(init)
                            config_show.append(final)
    if show:
        view(config_show)
    return out_pairs

def copy_stable_pairs(pairsdir,sites,site_adjacency,nslab,max_dist=3.0):
    out_pairs = []
    count = 0
    coadname_dict = {"O=[Pt]": 1, "N#[Pt]": 1, "O[Pt]": 2, "[Pt]": 1}
    for pair in os.listdir(pairsdir):
        coadname = pair.split("_")[1]
        for num in os.listdir(os.path.join(pairsdir,pair)):
            p = os.path.join(pairsdir,pair,num) 
            if num.isdigit() and os.path.isdir(p):
                init = read(os.path.join(p,"init.xyz"))
                with open(os.path.join(p,"info.json"),'r') as f:
                    initinfo = json.load(f)
                try:
                    m = Molecule().from_adjacency_list(initinfo["adjlist"])
                    #m.update(sort_atoms=False)
                except AtomTypeError:
                    m = Molecule().from_adjacency_list(initinfo["adjlist"],raise_atomtype_exception=False,
                            raise_charge_exception=False,check_consistency=False)
                    for a in m.atoms:
                        if a.charge != 0:
                            a.charge = 0
                            v = 0
                            for a2,edge in a.edges.items():
                                v += edge.order*2
                            if not a.is_hydrogen():
                                a.lone_pairs = int(np.round((8-v)/2))
#                     try:
#                         m.update(sort_atoms=False)
#                     except Exception as e:
#                         print(m.to_adjacency_list())
#                         raise e
                
                g = extract_pair_graph(init,sites,site_adjacency,nslab,max_dist=max_dist,
                                       cut_coads_off_num=coadname_dict[coadname])
                    
                
                #g.update(sort_atoms=False)
                outpath = os.path.join(p,"out.xyz")
                if not os.path.exists(outpath):
                    out_pairs.append(g.to_group())
                else:
                    final = read(outpath)
                    try:
                        gout = extract_pair_graph(final,sites,site_adjacency,nslab,max_dist=max_dist)
                        if len(gout.atoms) != len(g.atoms):
                            out_pairs.append(g.to_group())
                            continue
                        #gout.update(sort_atoms=False)
                    except FindingPathError:
                        continue
                    
                    if not gout.is_isomorphic(g,save_order=True,strict=False):
                        out_pairs.append(g.to_group())
                    else: #matches
                        os.makedirs(os.path.join(os.path.split(pairsdir)[0],"sampling",str(count)))
                        shutil.copyfile(os.path.join(p,"out.xyz"),
                                        os.path.join(os.path.split(pairsdir)[0],"sampling",str(count),"out.xyz")
                                       )
                        shutil.copyfile(os.path.join(p,"info.json"),
                                        os.path.join(os.path.split(pairsdir)[0],"sampling",str(count),"info.json")
                                       )
                        shutil.copyfile(os.path.join(p,"init.xyz"),
                                        os.path.join(os.path.split(pairsdir)[0],"sampling",str(count),"init.xyz")
                                       )
                        count += 1


    return None

def generate_pair_geometries(adpath1,adpath2,slabpath,metal,facet,sites,site_adjacency,adinfo1=None,adinfo2=None,
                             max_dist=3.0,imag_freq_max=150.0,symmetric=None):
    """
    adpath1 can be bidentate
    adpath2 must be monodentate
    """
    #slab information
    slab = read(slabpath)
    nslab = len(slab)
    
    adsorbate_dir = os.path.split(adpath2)[0]
    allowed_structure_site_structures = generate_allowed_structure_site_structures(adsorbate_dir,sites,site_adjacency,nslab,max_dist=max_dist)

    if symmetric is None:
        symmetric = adpath1 == adpath2
    
    #extract information about adsorbates and valid adsorbate geometries
    if os.path.isfile(os.path.join(adpath1,"opt.xyz")): #TS
        is_ts = True
        if adinfo1 is None:
            adinfo1 = os.path.join(os.path.split(adpath1)[0],"info.json")
        
        with open(adinfo1,"r") as f:
            info1 = json.load(f)

        reactants = Molecule().from_adjacency_list(info1["reactants"])
        products = Molecule().from_adjacency_list(info1["products"])
        keep_binding_vdW_bonds_in_reactants=False
        keep_vdW_surface_bonds_in_reactants=False
        mol = reactants
        for bd in mol.get_all_edges():
            if bd.order == 0:
                if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
                    keep_binding_vdW_bonds_in_reactants = True
                    m = mol.copy(deep=True)
                    b = m.get_bond(m.atoms[mol.atoms.index(bd.atom1)],m.atoms[mol.atoms.index(bd.atom2)])
                    m.remove_bond(b)
                    out = m.split()
                    if len(out) == 1: #vdW bond is not only thing connecting adsorbate to surface
                        keep_vdW_surface_bonds_in_reactants = True
        keep_binding_vdW_bonds_in_products=False
        keep_vdW_surface_bonds_in_products=False
        mol = products
        for bd in mol.get_all_edges():
            if bd.order == 0:
                if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
                    keep_binding_vdW_bonds_in_products = True
                    m = mol.copy(deep=True)
                    b = m.get_bond(m.atoms[mol.atoms.index(bd.atom1)],m.atoms[mol.atoms.index(bd.atom2)])
                    m.remove_bond(b)
                    out = m.split()
                    if len(out) == 1: #vdW bond is not only thing connecting adsorbate to surface
                        keep_vdW_surface_bonds_in_products = True
        
        keep_binding_vdW_bonds = keep_binding_vdW_bonds_in_reactants and keep_binding_vdW_bonds_in_products
        keep_vdW_surface_bonds = keep_vdW_surface_bonds_in_reactants and keep_vdW_surface_bonds_in_products
        
        admol1,neighbor_sites1,ninds1 = generate_TS_2D(read(os.path.join(adpath1,"opt.xyz")), adinfo1, metal, facet, sites, site_adjacency, nslab, max_dist=np.inf,
                                                    imag_freq_path=os.path.join(adpath1,"vib.json_vib.json"),
                                                    allowed_structure_site_structures=allowed_structure_site_structures,
                                                    keep_binding_vdW_bonds=keep_binding_vdW_bonds,
                                                    keep_vdW_surface_bonds=keep_vdW_surface_bonds,
        )
        
        aseinds1 = []
        for i,at in enumerate(admol1.atoms):
            if (not at.is_surface_site()) and at.is_bonded_to_surface():
                aseinds1.append(i-len(neighbor_sites1)+nslab)
        ad1s = [read(os.path.join(adpath1,"opt.xyz"))]
        ad1xyzs = [os.path.join(adpath1,"opt.xyz")]
        ad12Ds = [admol1]
        ad12Dneighbors = [neighbor_sites1]
        ad12Dninds = [ninds1]
    else:
        is_ts = False
        if adinfo1 is None:
            adinfo1 = os.path.join(adpath1,"info.json")
        
        with open(adinfo1,"r") as f:
            info1 = json.load(f)
        
        keep_binding_vdW_bonds = False
        keep_vdW_surface_bonds = False
        mol1 = Molecule().from_adjacency_list(info1["adjlist"])
        for bd in mol1.get_all_edges():
            if bd.order == 0:
                if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
                    keep_binding_vdW_bonds = True
                    m = mol1.copy(deep=True)
                    b = m.get_bond(m.atoms[mol1.atoms.index(bd.atom1)],m.atoms[mol1.atoms.index(bd.atom2)])
                    m.remove_bond(b)
                    out = m.split()
                    if len(out) == 1: #vdW bond is not only thing connecting adsorbate to surface
                        keep_vdW_surface_bonds = True
                            
        atom_to_molecule_surface_atom_map1 = {int(key):int(val) for key,val in info1["gratom_to_molecule_surface_atom_map"].items()}
        ad1s,ad1xyzs = get_unique_adsorbate_geometries(adpath1,mol1,sites,site_adjacency,atom_to_molecule_surface_atom_map1,
                                    nslab,imag_freq_max=imag_freq_max,return_xyzs=True)
        ad12Ds = []
        ad12Dneighbors = []
        ad12Dninds = []
        for a in ad1s:
            admol1,neighbor_sites1,ninds1 = generate_adsorbate_2D(a, sites, site_adjacency, nslab, max_dist=np.inf,keep_binding_vdW_bonds=keep_binding_vdW_bonds,
                                                                  keep_vdW_surface_bonds=keep_vdW_surface_bonds)
            ad12Ds.append(admol1)
            ad12Dneighbors.append(neighbor_sites1)
            ad12Dninds.append(ninds1)
        
        aseinds1 = [x+nslab for x in atom_to_molecule_surface_atom_map1.keys()]
        
    if os.path.isfile(adpath2):
        raise ValueError
    else:
        if adinfo2 is None:
            adinfo2 = os.path.join(adpath2,"info.json")
        
        with open(adinfo2,"r") as f:
            info2 = json.load(f)
        
        mol2 = Molecule().from_adjacency_list(info2["adjlist"])
        assert len(mol2.get_adatoms()) == 1
        atom_to_molecule_surface_atom_map2 = { int(key):int(val) for key,val in info2["gratom_to_molecule_surface_atom_map"].items()}
        ad2s = get_unique_adsorbate_geometries(adpath2,mol2,sites,site_adjacency,atom_to_molecule_surface_atom_map2,
                                    nslab,imag_freq_max=imag_freq_max)
    

    #generate pairs
    pairs = []
    ad1_to_ad1_sites = dict()
    ad1_to_ad2_sites = dict()
    ad1_to_ad2_heights = dict()
    ad1_to_actual_occ = dict()
    #go through all adsorbate 1 configurations
    for j,ad1 in enumerate(ad1s):
        #find occupied sites
        ad1_to_ad2_sites[j] = dict()
        ad1_to_ad2_heights[j] = dict()
        occ = get_occupied_sites(ad1,sites,nslab)
        actual_occ = []
        for site in occ:
            if site["bonding_index"] in aseinds1:
                actual_occ.append(site)
        
        if len(actual_occ) == 0:
#             print(adinfo1)
#             print(aseinds1)
#             print(atom_to_molecule_surface_atom_map1)
#             print(mol1.to_adjacency_list())
#             print(occ)
#             print(j)
            raise ValueError
        
        ad1_to_ad1_sites[j] = [x for x in actual_occ]
        
        #find all sites to be resolved with pairs
        neighbor_sites = []
        
        for site in sites:
            if any(sites_match(site,s,slab) for s in actual_occ):
                continue
            for occ_site in actual_occ:
                bd,dist = get_distances([site["position"]], [occ_site["position"]], cell=slab.cell, pbc=slab.pbc)
                v = np.linalg.norm(bd[:2])
                if v < max_dist:
                    neighbor_sites.append(site)
                    break
        
        #estimate heights for all placements
        stable_neighbor_sites = []
        heights = []
        for site in neighbor_sites:
            for i,ad2 in enumerate(ad2s):
                if i not in ad1_to_ad2_sites[j].keys():
                    ad1_to_ad2_sites[j][i] = []
                    ad1_to_ad2_heights[j][i] = []
        
                occ = get_occupied_sites(ad2,sites,nslab)
                occsite = None
                
                for s in occ:
                    if (s['bonding_index'] - nslab) in atom_to_molecule_surface_atom_map2.keys():
                        occsite = s
                        break
                else:
#                     print(occ)
#                     print(atom_to_molecule_surface_atom_map2.keys())
#                     print(s['bonding_index'])
                    raise ValueError

                if occsite["site"] == site["site"] and occsite["morphology"] == site["morphology"]:
                    ad1_to_ad2_sites[j][i].append(site)
                    ad1_to_ad2_heights[j][i].append(occsite["bond_length"])
                    stable_neighbor_sites.append(site)
                    heights.append(occsite["bond_length"])
                    break

    
    adpairs = []
    pairmols = []
    adpairs_orig_xyzs = []
    
    if not symmetric:
        for j in range(len(ad1s)):
            ad1_sites = ad1_to_ad1_sites[j]
            ad1xyz = ad1xyzs[j]
            ad2_sites = []
            ad2_geoms = []
            heights = []
            for i,sites in ad1_to_ad2_sites[j].items():
                for k,site in enumerate(sites):
                    heights.append(ad1_to_ad2_heights[j][i][k])
                    ad2_sites.append(ad1_to_ad2_sites[j][i][k])
                    ad2_geoms.append(ad2s[i])


            inds = get_unique_site_inds(ad2_sites,slab,fixed_points=[x["position"] for x in ad1_sites])

            for i in inds:
#                 if any(sites_match(ad2_sites[i],s,slab) for s in ad1_to_ad1_sites[j]):
#                     continue
                
                ad = deepcopy(ad1s[j])
                surf_ind = list(atom_to_molecule_surface_atom_map2.keys())[0]
                add_adsorbate_to_site(ad, ad2_geoms[i][nslab:], surf_ind, ad2_sites[i], height=heights[i])
                nsites = ad12Dneighbors[j]
                ind = [k for k,x in enumerate(nsites) if sites_match(ad2_sites[i],x,slab)][0]
                amol = deepcopy(ad12Ds[j])
                satom = amol.atoms[ind]
                m2 = deepcopy(mol2)
                admol2 = m2.get_desorbed_molecules()[0]
                label_atom_dict = admol2.get_all_labeled_atoms()
                for label,at in label_atom_dict.items(): #one iteration only
                    amol = amol.merge(admol2)
                    if label == "*1":
                        at.decrement_radical()
                        order = 1
                    elif label == "*2":
                        at.decrement_radical()
                        at.decrement_radical()
                        order = 2
                    elif label == "*3":
                        at.decrement_radical()
                        at.decrement_lone_pairs()
                        order = 3
                    elif label == "*4":
                        at.decrement_radical()
                        at.decrement_radical()
                        at.decrement_lone_pairs()
                        order = 4
                    else:
                        raise ValueError(label)
                    try:
                        bd = Bond(satom,at,order=order)
                    except:
#                         print(ad12Ds[j].to_adjacency_list())
#                         print(mol2.to_adjacency_list())
#                         print(amol.to_adjacency_list())
#                         print(satom)
#                         print(at)
#                         print(order)
                        raise ValueError
                    amol.add_bond(bd)
                    amol.clear_labeled_atoms()
                    if amol.multiplicity == -187: #handle surface molecules
                        amol.multiplicity = amol.get_radical_count() + 1
                for at in amol.atoms:
                    at.update_charge()
                for pmol in pairmols:
                    if pmol.is_isomorphic(amol,save_order=True): #duplicate
                        break
                else:
                    adpairs.append(ad)
                    pairmols.append(amol)
                    adpairs_orig_xyzs.append(ad1xyz)
    else: #symmetric case, monodentate-monodentate since ad2 must be monodentate
        ad2_site_pairs = []
        ad2_geoms = []
        heights = []
        ad1_inds = []
        for j in range(len(ad1s)):
            ad1_sites = ad1_to_ad1_sites[j]
            ad1xyz = ad1xyzs[j]
            for i,sites in ad1_to_ad2_sites[j].items():
                for k,site in enumerate(sites):
                    heights.append(ad1_to_ad2_heights[j][i][k])
                    ad2_site_pairs.append((ad1_sites[0],site))
                    ad2_geoms.append(ad2s[i])
                    ad1_inds.append(j)

        inds = get_unique_site_pair_inds(ad2_site_pairs,slab)

        for i in inds:
#             if any(sites_match(ad2_site_pairs[i][1],s,slab) for s in ad1_to_ad1_sites[j]):
#                 print("matched")
#                 continue
            ad = deepcopy(ad1s[ad1_inds[i]])
            surf_ind = list(atom_to_molecule_surface_atom_map2.keys())[0]
            add_adsorbate_to_site(ad, ad2_geoms[i][nslab:], surf_ind, ad2_site_pairs[i][1], 
                                  height=heights[i])
            nsites = ad12Dneighbors[ad1_inds[i]]
            try:
                ind = [k for k,x in enumerate(nsites) if sites_match(ad2_site_pairs[i][1],x,slab)][0]
            except Exception as e:
#                 print("site")
#                 print(ad2_site_pairs[i][1])
#                 print("occ")
#                 print(ad1_to_ad1_sites[j])
#                 print("nsites")
#                 print(nsites)
#                 bd,dist = get_distances([ad2_site_pairs[i][1]["position"]], [ad1_to_ad1_sites[j][0]["position"]], cell=slab.cell, pbc=slab.pbc)
#                 print(bd)
#                 print(dist)
                raise e
            amol = deepcopy(ad12Ds[j])
            satom = amol.atoms[ind]
            m2 = deepcopy(mol2)
            admol2 = m2.get_desorbed_molecules()[0]
            label_atom_dict = admol2.get_all_labeled_atoms()
            for label,at in label_atom_dict.items(): #one iteration only
                amol = amol.merge(admol2)
                if label == "*1":
                    at.decrement_radical()
                    order = 1
                elif label == "*2":
                    at.decrement_radical()
                    at.decrement_radical()
                    order = 2
                elif label == "*3":
                    at.decrement_radical()
                    at.decrement_lone_pairs()
                    order = 3
                elif label == "*4":
                    at.decrement_radical()
                    at.decrement_radical()
                    at.decrement_lone_pairs()
                    order = 4
                else:
                    raise ValueError(label)
                try:
                    bd = Bond(satom,at,order=order)
                except:
#                     print(ad12Ds[j].to_adjacency_list())
#                     print(mol2.to_adjacency_list())
#                     print(amol.to_adjacency_list())
#                     print(satom)
#                     print(at)
#                     print(order)
                    raise ValueError
                amol.add_bond(bd)
                amol.clear_labeled_atoms()
                if amol.multiplicity == -187: #handle surface molecules
                    amol.multiplicity = amol.get_radical_count() + 1
            for at in amol.atoms:
                at.update_charge()
            for pmol in pairmols:
                if pmol.is_isomorphic(amol,save_order=True): #duplicate
                    break
            else:
                adpairs.append(ad)
                pairmols.append(amol)
                adpairs_orig_xyzs.append(ad1xyz)

    return adpairs,pairmols,adpairs_orig_xyzs

def get_unique_site_inds(sites,slab,fixed_points=None,tol=0.15):
    fingerprints = []
    for k,site in enumerate(sites):
        if fixed_points is None:
            fingerprints.append((site["morphology"],site["site"]))
        else:
            dists = [site["morphology"],site["site"]]
            for fixed_point in fixed_points:
                bd,d = get_distances([site["position"]], [fixed_point], cell=slab.cell, pbc=(True,True,False))
                xydist = np.linalg.norm(bd[0][0][:2])
                zdist = bd[0][0][2]
                dists.extend([xydist,zdist])
            fingerprints.append(tuple(dists))
    
    unique_sites = []
    unique_inds = []
    for i,f in enumerate(fingerprints):
        boo = False
        for uf in unique_sites:  
            if fingerprints_match(f,uf,tol=tol):
                boo = True
                break
        if boo:
            continue
        else:
            unique_sites.append(f)
            unique_inds.append(i)

    return unique_inds

def fingerprints_match(f1,f2,tol=0.15):
    for i in range(len(f1)):
        if isinstance(f1[i],str) or isinstance(f1[i],frozenset):
            if f1[i] != f2[i]:
                return False
        elif isinstance(f1[i],float):
            if abs(f1[i] - f2[i]) > tol:
                return False
        else:
            raise ValueError
    else:
        return True

def get_unique_site_pair_inds(site_pairs,slab,tol=0.15):
    fingerprints = []
    for k,sites in enumerate(site_pairs):
        site1 = sites[0]
        site2 = sites[1]
        bd,d = get_distances([site1["position"]], [site2["position"]], cell=slab.cell, pbc=(True,True,False))
        xydist = np.linalg.norm(bd[0][0][:2])
        zdist = bd[0][0][2]
        fingerprints.append((frozenset([(site1["morphology"],site1["site"]),(site2["morphology"],site2["site"])]),
                             xydist,zdist,))
    
    unique_sites = []
    unique_inds = []
    for i,f in enumerate(fingerprints):
        boo = False
        for uf in unique_sites:  
            if fingerprints_match(f,uf,tol=tol):
                boo = True
                break
        if boo:
            continue
        else:
            unique_sites.append(f)
            unique_inds.append(i)

    return unique_inds

def setup_pair_opts_for_rxns(targetdir,adsorbates,tsdirs,coadnames,metal,facet,sites,site_adjacency,max_dist=3.0,imag_freq_max=150.0):
    pairdir = os.path.join(targetdir,"pairs")
    addir = os.path.join(os.path.split(os.path.split(tsdirs[0])[0])[0],"Adsorbates")
    slabpath = os.path.join(os.path.split(os.path.split(tsdirs[0])[0])[0],"slab.xyz")
    if not os.path.exists(pairdir):
        os.makedirs(pairdir)
    
    ads = set(adsorbates)
    combs = []
    for tsdir in tsdirs:
        for coadname in coadnames:
            tp = (tsdir,coadname)
            combs.append(tp)
        
        with open(os.path.join(os.path.split(tsdir)[0],"info.json"),"r") as f:
            info = json.load(f)
        for molname in info["species_names"]+info["reverse_names"]:
            with open(os.path.join(addir,molname,"info.json"),"r") as f:
                molinfo = json.load(f)
            m = Molecule().from_adjacency_list(molinfo["adjlist"])
            if m.contains_surface_site():
                ads.add(molname)
    
    
    for adname in ads:
        for coadname in coadnames:
            tp = (adname,coadname)
            revtp = (coadname,adname)
            if (revtp not in combs) and (tp not in combs):
                combs.append(tp)
                
    outdirs_ad = []
    outdirs_ts = []
    for s in combs:
        if os.path.exists(s[0]):
            is_ts = True
        else: 
            is_ts = False
        if not is_ts:
            name = "_".join(s)
        else:
            name = "_".join([os.path.split(os.path.split(s[0])[0])[1],s[1]])
        namedir = os.path.join(pairdir,name)
        if not os.path.exists(namedir):
            os.makedirs(namedir)
            if not is_ts:
                ds = [os.path.join(addir,x) for x in s]
            else:
                ds = [s[0],os.path.join(addir,s[1])]
            pairs,pairmols,pairxyzs = generate_pair_geometries(ds[0],ds[1],slabpath,metal,facet,sites,site_adjacency,
                                 max_dist=max_dist,imag_freq_max=imag_freq_max)
            for i,pair in enumerate(pairs):
                os.makedirs(os.path.join(namedir,str(i)))
                write(os.path.join(namedir,str(i),"init.xyz"), pair)
                if not is_ts:
                    moldict = {"adjlist": pairmols[i].to_adjacency_list(),"xyz": pairxyzs[i], "coadname": s[1]}
                else:
                    moldict = {"adjlist": pairmols[i].to_adjacency_list(),"tsdir": s[0], "xyz": pairxyzs[i], "coadname": s[1]}
                with open(os.path.join(namedir,str(i),"info.json"),'w') as f:
                    json.dump(moldict,f)
                if not is_ts:
                    outdirs_ad.append(os.path.join(namedir,str(i),"init.xyz"))
                else:
                    outdirs_ts.append(os.path.join(namedir,str(i),"init.xyz"))
                    
    return outdirs_ad,outdirs_ts

def get_adsorbate_ts_information(xyz,slabxyz,is_ts,
                               metal,facet,sites,site_adjacency,max_dist=3.0):
    slab = read(slabxyz)
    
        
    if is_ts:
        ad = read(xyz)
        xyzinfo = os.path.join(os.path.split(xyz)[0],"..","info.json")
        with open(xyzinfo,"r") as f:
            info = json.load(f)

        reactants = Molecule().from_adjacency_list(info["reactants"])
        products = Molecule().from_adjacency_list(info["products"])
        template_mol_map = [{ int(key):int(val) for key,val in x.items()} for x in info["template_mol_map"]]
        molecule_to_atom_maps = [{ int(key):int(val) for key,val in x.items()} for x in info["molecule_to_atom_maps"]]
        ads_sizes = info["ads_sizes"]
        nslab = info["nslab"]
        forward = info["forward"]

        broken_bonds,formed_bonds = get_broken_formed_bonds(reactants,products)

        extra_bonds = formed_bonds if forward else broken_bonds
        template = reactants if forward else products
        adatoms = template.get_adatoms()
        adinds = []
        for ind,atom in enumerate(template.atoms):
            if atom.is_surface_site():
                if len(atom.bonds) == 0:
                    s = [bd for bd in extra_bonds if atom.label in bd]
                    if len(s) > 0:
                        labels = list(s[0])
                        labels.remove(atom.label)
                        alabel = labels[0]
                        a = template.get_labeled_atoms(alabel)[0]
                        adinds.append(template.atoms.index(a))
                else:
                    a = list(atom.bonds.keys())[0]
                    adinds.append(template.atoms.index(a))

        aseinds = []
        for ind in adinds:
            aseind = get_ase_index(ind,template_mol_map,molecule_to_atom_maps,nslab,ads_sizes)
            aseinds.append(aseind)
            
        
    else:
        ad = read(xyz)
        with open(os.path.join(os.path.split(os.path.split(xyz)[0])[0],"info.json"),"r") as f:
            info = json.load(f)
        
        mol2D = Molecule().from_adjacency_list(info["adjlist"])
        atom_to_molecule_surface_map = { int(key):int(val) for key,val in info["gratom_to_molecule_surface_atom_map"].items()}
        nslab = info["nslab"]
        aseinds = [x+nslab for x in atom_to_molecule_surface_map.keys()]
        
    occ = get_occupied_sites(ad,sites,nslab)
    
    
    actual_occ = []
    for site in occ:
        if site["bonding_index"] in aseinds:
            actual_occ.append(site)

    neighbor_sites = []

    for site in sites:
        if any(sites_match(site,s,slab) for s in actual_occ):
            continue
        for occ_site in actual_occ:
            v,dist = get_distances([site["position"]], [occ_site["position"]], cell=slab.cell, pbc=slab.pbc)
            if np.linalg.norm(v[:2]) < max_dist:
                neighbor_sites.append(site)
                break
    
    if is_ts:
        admol,neighbor_sites_2D,ninds = generate_TS_2D(ad, xyzinfo, metal, facet, sites, site_adjacency, nslab, max_dist=None)
    else:
        admol,neighbor_sites_2D,ninds = generate_adsorbate_2D(ad, sites, site_adjacency, nslab, max_dist=None)
    
    return ad,admol,neighbor_sites_2D,ninds,actual_occ,neighbor_sites,aseinds,slab,nslab

def get_coadsorbate_information(coadnames,ads_dir,neighbor_sites,sites,site_adjacency,nslab,slab,imag_freq_max=150.0):
    coad_to_stable_neighbor_sites = dict()
    coad_site_stable_parameters = {name:[] for name in coadnames}
    coad_atom_to_molecule_surface_atom_map = dict()
    infocoad_dict = dict()
    stable_neighbor_sites_total = []
    coads_dict = dict()
    coad2Ds = dict()
    for coadname in coadnames:
        coaddir = os.path.join(ads_dir,coadname)
        with open(os.path.join(coaddir,"info.json"),"r") as f:
            infocoad = json.load(f)

        coad2D = Molecule().from_adjacency_list(infocoad["adjlist"])
        coad2Ds[coadname] = coad2D
        atom_to_molecule_surface_atom_map = {int(key):int(val) for key,val in infocoad["gratom_to_molecule_surface_atom_map"].items()}
        coads = get_unique_adsorbate_geometries(coaddir,Molecule().from_adjacency_list(infocoad["adjlist"]),
                               sites,site_adjacency,atom_to_molecule_surface_atom_map,
                               nslab,imag_freq_max=imag_freq_max)
        atom_to_molecule_surface_atom_map = {int(key):int(val) for key,val in infocoad["gratom_to_molecule_surface_atom_map"].items()}
        coad_atom_to_molecule_surface_atom_map[coadname] = atom_to_molecule_surface_atom_map
        coads_dict[coadname] = coads
        infocoad_dict[coadname] = infocoad
        stable_neighbor_sites = []

        for site in neighbor_sites:

            for coad in coads:
                occ = get_occupied_sites(coad,sites,nslab)
                occsite = None
                for s in occ:
                    if (s['bonding_index'] - nslab) in atom_to_molecule_surface_atom_map.keys():
                        occsite = s
                        break
                else:
                    raise ValueError

                if occsite["site"] == site["site"] and occsite["morphology"] == site["morphology"]:
                    stable_neighbor_sites.append(site)
                    if not any(sites_match(s,site,slab) for s in stable_neighbor_sites_total):
                        stable_neighbor_sites_total.append(site)
                    if not (site["site"],site["morphology"]) in coad_site_stable_parameters[coadname]:
                        coad_site_stable_parameters[coadname].append((site["site"],site["morphology"]))
                    break

        coad_to_stable_neighbor_sites[coadname] = stable_neighbor_sites

    return coad_to_stable_neighbor_sites, coad_site_stable_parameters,coad_atom_to_molecule_surface_atom_map,infocoad_dict,stable_neighbor_sites_total,coads_dict,coad2Ds

def generate_coadsorbed_xyzs(outdir,ad_xyzs,ts_xyzs,slabxyz,pairsdir,ads_dir,
                             coadnames,metal,facet,sites,site_adjacency,max_dist=3.0):
    slab = read(slabxyz)
    nslab = len(slab)

    unstable_pairs = get_unstable_pairs(pairsdir,ads_dir,sites,site_adjacency,nslab,max_dist=None)
    
    ad_dict = dict()
    ad_admol_dict = dict()
    ad_neighbor_sites_2D_dict = dict()
    ad_ninds_dict = dict()
    ad_neighbor_sites_dict = dict()
    ad_actual_occ_dict = dict()
    ad_aseinds_dict = dict()
    
    for ad_xyz in ad_xyzs:
        ad,admol,neighbor_sites_2D,ninds,actual_occ,neighbor_sites,aseinds,slab,nslab = get_adsorbate_ts_information(ad_xyz,slabxyz,False,metal,facet,sites,site_adjacency,max_dist=max_dist)
        ad_dict[ad_xyz] = ad
        ad_admol_dict[ad_xyz] = admol
        ad_neighbor_sites_2D_dict[ad_xyz] = neighbor_sites_2D
        ad_ninds_dict[ad_xyz] = ninds
        ad_neighbor_sites_dict[ad_xyz] = neighbor_sites
        ad_actual_occ_dict[ad_xyz] = actual_occ
        ad_aseinds_dict[ad_xyz] = aseinds
        
    ts_dict = dict()
    ts_admol_dict = dict()
    ts_neighbor_sites_2D_dict = dict()
    ts_ninds_dict = dict()
    ts_neighbor_sites_dict = dict()
    ts_actual_occ_dict = dict()
    ts_aseinds_dict = dict()
    ts_xyz = None
    for ts_xyz in ts_xyzs:
        ad,admol,neighbor_sites_2D,ninds,actual_occ,neighbor_sites,aseinds,slab,nslab = get_adsorbate_ts_information(ts_xyz,slabxyz,True,metal,facet,sites,site_adjacency,max_dist=max_dist)
        ts_dict[ts_xyz] = ad
        ts_admol_dict[ts_xyz] = admol
        ts_neighbor_sites_2D_dict[ts_xyz] = neighbor_sites_2D
        ts_ninds_dict[ts_xyz] = ninds
        ts_neighbor_sites_dict[ts_xyz] = neighbor_sites
        ts_actual_occ_dict[ts_xyz] = actual_occ
        ts_aseinds_dict[ts_xyz] = aseinds

    
    outatoms = []
    outmol2Dsts = []
    outmol2Dsad = []
    outxyzsts = []
    outxyzsad = []
    
    for coadname in coadnames:
        for j,ts_xyz in enumerate(ts_xyzs):
            ts_name = "TS" + str(j)
            print(ts_name)
            if not os.path.exists(os.path.join(outdir,ts_name)):
                os.makedirs(os.path.join(outdir,ts_name))
            shutil.copyfile(os.path.join(os.path.split(os.path.split(ts_xyz)[0])[0],"info.json"),
                            os.path.join(outdir,ts_name,"info.json"))
            coad_to_stable_neighbor_sites,coad_site_stable_parameters,coad_atom_to_molecule_surface_atom_map,infocoad_dict,stable_neighbor_sites_total,coads_dict,coad2Ds = get_coadsorbate_information(coadnames,
                                                                        ads_dir,ts_neighbor_sites_dict[ts_xyz],sites,site_adjacency,nslab,slab)
            atoms,mol2Ds = generate_coadsorbed_geoms(ts_dict[ts_xyz],
                                                          ts_admol_dict[ts_xyz],
                                                          ts_neighbor_sites_2D_dict[ts_xyz],ts_ninds_dict[ts_xyz],
                                                          ts_actual_occ_dict[ts_xyz],ts_neighbor_sites_dict[ts_xyz],
                               ts_aseinds_dict[ts_xyz],slab,nslab,True,ads_dir,unstable_pairs,coadname,
                               metal,facet,sites,site_adjacency,coad_to_stable_neighbor_sites, coad_site_stable_parameters,
                               coad_atom_to_molecule_surface_atom_map,infocoad_dict,
                               stable_neighbor_sites_total,coads_dict,coad2Ds,max_dist=max_dist)
            for i,mol2D in enumerate(mol2Ds): #check if we already have it
                for mol2Dout in outmol2Dsts:
                    if mol2Dout.is_isomorphic(mol2D,save_order=True):
                        break
                else:
                    outatoms.append(atoms[i])
                    outmol2Dsts.append(mol2Ds[i])
                    if not os.path.exists(os.path.join(outdir,ts_name,coadname,str(i))):
                        os.makedirs(os.path.join(outdir,ts_name,coadname,str(i)))
                    with open(os.path.join(outdir,ts_name,coadname,str(i),"info.json"),"w") as f:
                        d = {"adjlist": mol2D.to_adjacency_list(),
                            "xyz": ts_xyz}
                        json.dump(d,f)
                    write(os.path.join(outdir,ts_name,coadname,str(i),"init.xyz"),atoms[i])
                    outxyzsts.append(os.path.join(outdir,ts_name,coadname,str(i),"init.xyz"))
        
        for j,ad_xyz in enumerate(ad_xyzs):
            with open(os.path.join(os.path.split(os.path.split(ad_xyz)[0])[0],"info.json"),"r") as f:
                info = json.load(f)
            ad_name = info["name"]
            print(ad_name)
            if not os.path.exists(os.path.join(outdir,ad_name)):
                os.makedirs(os.path.join(outdir,ad_name))
            shutil.copyfile(os.path.join(os.path.split(os.path.split(ad_xyz)[0])[0],"info.json"),
                                         os.path.join(outdir,ad_name,"info.json"))
            coad_to_stable_neighbor_sites,coad_site_stable_parameters,coad_atom_to_molecule_surface_atom_map,infocoad_dict,stable_neighbor_sites_total,coads_dict,coad2Ds = get_coadsorbate_information(coadnames,
                                                                        ads_dir,ad_neighbor_sites_dict[ad_xyz],sites,site_adjacency,nslab,slab)
            atoms,mol2Ds = generate_coadsorbed_geoms(ad_dict[ad_xyz],
                                                          ad_admol_dict[ad_xyz],
                                                          ad_neighbor_sites_2D_dict[ad_xyz],ad_ninds_dict[ad_xyz],
                                                          ad_actual_occ_dict[ad_xyz],ad_neighbor_sites_dict[ad_xyz],
                               ad_aseinds_dict[ad_xyz],slab,nslab,False,ads_dir,unstable_pairs,coadname,
                               metal,facet,sites,site_adjacency,coad_to_stable_neighbor_sites, coad_site_stable_parameters,
                               coad_atom_to_molecule_surface_atom_map,infocoad_dict,
                               stable_neighbor_sites_total,coads_dict,coad2Ds,max_dist=max_dist)
            for i,mol2D in enumerate(mol2Ds): #check if we already have it
                for mol2Dout in outmol2Dsad:
                    if mol2Dout.is_isomorphic(mol2D):
                        break
                else:
                    outatoms.append(atoms[i])
                    outmol2Dsad.append(mol2Ds[i])
                    if not os.path.exists(os.path.join(outdir,ad_name,coadname,str(i))):
                        os.makedirs(os.path.join(outdir,ad_name,coadname,str(i)))
                    with open(os.path.join(outdir,ad_name,coadname,str(i),"info.json"),"w") as f:
                        d = {"adjlist": mol2D.to_adjacency_list(),
                            "xyz": ad_xyz}
                        json.dump(d,f)
                    write(os.path.join(outdir,ad_name,coadname,str(i),"init.xyz"),atoms[i])
                    outxyzsad.append(os.path.join(outdir,ad_name,coadname,str(i),"init.xyz"))
        
    return outatoms,outmol2Dsad,outmol2Dsts,outxyzsad,outxyzsts

def generate_coadsorbed_geoms(ad,admol,neighbor_sites_2D,ninds,actual_occ,neighbor_sites,
                               aseinds,slab,nslab,is_ts,ads_dir,unstable_pairs,coadname,
                               metal,facet,sites,site_adjacency,coad_to_stable_neighbor_sites, coad_site_stable_parameters,
                               coad_atom_to_molecule_surface_atom_map,infocoad_dict,
                               stable_neighbor_sites_total,coads_dict,coad2Ds,max_dist=3.0):
    
    logging.error("stable neighbor sites: {}".format(len(stable_neighbor_sites_total)))
    
    coads = coads_dict[coadname]
    atom_to_molecule_surface_atom_map = coad_atom_to_molecule_surface_atom_map[coadname]
    infocoad = infocoad_dict[coadname]
    coad2D = coad2Ds[coadname]
    site_stable_parameters = coad_site_stable_parameters[coadname]
    coad_occ_dict = {coads.index(coad): get_occupied_sites(coad,sites,nslab) for coad in coads}
    coad_height_map = {coads.index(coad): list(get_bond_lengths_sites(Molecule().from_adjacency_list(infocoad["adjlist"]),
                                                                       coad,
                                         {int(x):y for x,y in infocoad["atom_to_molecule_atom_map"].items()},
                                        {int(x):y for x,y in infocoad["gratom_to_molecule_surface_atom_map"].items()},
                                                                       infocoad["nslab"],sites,site_adjacency,
                                                                facet=facet,metal=metal)[2].values())[0] for coad in coads}
    outgeoms = [ad]
    outmol2Ds = [admol]
    geo_fails = 0
    mol2D_fails = 0
    config_fails = 0
    site_fails = 0
    unique_fails = 0
    for i,site in enumerate(stable_neighbor_sites_total):
        logging.error("doing site {}".format(i))
        newoutgeoms = []
        newoutmol2Ds = []
        site_2D_inds = [i for i,x in enumerate(neighbor_sites_2D) if sites_match(site,x,slab)]
        if not site_2D_inds:
            site_fails += 1
            continue
        
        for j,geo in enumerate(outgeoms):
            mol2D = outmol2Ds[j]
            geo = deepcopy(geo)
            geo,coad2D = add_coadsorbate_3D(geo,site,ad,coads,site_stable_parameters,
                        atom_to_molecule_surface_atom_map,infocoad,coad_occ_dict,coad_height_map,coad2D,
                       metal,facet,sites,site_adjacency,nslab)
            if geo is None:
                geo_fails += 1
                continue
            mol2D = mol2D.copy(deep=True)
            try:
                mol2D = add_coadsorbate_2D(mol2D,site,coad2D,slab,neighbor_sites_2D,site_2D_inds)
            except Exception as e:
                print(mol2D.to_adjacency_list())
                print(site)
                print(site_2D_inds)
                raise e
            
            if mol2D is None:
                mol2D_fails += 1
                continue
            
            if configuration_is_valid(mol2D,admol,is_ts,unstable_pairs):
                for m in outmol2Ds:
                    if mol2D.is_isomorphic(m,save_order=True):
                        unique_fails += 1
                        break
                else:
                    assert len(geo) - nslab == len(mol2D.atoms) - len([a for a in mol2D.atoms if a.is_surface_site()])
                    newoutgeoms.append(geo)
                    newoutmol2Ds.append(mol2D)

        outgeoms.extend(newoutgeoms)
        outmol2Ds.extend(newoutmol2Ds)
        logging.error("added so far: {}".format(len(outgeoms)))
    
    
    outgeoms.remove(ad) #do not include the configuration with no coadsorbates in output
    outmol2Ds.remove(admol)
    return outgeoms,outmol2Ds



def configuration_is_valid(mol2D,admol,is_ts,unstable_pairs):
    unstable_ind_pairs = set()
    if is_ts:
        admol_splits = split_ts_to_reactants(admol,tagatoms=False)
        for asplit in admol_splits:
            snum = len([a for a in asplit.atoms if a.is_surface_site()])
            for unstable_pair in unstable_pairs:
                iso = asplit.find_subgraph_isomorphisms(unstable_pair,save_order=True)
                if iso:
                    inds = []
                    for a in iso[0].keys():
                        assert a in asplit.atoms
                        if a.is_bonded_to_surface() and not a.is_surface_site():
                            inds.append(asplit.atoms.index(a)-snum)
                    unstable_ind_pairs.add(frozenset(inds)) #groups of atom inds that if they are in a isomorphism indicate to allow anything for that ts split

    struct = mol2D
    if is_ts:
        structspl = split_ts_to_reactants(struct,tagatoms=False)
    else:
        structspl = [struct]
    
    validity_judgements = []
    failed = False
    for st in structspl:
        snum = len([a for a in st.atoms if a.is_surface_site()])
        failed = False
        for unstable_pair in unstable_pairs:
            iso = st.find_subgraph_isomorphisms(unstable_pair,save_order=True)

            if iso:
                inds = []
                for a in iso[0].keys():
                    if a.is_bonded_to_surface() and not a.is_surface_site():
                        inds.append(st.atoms.index(a)-snum)
                if frozenset(inds) in unstable_ind_pairs:
                    pass
                else:
                    failed = True
        if not failed:
            validity_judgements.append(True)
        else:
            validity_judgements.append(False)

    return all(validity_judgements)


def adsorbate_interaction_decomposition(mol,local_radius=5):
    surface_bonded_inds = []
    for i,at in enumerate(mol.atoms):
        if at.is_bonded_to_surface() and not at.is_surface_site():
            surface_bonded_inds.append(i)
    
    structs = []
    for i,indi in enumerate(surface_bonded_inds):
        for j,indj in enumerate(surface_bonded_inds):
            if i > j:
                st = mol.copy(deep=True)
                st.atoms[indi].label = "*"
                st.atoms[indj].label = "*"
                paths = find_shortest_paths_sites(st.atoms[indi],st.atoms[indj])
                ad_atoms_i,ad_surf_sites_i = find_adsorbate_atoms_surface_sites(st.atoms[indi],st)
                ad_atoms_j,ad_surf_sites_j = find_adsorbate_atoms_surface_sites(st.atoms[indj],st)
                if st.atoms[indi] in ad_atoms_j: #i & j are part of same adsorbate/TS structure
                    continue
                for p in paths:
                    for a in p:
                        if a.is_surface_site():
                            a.label = "*"
                            
                for ind in [indi,indj]:
                    for a in st.atoms:
                        if a.is_surface_site() and len(find_shortest_paths_sites(st.atoms[indi],a)[0]) - 2 <= local_radius:
                            a.label = "*"
                
                ats_to_delete = [a for a in st.atoms if a.is_surface_site() and a.label != "*" and a not in ad_surf_sites_i and a not in ad_surf_sites_j]
                for a in ats_to_delete:
                    st.remove_atom(a)
                try:
                    stout = [x for x in st.split() if x.contains_surface_site()][0]
                except IndexError as e:
                    print(mol.to_adjacency_list())
                    print(st.to_adjacency_list())
                    raise e
                for a in stout.atoms:
                    if a.is_surface_site():
                        a.label = ""
                
                structs.append(stout)
    
    return structs
    
def adsorbate_site_decomposition(mol):
    surface_bonded_inds = []
    for i,at in enumerate(mol.atoms):
        if at.is_bonded_to_surface() and not at.is_surface_site():
            surface_bonded_inds.append(i)
    
    structs = []
    for i,indi in enumerate(surface_bonded_inds):
        st = mol.copy(deep=True)
        st.atoms[indi].label = "*"
        structs.append(st)
    
    return structs

def get_adsorbed_atom_pairs(length=7, r_bonds=None):
    """
    length is number of site atoms between adsorbates
    """
    if r_bonds is None:
        r_bonds = [1, 2, 3, 0.05]
    groups = []
    for j in range(2,length+1):
        g = Group().from_adjacency_list("""1 * R u0 px cx""")
        a2 = g.atoms[0]
        for i in range(j):
            a = GroupAtom(atomtype=["X"],radical_electrons=[0],lone_pairs=[0],charge=[0])
            if i == 0:
                b = GroupBond(a2, a, order=r_bonds)
            else:
                b = GroupBond(a2, a, order=["S"])
            g.add_atom(a)
            g.add_bond(b)
            a2 = a
        a = GroupAtom(atomtype=["R"], label="*", radical_electrons=[0])
        b = GroupBond(a2, a, order=r_bonds)
        g.add_atom(a)
        g.add_bond(b)
        groups.append(g)

    return groups

import itertools
def get_adsorbed_atom_groups(Nad=3, length=7, r_bonds=None):
    """
    length is number of site atoms between adsorbates
    """
    assert Nad == 3, "Doesn't work for Nad=2 and the combination ordering seems to matter for Nad > 3, so only Nad=3"
    if r_bonds is None:
        r_bonds = [1, 2, 3, 0.05]
    groups = []
    lengths = list(range(2,length+1))
    for comb in itertools.combinations(lengths,Nad):
        g = Group().from_adjacency_list("""1 * R u0 px cx""")
        a2 = g.atoms[0]
        for k,j in enumerate(comb):
            for i in range(j):
                a = GroupAtom(atomtype=["X"],radical_electrons=[0],lone_pairs=[0],charge=[0])
                if i == 0:
                    b = GroupBond(a2, a, order=r_bonds)
                else:
                    b = GroupBond(a2, a, order=["S"])
                g.add_atom(a)
                g.add_bond(b)
                a2 = a
            if k+1 < len(comb):
                a = GroupAtom(atomtype=["R"], label="*", radical_electrons=[0])
                b = GroupBond(a2, a, order=r_bonds)
                g.add_atom(a)
                g.add_bond(b)
            else:
                b = GroupBond(g.atoms[1],a2, order=["S"])
                g.add_bond(b)
        
        groups.append(g)

    return groups

def get_atom_centered_correction(m,coadmol_E_dict):
    out_structs = split_adsorbed_structures(m,clear_site_info=False)
    correction = 0.0
    minE = min(coadmol_E_dict.values())
    for struct in out_structs:
        for coadmol,E in coadmol_E_dict.items():
            if struct.is_isomorphic(coadmol,save_order=True):
                correction += E - minE
                break
    return correction

def get_atom_center_stability(m,coadmol_stability_dict):
    out_structs = split_adsorbed_structures(m,clear_site_info=False)
    for struct in out_structs:
        for coadmol,v in coadmol_stability_dict.items():
            if struct.is_isomorphic(coadmol,save_order=True):
                if not v:
                    return False

    return True

def break_train_val_test(datums,test_fract=0.1,val_fract=0.1):
    ds = datums[:]
    N = len(datums)
    Ntest = round(N*test_fract)
    Nval = round(N*val_fract)
    np.random.shuffle(ds)
    test = ds[:Ntest]
    val = ds[Ntest:Ntest+Nval]
    train = ds[Nval+Ntest:]
    return train,val,test

def split_triad(tagged_grp):
    labeled_atoms = tagged_grp.get_labeled_atoms("*")
    
    pair_atoms_inds = []
    a = labeled_atoms[0]
    path = [a]
    notatom = None
    c = 0
    while c < 3:
        if len(path) == 1:
            anew = list(path[-1].bonds.keys())[0]
        elif len(path[-1].bonds) == 3 and len(path) > 2:
            anew = [x for x in path[-1].bonds.keys() if not x.is_surface_site()][0]
        elif len(path[-1].bonds) == 3:
            ats = list(path[-1].bonds.keys())
            for a in ats:
                if a not in path and a is not notatom:
                    anew = a
                    break
            else:
                raise ValueError
        elif len(path[-1].bonds) == 2:
            ats = list(path[-1].bonds.keys())
            if ats[0] in path:
                anew = ats[1]
            else:
                anew = ats[0]
        path.append(anew)
        
        if not anew.is_surface_site():
            pair_atoms_inds.append([tagged_grp.atoms.index(a) for a in path])
            c += 1
            notatom = path[-3]
            path = path[-2:][::-1]
            
    out_pairs = []
    for pair_inds in pair_atoms_inds:
        g = tagged_grp.copy(deep=True)
        ats = [g.atoms[ind] for ind in pair_inds]
        atom_to_remove = []
        for a in g.atoms:
            if a not in ats:
                atom_to_remove.append(a)
        for a in atom_to_remove:
            g.remove_atom(a)
        g.update()
        out_pairs.append(g)
        
    return out_pairs

def is_descendent_of_or_is(node,ancestor_node):
    n = node
    while n.parent is not None:
        if n is ancestor_node:
            return True
        else:
            n = n.parent
    return False

#currently not used anymore as the default algorithm is sufficient for pair-wise decompositions
class CoverageDependenceRegressor(MultiEvalSubgraphIsomorphicDecisionTreeRegressor):
    def fit_rule(self, alpha=0.1):
        max_depth = max([node.depth for node in self.nodes.values()])
        y = np.array([datum.value for datum in self.datums])
        preds = np.zeros(len(self.datums))
        self.node_uncertainties = dict()
        weights = self.weights
        W = self.W
        triad_node = self.nodes["Root_Triad"]
        for pair in [True,False]:
            for depth in range(max_depth + 1):
                if depth == 0:
                    self.nodes["Root"].rule = Rule(value=0.0,num_data=0)
                    continue #skip Root node
                else:
                    if pair:
                        nodes = [node for node in self.nodes.values() if node.depth == depth and not is_descendent_of_or_is(node,triad_node)]
                    else:
                        nodes = [node for node in self.nodes.values() if node.depth == depth and is_descendent_of_or_is(node,triad_node)]
                
                if len(nodes) == 0:
                    continue
                    
                # generate matrix
                A = sp.lil_matrix((len(self.datums), len(nodes)))
                y -= preds
    
                for i, datum in enumerate(self.datums):
                    for node in self.mol_node_maps[datum]["nodes"]:
                        while node is not None:
                            if node in nodes:
                                j = nodes.index(node)
                                A[i, j] += 1.0
                            node = node.parent
    
                clf = linear_model.Lasso(
                    alpha=alpha,
                    fit_intercept=False,
                    tol=1e-4,
                    max_iter=1000000000,
                    selection="random",
                )
                if weights is not None:
                    lasso = clf.fit(A, y, sample_weight=weights)
                else:
                    lasso = clf.fit(A, y)
                
                preds = A * clf.coef_
                self.data_delta = preds - y
    
                for i, val in enumerate(clf.coef_):
                    nodes[i].rule = Rule(value=val, num_data=np.sum(A[:, i]))

        train_error = [self.evaluate(d.mol, estimate_uncertainty=False) - d.value for d in self.datums]

        logging.info("training MAE: {}".format(np.mean(np.abs(np.array(train_error)))))

        if self.validation_set:
            val_error = [self.evaluate(d.mol, estimate_uncertainty=False) - d.value for d in self.validation_set]
            val_mae = np.mean(np.abs(np.array(val_error)))
            if val_mae < self.min_val_error:
                self.min_val_error = val_mae
                self.best_tree_nodes = list(self.nodes.keys())
                self.bestA = A
                self.best_nodes = {k: v for k, v in self.nodes.items()}
                self.best_mol_node_maps = {
                    k: {"mols": v["mols"][:], "nodes": v["nodes"][:]}
                    for k, v in self.mol_node_maps.items()
                }
                self.best_rule_map = {name:self.nodes[name].rule for name in self.best_tree_nodes}
            self.val_mae = val_mae
            logging.info("validation MAE: {}".format(self.val_mae))

        if self.test_set:
            test_error = [self.evaluate(d.mol) - d.value for d in self.test_set]
            test_mae = np.mean(np.abs(np.array(test_error)))
            logging.info("test MAE: {}".format(test_mae))
            
        logging.info("# nodes: {}".format(len(self.nodes)))
        
def train_sidt_cov_dep_regressor(pairs_datums,sampling_datums,r_site=None,
                                r_atoms=None,node_fract_training=0.7):

    if r_site is None:
        r_site = ["","ontop","bridge","hcp","fcc"]

    if r_atoms is None:
        r_atoms = ["C","O","N","H","X"]
    
    root_pair = Group().from_adjacency_list("""
    1 * R u0 px cx {2,[S,D,T,Q,R,vdW]}
    2   X u0 p0 c0 {1,[S,D,T,Q,R,vdW]}
    3 * R u0 px cx {4,[S,D,T,Q,R,vdW]}
    4   X u0 p0 c0 {3,[S,D,T,Q,R,vdW]}
    """)

    root_pair_node = Node(group=root_pair,name="Root",parent=None,depth=0)

    pair_nodes = []
    initial_pair_root_splits = get_adsorbed_atom_pairs(length=7,r_bonds=[1,2,3,0.05,0])
    for i,g in enumerate(initial_pair_root_splits):
        gcleared = g.copy(deep=True)
        gcleared.clear_labeled_atoms()
        for d in pairs_datums+sampling_datums:
            if d.mol.is_subgraph_isomorphic(gcleared,save_order=True): #only add if the group exists
                n = Node(group=g,name=root_pair_node.name+"_"+str(i),parent=root_pair_node,depth=1)
                pair_nodes.append(n)
                root_pair_node.children.append(n)
                break

    pairnodes = {n.name:n for n in [root_pair_node]+pair_nodes}

    Nfullnodes = (len(pairs_datums)+len(sampling_datums))*node_fract_training

    try: #pysidt >1.0
        treepair = MultiEvalSubgraphIsomorphicDecisionTreeRegressor([adsorbate_interaction_decomposition],
                                                   nodes=pairnodes,
                                                   r=[ATOMTYPES[x] for x in r_atoms],
                                                   r_bonds=[1,2,3,0.05],
                                                   r_un=[0],
                                                   r_site=r_site,
                                                   max_structures_to_generate_extensions=100,
                                                   fract_nodes_expand_per_iter=0.025,
                                                   iter_max=2,
                                                   iter_item_cap=100,
                                                   weigh_node_selection_by_occurrence=False,
                                                  )
    except: #pysidt 1.0.0
        treepair = MultiEvalSubgraphIsomorphicDecisionTreeRegressor([adsorbate_interaction_decomposition],
                                                   nodes=pairnodes,
                                                   r=[ATOMTYPES[x] for x in r_atoms],
                                                   r_bonds=[1,2,3,0.05],
                                                   r_un=[0],
                                                   r_site=r_site,
                                                   max_structures_to_generate_extensions=100,
                                                   fract_nodes_expand_per_iter=0.025,
                                                   iter_max=2,
                                                   iter_item_cap=100,
                                                  )
        
    treepair.generate_tree(data=pairs_datums+sampling_datums,max_nodes=Nfullnodes)
    
    return treepair

def train_sidt_cov_dep_classifier(datums,r_site=None,
                                r_atoms=None,node_fract_training=0.5,node_min=50):

    if r_site is None:
        r_site = ["","ontop","bridge","hcp","fcc"]

    if r_atoms is None:
        r_atoms = ["C","O","N","H","X"]
    
    root_node = Node(group=None,name="Root",depth=0)
    
    root_pair = Group().from_adjacency_list("""
    1 * R u0 px cx {2,[S,D,T,Q,R]}
    2   X u0 p0 c0 {1,[S,D,T,Q,R]}
    3 * R u0 px cx {4,[S,D,T,Q,R]}
    4   X u0 p0 c0 {3,[S,D,T,Q,R]}
    """)
    
    root_triad = Group().from_adjacency_list("""
    1 * R u0 px cx {2,[S,D,T,Q,R]}
    2   X u0 p0 c0 {1,[S,D,T,Q,R]}
    3 * R u0 px cx {4,[S,D,T,Q,R]}
    4   X u0 p0 c0 {3,[S,D,T,Q,R]}
    5 * R u0 px cx {6,[S,D,T,Q,R]}
    6   X u0 p0 c0 {5,[S,D,T,Q,R]}
    """)

    root_pair_node = Node(group=root_pair,name="Root_Pair",parent=root_node,depth=1)
    root_triad_node = Node(group=root_triad,name="Root_Triad",parent=root_node,depth=1)
    root_node.children.extend([root_pair_node,root_triad_node])

    pair_nodes = []
    initial_root_splits = get_adsorbed_atom_pairs(length=7,r_bonds=[1,2,3,0.05])
    for i,g in enumerate(initial_root_splits):
        n = Node(group=g,name=root_pair_node.name+"_"+str(i),parent=root_pair_node,depth=2)
        pair_nodes.append(n)
        root_pair_node.children.append(n)

    triad_nodes = []
    initial_root_splits = get_adsorbed_atom_groups(Nad=3,length=7,r_bonds=[1,2,3,0.05])
    for i,g in enumerate(initial_root_splits):
        n = Node(group=g,name=root_triad_node.name+"_"+str(i),parent=root_triad_node,depth=2)
        triad_nodes.append(n)
        root_triad_node.children.append(n)

    nodes = {n.name:n for n in [root_node,root_pair_node,root_triad_node]+pair_nodes+triad_nodes}
    
    node_len = len(nodes)

    tree = MultiEvalSubgraphIsomorphicDecisionTreeBinaryClassifier([adsorbate_interaction_decomposition,adsorbate_triad_interaction_decomposition],
                                                                   nodes=nodes,
                                               r=[ATOMTYPES[x] for x in r_atoms],
                                               r_bonds=[1,2,3,0.05],
                                                         r_un=[0],
                                               r_site=r_site,
                                                iter_max=2,
                                                iter_item_cap=100,
                                                max_structures_to_generate_extensions=100,
                                              )

    train_sample,val,test = break_train_val_test(datums,test_fract=0.0,val_fract=0.1)

    if len(train_sample)*node_fract_training > node_min:
        tree.generate_tree(data=train_sample,validation_set=val,max_nodes=len(train_sample)*node_fract_training,
                       postpruning_based_on_val=True)
    else:
        tree.generate_tree(data=train_sample,validation_set=val,max_nodes=node_min)

    return tree

def add_ad_to_site(admol,coad,site):
    c = coad.copy(deep=True)
    at = [a for a in c.atoms if a.is_bonded_to_surface()][0]
    bd = [bd for a,bd in at.bonds.items() if a.is_surface_site()][0]
    order = bd.order
    for a in c.atoms:
        if a.is_surface_site():
            c.remove_atom(a)
    
    atind = c.atoms.index(at)
    sind = admol.atoms.index(site)
    admol_length = len(admol.atoms)
    
    admolout = admol.copy(deep=True).merge(c)
    snew = admolout.atoms[sind]
    atnew = admolout.atoms[len(admol.atoms)+atind]
    assert site.site == snew.site and site.morphology == snew.morphology
    assert atnew.element == at.element
    
    bd = Bond(atnew,snew,order=order)
    admolout.add_bond(bd)
    
    admolout.clear_labeled_atoms()
    admolout.multiplicity=1
    try:
        admolout.update(sort_atoms=False)
    except Exception as e:
        raise ValueError((e,admol.to_adjacency_list(),coad.to_adjacency_list(),admolout.to_adjacency_list()))
    admolout.update_connectivity_values()
    
    return admolout
    
def get_configurations(admol, coad, coad_stable_sites, tree_interaction_classifier=None, coadmol_stability_dict=None, unstable_groups=None,
                       tree_interaction_regressor=None, tree_atom_regressor=None, coadmol_E_dict=None, energy_tol=None):
    empty_sites = [a for a in admol.atoms if a.is_surface_site() and a.site in coad_stable_sites and not any([not a2.is_surface_site() for a2 in a.bonds.keys()])]
    print("empty sites")
    print(len(empty_sites))
    empty_site_inds = [admol.atoms.index(s) for s in empty_sites]
    outmols = [admol]
    outmols_split = [[admol]]
    lowest_energy_surface_bond_dict = dict()
    print(len(outmols))
    for i in range(len(empty_sites)):
        newoutmols = []
        for m in outmols_split[i]: #all configurations with a fixed number of co-adsorbates (highest number hit yet)
            for sind in empty_site_inds:
                if any(not a.is_surface_site() for a in m.atoms[sind].bonds.keys()): #site is already occupied 
                    continue
                
                outmol = add_ad_to_site(m,coad,m.atoms[sind])
    
                if (tree_interaction_classifier and not tree_interaction_classifier.evaluate(outmol)) or \
                    (coadmol_stability_dict and not get_atom_center_stability(outmol,coadmol_stability_dict)): #unstable configuration
                    continue

                if unstable_groups:
                    is_unstable = False
                    for grp in unstable_groups:
                        if outmol.is_subgraph_isomorphic(grp,save_order=True):
                            is_unstable = True 
                            break
                    if is_unstable:
                        continue
                    
                for mol in newoutmols:
                    if outmol.is_isomorphic(mol,save_order=True):
                        break
                else:
                    if energy_tol: #skip configurations that are more than energy_tol higher than the lowest energy configuraiton for a given coverage found
                        if tree_atom_regressor is not None:
                            E = tree_atom_regressor.evaluate(m) + tree_interaction_regressor.evaluate(m)
                        elif coadmol_E_dict is not None:
                            E = get_atom_centered_correction(m,coadmol_E_dict) + tree_interaction_regressor.evaluate(m)
                        else:
                            raise ValueError
    
                        Nsurfbonds = len([b for b in m.get_all_edges() if (b.atom1.is_surface_site() and not b.atom2.is_surface_site()) or (b.atom2.is_surface_site() and not b.atom1.is_surface_site())])
                        if Nsurfbonds not in lowest_energy_surface_bond_dict.keys() or lowest_energy_surface_bond_dict[Nsurfbonds] > E:
                            lowest_energy_surface_bond_dict[Nsurfbonds] = E
                            newoutmols.append(outmol)
                        elif lowest_energy_surface_bond_dict[Nsurfbonds] + energy_tol < E:
                            pass #don't add configuration
                        else:
                            newoutmols.append(outmol)
                    else:
                        newoutmols.append(outmol)
                        
        outmols.extend(newoutmols)
        outmols_split.append(newoutmols)
        print(len(outmols))

    return outmols

def check_stable(config,stability_datums):
    for d in stability_datums:
        if d.value == False and config.is_isomorphic(d.mol,save_order=True):
            return False
    else:
        return True

def evaluate_from_datums(config,energy_datums):
    for d in energy_datums:
        if config.is_isomorphic(d.mol,save_order=True):
            return d.value
    else:
        return None

def get_cov_energies_configs_concern_tree(tree_interaction_regressor, configs, coad_stable_sites, Ncoad_isolated, concern_energy_tol=None, tree_atom_regressor=None, coadmol_E_dict=None, 
                     stability_datums=None):
    Ncoad_energy_dict = dict()
    Ncoad_config_dict = dict()
    config_to_Eunctr = dict()
    configs_of_concern = {}
    Nempty = len([a for a in configs[0].atoms if a.is_surface_site() and a.site in coad_stable_sites])
    for i,m in enumerate(configs):
        Ncoad = len([a for a in m.atoms if a.is_surface_site() and any(not a2.is_surface_site() for a2 in a.bonds.keys())]) - Ncoad_isolated

        if stability_datums:
            stable = check_stable(m,stability_datums)
        else:
            stable = True
            
        if not stable:
            continue
            
        if tree_atom_regressor is not None:
            Einteraction,std,tr = tree_interaction_regressor.evaluate(m,trace=True, estimate_uncertainty=True)
            E = tree_atom_regressor.evaluate(m) + Einteraction
        elif coadmol_E_dict is not None:
            Einteraction,std,tr = tree_interaction_regressor.evaluate(m,trace=True, estimate_uncertainty=True)
            E = get_atom_centered_correction(m,coadmol_E_dict) + Einteraction
        else:
            raise ValueError
        
        config_to_Eunctr[i] = (m,E,std,tr)
        
        if Ncoad not in Ncoad_energy_dict.keys():
            Ncoad_energy_dict[Ncoad] = E
            Ncoad_config_dict[Ncoad] = m.to_adjacency_list()
        elif E < Ncoad_energy_dict[Ncoad]:
            Ncoad_energy_dict[Ncoad] = E
            Ncoad_config_dict[Ncoad] = m.to_adjacency_list()
        elif E == Ncoad_energy_dict[Ncoad]:
            if isinstance(Ncoad_config_dict[Ncoad],list):
                Ncoad_config_dict[Ncoad].append(m.to_adjacency_list())
            else:
                Ncoad_config_dict[Ncoad] = [Ncoad_config_dict[Ncoad], m.to_adjacency_list()]
    
    for i in config_to_Eunctr.keys():
        Ncoad = len([a for a in m.atoms if a.is_surface_site() and any(not a2.is_surface_site() for a2 in a.bonds.keys())]) - Ncoad_isolated
        m,E,std,tr = config_to_Eunctr[i]
        if ((concern_energy_tol is None) or (Ncoad_energy_dict[Ncoad] + concern_energy_tol > E)):
            configs_of_concern[i] = (m,E,tr,std)
    
    return Ncoad_energy_dict,Ncoad_config_dict,configs_of_concern
    
def get_cov_energies(tree_interaction_regressor, configs, coad_stable_sites, Ncoad_isolated, tree_atom_regressor=None, coadmol_E_dict=None, 
                     stability_datums=None):
    Ncoad_energy_dict = dict()
    Ncoad_config_dict = dict()
    Nempty = len([a for a in configs[0].atoms if a.is_surface_site() and a.site in coad_stable_sites])
    for m in configs:
        Ncoad = len([a for a in m.atoms if a.is_surface_site() and any(not a2.is_surface_site() for a2 in a.bonds.keys())]) - Ncoad_isolated

        if tree_atom_regressor is not None:
            E = tree_atom_regressor.evaluate(m) + tree_interaction_regressor.evaluate(m)
        elif coadmol_E_dict is not None:
            E = get_atom_centered_correction(m,coadmol_E_dict) + tree_interaction_regressor.evaluate(m)
        else:
            raise ValueError
        if Ncoad not in Ncoad_energy_dict.keys():
            if (not stability_datums) or check_stable(m,stability_datums):
                Ncoad_energy_dict[Ncoad] = E
                Ncoad_config_dict[Ncoad] = m.to_adjacency_list()
        elif E < Ncoad_energy_dict[Ncoad]:
            if (not stability_datums) or check_stable(m,stability_datums):
                Ncoad_energy_dict[Ncoad] = E
                Ncoad_config_dict[Ncoad] = m.to_adjacency_list()
        elif E == Ncoad_energy_dict[Ncoad]:
            if (not stability_datums) or check_stable(m,stability_datums):
                if isinstance(Ncoad_config_dict[Ncoad],list):
                    Ncoad_config_dict[Ncoad].append(m.to_adjacency_list())
                else:
                    Ncoad_config_dict[Ncoad] = [Ncoad_config_dict[Ncoad], m.to_adjacency_list()]
    return Ncoad_energy_dict,Ncoad_config_dict

def get_configs_of_concern(tree_interaction_regressor,configs,coad_stable_sites,Ncoad_energy_dict,Nocc_isolated,concern_energy_tol,tree_atom_regressor=None,
                           coadmol_E_dict=None):
    configs_of_concern = {}
    Nempty = len([a for a in configs[0].atoms if a.is_surface_site() and a.site in coad_stable_sites])
    for m in configs:
        Ncoad = len([a for a in m.atoms if a.is_surface_site() and any(not a2.is_surface_site() for a2 in a.bonds.keys())]) - Nocc_isolated
        if tree_atom_regressor is not None:
            Einteraction,std,tr = tree_interaction_regressor.evaluate(m,trace=True, estimate_uncertainty=True)
            E = tree_atom_regressor.evaluate(m) + Einteraction
        elif coadmol_E_dict is not None:
            Einteraction,std,tr = tree_interaction_regressor.evaluate(m,trace=True, estimate_uncertainty=True)
            E = get_atom_centered_correction(m,coadmol_E_dict) + Einteraction
        else:
            raise ValueError
        if Ncoad_energy_dict[Ncoad] + concern_energy_tol > E:
            configs_of_concern[m] = (E,tr,std)
        
    return configs_of_concern

def load_coverage_delta(d,ad_energy_dict,slab,metal,facet,sites,site_adjacency,ts_pynta_dir=None,allowed_structure_site_structures=None,
                       out_file_name="out",vib_file_name="vib",is_ad=None,keep_binding_vdW_bonds=False,keep_vdW_surface_bonds=False):

    try:
        info = json.load(open(os.path.join(d,"info.json")))
    except:
        info = None

    if is_ad is None:
        is_ts =  "xyz" in info.keys() and os.path.split(os.path.split(os.path.split(info["xyz"])[0])[0])[1][:2] == "TS"
    else:
        is_ts = not is_ad
    
    atoms = read(os.path.join(d,out_file_name+".xyz"))
    
    if not is_ts:
        try:
            admol,neighbor_sites,ninds = generate_adsorbate_2D(atoms, sites, site_adjacency, len(slab), max_dist=np.inf, cut_off_num=None, 
                                                               allowed_structure_site_structures=allowed_structure_site_structures,
                                                               keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
        except (SiteOccupationException,TooManyElectronsException,ValueError,FailedFixBondsException) as e:
            return None,None,None,None

        logging.error(admol.to_adjacency_list())
        split_structs = split_adsorbed_structures(admol)
        try:
            vibdata = get_vibdata(os.path.join(d,out_file_name+".xyz"),os.path.join(d,vib_file_name+".json"),len(slab))

            Ecad = atoms.get_potential_energy() - slab.get_potential_energy() + vibdata.get_zero_point_energy()
    
            Esep = 0.0
            for split_struct in split_structs:
                # if not split_struct.contains_surface_site(): #if gas phase species don't extract
                #     return None,None,None,None
                for ad,E in ad_energy_dict.items():
                    if split_struct.is_isomorphic(ad,save_order=True):
                        Esep += E
                        break
                else:
                    logging.error("no matching adsorbate for {}".format(split_struct.to_smiles()))
                    return None,None,None,None
            
            dE = Ecad - Esep
        except:
            dE = None
    else:
        if ts_pynta_dir is None:
            xyz = info["xyz"]
            ts_path = os.path.split(os.path.split(info["xyz"])[0])[0]
        else:
            ts_path = os.path.join(ts_pynta_dir,os.path.split(os.path.split(os.path.split(info["xyz"])[0])[0])[1])
            xyz = os.path.join(ts_path,os.path.split(os.path.split(info["xyz"])[0])[1],os.path.split(info["xyz"])[1])
        
        ts_info_path = os.path.join(ts_path,"info.json")

        try:
            admol,neighbor_sites,ninds = generate_TS_2D(atoms, ts_info_path, metal, facet, sites, site_adjacency, len(slab), imag_freq_path=os.path.join(d,"vib.0.traj"), 
                                                        max_dist=np.inf, allowed_structure_site_structures=allowed_structure_site_structures,
                                                        keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
        except (SiteOccupationException,TooManyElectronsException, ValueError, FailedFixBondsException) as e:
            return None,None,None,None
        try:
            vibdata = get_vibdata(os.path.join(d,out_file_name+".xyz"),os.path.join(d,vib_file_name+".json"),len(slab))
        except FileNotFoundError:
            return admol,neighbor_sites,ninds,None
        
        Ecad = atoms.get_potential_energy() - slab.get_potential_energy() + vibdata.get_zero_point_energy()
        
        
        ts = read(xyz)
        ts_vibdata = get_vibdata(xyz,os.path.join(os.path.split(xyz)[0],"vib.json_vib.json"),len(slab))
        Ets = ts.get_potential_energy() - slab.get_potential_energy()  + ts_vibdata.get_zero_point_energy()
        
        num_ts_atoms = len(ts) - len(slab)
        
        coatoms = atoms.copy()
        for i in range(num_ts_atoms): #remove ts atoms before we identify coadsorbates
            coatoms.pop(len(slab))
        try:
            coadmol,coneighbor_sites,coninds = generate_adsorbate_2D(coatoms, sites, site_adjacency, len(slab), max_dist=np.inf, cut_off_num=None)
        except (SiteOccupationException,TooManyElectronsException,ValueError,FailedFixBondsException) as e:
            return None,None,None,None
            
        split_structs = split_adsorbed_structures(coadmol)
        Esep = Ets
        for split_struct in split_structs:
            for ad,E in ad_energy_dict.items():
                if split_struct.is_isomorphic(ad,save_order=True):
                    Esep += E
                    break
            else:
                logging.error("no matching adsorbate for {}".format(split_struct.to_smiles()))
                return None,None,None,None
                
        dE = Ecad - Esep
    
    return admol,neighbor_sites,ninds,dE



def check_TS_coadsorbate_disruption(atoms,admol,nslab3D,nslab2D,ntsatoms,mult=1.3):
    ad = atoms[nslab3D:]
    cutoffs = natural_cutoffs(ad,mult=mult)
    adanalysis = Analysis(ad,cutoffs=cutoffs)
    adadj = adanalysis.adjacency_matrix[0]
    for i in range(nslab2D+ntsatoms,len(admol.atoms)):
        at = admol.atoms[i]
        ind = i - nslab2D
        if adadj[ind,:ntsatoms].sum() > 0 or adadj[:ntsatoms,ind].sum() > 0 :
            for k in range(ntsatoms):
                if (not adadj[ind,k]) and (not adadj[k,ind]):
                    continue
                else:
                    at = admol.atoms[nslab2D+k]
                    for bd in admol.get_bonds(at).values():
                        if bd.get_order_str() == 'R':
                            return True
                        
    return False
    
def tagsites(atoms,sites):
    aview = deepcopy(atoms)
    anames = ['He','Ne','Ar','Kr','Xe','Rn']
    for i,site in enumerate(sites):
        add_adsorbate_to_site(aview,Atoms(anames[i], [(0, 0, 0)]), 0, sites[i], height=1.0)
    return aview

def extract_sample(d,ad_energy_dict,slab,metal,facet,sites,site_adjacency,pynta_dir,max_dist=3.0,rxn_alignment_min=0.7,
                   coad_disruption_tol=1.1,out_file_name="out",init_file_name="init",vib_file_name="vib",is_ad=None,
                  use_allowed_site_structures=True):
    out_dict = dict()
    
    atoms_init = read(os.path.join(d,init_file_name+".xyz"))
    if not os.path.exists(os.path.join(d,out_file_name+".xyz")):
        atoms = None
    else:
        atoms = read(os.path.join(d,out_file_name+".xyz"))

    nslab = len(slab)
    
    if len(atoms_init) < len(slab): #gas phase
        return None
    #view(atoms_init)
    #view(atoms)
    # try:
    #     vibdata = get_vibdata(os.path.join(d,out_file_name+".xyz"),os.path.join(d,vib_file_name+".json"),len(slab))
    # except:
    #     vibdata = None
    keep_binding_vdW_bonds=False 
    keep_vdW_surface_bonds=False
    
    try:
        with open(os.path.join(d,"info.json"),'r') as f:
            info = json.load(f)
            mol = Molecule().from_adjacency_list(info["adjlist"],check_consistency=False)
            for bd in mol.get_all_edges():
                if bd.order == 0:
                    if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
                        keep_binding_vdW_bonds = True
                        m = mol.copy(deep=True)
                        b = m.get_bond(m.atoms[mol.atoms.index(bd.atom1)],m.atoms[mol.atoms.index(bd.atom2)])
                        m.remove_bond(b)
                        out = m.split()
                        if len(out) == 1: #vdW bond is not only thing connecting adsorbate to surface
                            keep_vdW_surface_bonds = True
    except FileNotFoundError:
        info = None
        
    if is_ad is None:
        is_ad = os.path.split(os.path.split(os.path.split(os.path.split(info['xyz'])[0])[0])[0])[1] == 'Adsorbates'

    if is_ad and info:
        orig_xyz = os.path.join(pynta_dir,os.sep.join(os.path.normpath(info['xyz']).split(os.sep)[-4:]))
    else:
        orig_xyz = os.path.join(pynta_dir,os.sep.join(os.path.normpath(info['xyz']).split(os.sep)[-3:]))
    
    if atoms is None:
        out_dict["isomorphic"] = False 
        out_dict["init_info"] = info["adjlist"]
        out_dict["out"] = None 
        out_dict["init_extracted"] = None 
        out_dict["sample_dir"] = d
        out_dict["orig_xyz"] = orig_xyz
        out_dict["valid"] = False
        out_dict["dE"] = None
        return out_dict
        
    if info:
        with open(os.path.join(os.path.split(os.path.split(orig_xyz)[0])[0],"info.json"),'r') as f:
            info_clean = json.load(f)
    if orig_xyz:
        if not is_ad:
            reactants_clean = Molecule().from_adjacency_list(info_clean["reactants"])
            cut_off_num = len(atoms) - nslab - len([a for a in reactants_clean.atoms if not a.is_surface_site()])
        else:
            ad_clean = Molecule().from_adjacency_list(info_clean["adjlist"])
            cut_off_num = len(atoms) - nslab - len([a for a in ad_clean.atoms if not a.is_surface_site()])
    else:
        cut_off_num = None

    if use_allowed_site_structures:
        allowed_structure_site_structures = generate_allowed_structure_site_structures(os.path.join(pynta_dir,"Adsorbates"),
                                                                                         sites,site_adjacency,nslab,max_dist=max_dist,
                                                                                                 cut_off_num=None)
    else:
        allowed_structure_site_structures = None
    
    admol,neighbor_sites,ninds,dE = load_coverage_delta(d,ad_energy_dict,slab,metal,facet,sites,site_adjacency,
                                            ts_pynta_dir=pynta_dir,allowed_structure_site_structures=allowed_structure_site_structures,
                                                       out_file_name=out_file_name,vib_file_name=vib_file_name,is_ad=is_ad,
                                                       keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
    
    
    if is_ad:
        try:
            admol_init,neighbor_sites_init,ninds_init = generate_adsorbate_2D(atoms_init, sites, site_adjacency, nslab, 
                     max_dist=None, allowed_structure_site_structures=allowed_structure_site_structures,
                     keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
        except (TooManyElectronsException,FailedFixBondsException,ValueError):
            admol_init = None
            neighbor_sites_init = None
            ninds_init = None
            
        if info:
            admol_info = Molecule().from_adjacency_list(info["adjlist"],check_consistency=False)
            
        valid = admol is not None #and len(admol.split()) == 1

    else:
        valid = True
        info_path = os.path.join(os.path.split(os.path.split(orig_xyz)[0])[0],"info.json")
        imag_freq_path = os.path.join(os.path.split(orig_xyz)[0],"vib.0.traj")

        try:
            admol_init,neighbor_sites_init,ninds_init = generate_TS_2D(atoms_init, info_path, metal, facet, sites, site_adjacency, nslab, 
                     imag_freq_path=imag_freq_path, max_dist=None, 
                    allowed_structure_site_structures=allowed_structure_site_structures,
                    keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
        except (TooManyElectronsException,FailedFixBondsException,ValueError) as e:
            valid = False
            admol_init = None
        except SiteOccupationException as e:
            logging.error("raised SiteOccupationException on init structure")
            valid = False
            admol_init = None

        admol_info = Molecule().from_adjacency_list(info["adjlist"],check_consistency=False)
        
        if valid:
            valid = not (admol is None)
            if not valid:
                logging.error("TS admol has disconnected graph")

        if valid:
            #convert TS path
            ts_path = os.path.join(pynta_dir, os.path.split(os.path.split(os.path.split(orig_xyz)[0])[0])[1], os.path.split(os.path.split(orig_xyz)[0])[1], os.path.split(orig_xyz)[1])
            
            vibdata_ts = get_vibdata(ts_path,os.path.join(os.path.split(ts_path)[0],"vib.json_vib.json"),len(slab))
            with open(os.path.join(os.path.split(os.path.split(ts_path)[0])[0],"info.json"),'r') as f:
                ts_info = json.load(f)
        
            nslab3D = ts_info["nslab"]
            ntsatoms = len([a for a in Molecule().from_adjacency_list(ts_info['reactants']).atoms if not a.is_surface_site()])
            nslab2D = len(ninds)

            if os.path.exists(os.path.join(d,vib_file_name+".json")):
                coad_disrupt = check_TS_coadsorbate_disruption(atoms,admol,nslab3D,nslab2D,ntsatoms,mult=coad_disruption_tol)
                tr_sample = Trajectory(os.path.join(d,"vib.0.traj"))

                #view(tr_sample)
                tr_ts = Trajectory(os.path.join(os.path.split(ts_path)[0],"vib.0.traj"))
                v_ts = tr_ts[15].positions[nslab:] - tr_ts[16].positions[nslab:]
                v_ts /= np.linalg.norm(v_ts)
                v_sample = tr_sample[15].positions[nslab:] - tr_sample[16].positions[nslab:]
                v_sample /= np.linalg.norm(v_sample)
                v_prod = np.abs(np.sum(v_sample[:len(v_ts)]*v_ts))

                if len(admol.split()) > 1 or v_prod < rxn_alignment_min or coad_disrupt:
                    valid = False
                else:
                    valid = True
            else:
                coad_disrupt = None
                valid = None
            
    out_dict["valid"] = valid
    out_dict["dE"] = dE
    if admol is None or admol_info is None:
        out_dict["isomorphic"] = False
    else:
        out_dict["isomorphic"] = admol.is_isomorphic(admol_info,save_order=True)
    if admol:
        out_dict["out"] = admol.to_adjacency_list()
    else:
        out_dict["out"] = None
    if info:
        out_dict["init_info"] = admol_info.to_adjacency_list()
    if admol_init:
        out_dict["init_extracted"] = admol_init.to_adjacency_list()
    else:
        out_dict["init_extracted"] = None
    out_dict['orig_xyz'] = orig_xyz
    out_dict['sample_dir'] = d
    return out_dict

def process_calculation(d,ad_energy_dict,slab,metal,facet,sites,site_adjacency,pynta_dir,coadmol_E_dict,max_dist=3.0,rxn_alignment_min=0.7,
                    coad_disruption_tol=1.1,out_file_name="out",init_file_name="init",vib_file_name="vib",is_ad=None):
    
    datums_stability = []
    datum_E = None
    
    outdict = extract_sample(d,ad_energy_dict,slab,metal,facet,sites,site_adjacency,pynta_dir,max_dist=max_dist,rxn_alignment_min=rxn_alignment_min,
                    coad_disruption_tol=coad_disruption_tol,
                    out_file_name=out_file_name,init_file_name=init_file_name,vib_file_name=vib_file_name,is_ad=is_ad)
    
    if outdict["init_info"]:
        mol_init = Molecule().from_adjacency_list(outdict["init_info"],check_consistency=False)
    else:
        mol_init = Molecule().from_adjacency_list(outdict["init_extracted"],check_consistency=False)  
    
    if (outdict["valid"] is not None) and (not outdict["isomorphic"]):
        datums_stability.append(Datum(mol_init,False))

    if outdict["out"] is not None:
        mol_out = Molecule().from_adjacency_list(outdict["out"],check_consistency=False)
        
        if len(mol_out.split()) == 1: #ignore if desorbed species
            if outdict["valid"] and outdict["dE"] is not None: #there is an energy value
                skip = False
                for at in mol_out.atoms:
                    n = 0
                    if at.is_surface_site():
                        for k,v in mol_out.get_edges(at).items():
                            if not k.is_surface_site():
                                n += 1
                    if n > 1:
                        skip = True
                if not skip:
                    try:
                        datum_E = Datum(mol_out,outdict["dE"]*96.48530749925793*1000.0 - get_atom_centered_correction(mol_out,coadmol_E_dict)) #eV to J/mol
                        datums_stability.append(Datum(mol_out,True))
                    except KeyError:
                        pass
                    
    
    return datum_E,datums_stability

def get_configs_for_calculation(configs_of_concern_by_coad_admol,Ncoad_energy_by_coad_admol,admol_name_structure_dict,coadnames,
                                computed_configs,tree_regressor,Ncalc_per_iter,T=5000.0,calculation_selection_iterations=10):
    group_to_occurence = dict()
    configs_of_concern = []
    for coadname in coadnames:
        for admol_name in configs_of_concern_by_coad_admol[coadname].keys():
            admol = admol_name_structure_dict[admol_name]
            Nocc_isolated = len([a for a in admol.atoms if a.is_surface_site() and any(not a2.is_surface_site() for a2 in a.bonds.keys())])
            configs_of_concern_admol = configs_of_concern_by_coad_admol[coadname][admol_name]
            group_to_occurence_admol = dict()
            N = 0 #number of group contributions associated with given admol
            for v in configs_of_concern_admol:
                m = v[0]
                E = v[1]
                sigma = v[3]
                grps = v[2]
                Nocc = len([a for a in m.atoms if a.is_surface_site() and any(not a2.is_surface_site() for a2 in a.bonds.keys())])
                configs_of_concern.append(v)
                Emin = Ncoad_energy_by_coad_admol[coadname][admol_name][Nocc-Nocc_isolated]
                for grp in grps:
                    if grp in group_to_occurence_admol:
                        group_to_occurence_admol[grp] += np.exp(-(E-Emin)/(8.314*T))
                        N += 1
                    else:
                        group_to_occurence_admol[grp] = np.exp(-(E-Emin)/(8.314*T))
                        N += 1
            for g,n in group_to_occurence_admol.items():
                if grp in group_to_occurence:
                    group_to_occurence[grp] += n/N
                else:
                    group_to_occurence[grp] = n/N 
                    
    concern_groups = list(group_to_occurence.keys()) #selected ordering
        
    group_to_weight = np.array([group_to_occurence[g]*tree_regressor.nodes[g].rule.uncertainty for g in concern_groups])
    
    config_to_group_fract = dict()
    for j,v in enumerate(configs_of_concern):
        config,E,tr,std = v
        config_group_unc = np.array([tr.count(g)*tree_regressor.nodes[g].rule.uncertainty for g in concern_groups])
        s = config_group_unc.sum()
        if s == 0:
            config_to_group_fract[j] = config_group_unc
        else:
            config_to_group_fract[j] = config_group_unc/s
            
    
    configs_for_calculation = []
    group_fract_for_calculation = []
    maxval = 0.0
    config_list = [x[0] for x in configs_of_concern]
    #sort lower numbers of coadsorbates first, larger numbers of coadsorbates later
    sorted_config_list = sorted(config_list[:], key = lambda x: len([bd for bd in x.get_all_edges() if (bd.atom1.is_surface_site() and not bd.atom2.is_surface_site()) or (bd.atom2.is_surface_site() and not bd.atom1.is_surface_site())]))

    for q in range(calculation_selection_iterations): #loop the greedy optimization a number of times to get close to local min
        for config in sorted_config_list:
            ind = [i for i,c in enumerate(config_list) if c is config][0]
            if ind not in config_to_group_fract.keys():
                logging.error("config not in config_to_group_fract")
                continue
            for c in computed_configs:
                if config.is_isomorphic(c,save_order=True):
                    break
            else:
                if len(configs_for_calculation) < Ncalc_per_iter:
                    configs_for_calculation = configs_for_calculation + [config]
                    group_fract = config_to_group_fract[ind]
                    group_fract_for_calculation.append(group_fract)
                    maxval = np.linalg.norm(sum(group_fract_for_calculation) * group_to_weight, ord=1)
                else:
                    group_fract = config_to_group_fract[j]
                    g_old_sum = sum(group_fract_for_calculation)
                    maxarglocal = None
                    maxvallocal = maxval
                    for i in range(Ncalc_per_iter):
                        val = np.linalg.norm((g_old_sum - group_fract_for_calculation[i] + group_fract) * group_to_weight, ord=1)
                        if val > maxvallocal:
                            maxarglocal = i
                            maxvallocal = val
                            
                    if maxarglocal is not None:
                        group_fract_for_calculation[maxarglocal] = group_fract
                        configs_for_calculation[maxarglocal] = config
                        maxval = maxvallocal

    coad_admol_to_config_for_calculation = dict()
    for config in configs_for_calculation:
        found = False
        for coadname in coadnames:
            coad_admol_to_config_for_calculation[coadname] = dict()
            for admol_name,v in configs_of_concern_by_coad_admol[coadname].items():
                if any(x[0] is config for x in v):
                    if admol_name in coad_admol_to_config_for_calculation[coadname].keys():
                        coad_admol_to_config_for_calculation[coadname][admol_name].append(config)
                    else:
                        coad_admol_to_config_for_calculation[coadname][admol_name] = [config]
                    found = True
                    break
            if found:
                break
        else:
            raise ValueError

    return configs_for_calculation,coad_admol_to_config_for_calculation

