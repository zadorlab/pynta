from molecule.molecule import Molecule,Atom,Bond
from molecule.exceptions import AtomTypeError
from ase.io import read, write
from ase.geometry import get_distances
from ase.visualize import view
from acat.adsorption_sites import SlabAdsorptionSites
from pynta.utils import get_unique_sym, get_occupied_sites, sites_match
from pynta.mol import *
from pynta.geometricanalysis import *
from pynta.tasks import *
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

def get_unstable_pairs(pairsdir,adsorbate_dir,sites,site_adjacency,nslab,max_dist=3.0,show=False):
    out_pairs = []
    coadname_dict = {"O=[Pt]": 1, "N#[Pt]": 1, "O[Pt]": 2, "[Pt]": 1}
    allowed_structure_site_structure_map = generate_allowed_structure_site_structure_map(adsorbate_dir,sites,site_adjacency,nslab,max_dist=max_dist,cut_coads_off_num=None)
    if show:
        config_show = []
    for pair in os.listdir(pairsdir):
        if not "_" in pair or pair[0] == ".":
            continue
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
                                       cut_coads_off_num=coadname_dict[coadname],allowed_structure_site_structure_map=allowed_structure_site_structure_map)
                    
                
                #g.update(sort_atoms=False)
                outpath = os.path.join(p,"out.xyz")
                if not os.path.exists(outpath):
                    if show:
                        pass
                        
                    out_pairs.append(g.to_group())
                else:
                    final = read(outpath)
                    try:
                        gout = extract_pair_graph(final,sites,site_adjacency,nslab,max_dist=max_dist,allowed_structure_site_structure_map=allowed_structure_site_structure_map)
                        if len(gout.atoms) != len(g.atoms):
                            out_pairs.append(g.to_group())
                            if show:
                                config_show.append(init)
                                config_show.append(final)
                            continue
                        #gout.update(sort_atoms=False)
                    except FindingPathError:
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

def generate_pair_geometries(adpath1,adpath2,slabpath,metal,facet,adinfo1=None,adinfo2=None,
                             max_dist=3.0,imag_freq_max=150.0,symmetric=None):
    """
    adpath1 can be bidentate
    adpath2 must be monodentate
    """
    #slab information
    slab = read(slabpath)
    nslab = len(slab)
    cas = SlabAdsorptionSites(slab,facet,allow_6fold=False,composition_effect=False,
                            label_sites=True,
                            surrogate_metal=metal)
    sites = cas.get_sites()
    site_adjacency = cas.get_neighbor_site_list()
    
    if symmetric is None:
        symmetric = adpath1 == adpath2
    
    #extract information about adsorbates and valid adsorbate geometries
    if os.path.isfile(adpath1):
        raise ValueError
    else:
        if adinfo1 is None:
            adinfo1 = os.path.join(adpath1,"info.json")
        
        with open(adinfo1,"r") as f:
            info1 = json.load(f)
        
        mol1 = Molecule().from_adjacency_list(info1["adjlist"])
        atom_to_molecule_surface_atom_map1 = { int(key):int(val) for key,val in info1["gratom_to_molecule_surface_atom_map"].items()}
        ad1s = get_unique_adsorbate_geometries(adpath1,mol1,cas,atom_to_molecule_surface_atom_map1,
                                    nslab,imag_freq_max=imag_freq_max)
        ad12Ds = []
        ad12Dneighbors = []
        ad12Dninds = []
        for a in ad1s:
            admol,neighbor_sites,ninds = generate_adsorbate_2D(a, sites, site_adjacency, nslab, max_dist=max_dist)
            ad12Ds.append(admol)
            ad12Dneighbors.append(neighbor_sites)
            ad12Dninds.append(ninds)
        
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
        ad2s = get_unique_adsorbate_geometries(adpath2,mol2,cas,atom_to_molecule_surface_atom_map2,
                                    nslab,imag_freq_max=imag_freq_max)
    

    #generate pairs
    pairs = []
    ad1_to_ad1_sites = dict()
    ad1_to_ad2_sites = dict()
    ad1_to_ad2_heights = dict()
    ad1_to_actual_occ = dict()
    for j,ad1 in enumerate(ad1s):
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
        
        
        stable_neighbor_sites = []
        heights = []
        for site in neighbor_sites:
            for i,ad2 in enumerate(ad2s):
                if i not in ad1_to_ad2_sites[j].keys():
                    ad1_to_ad2_sites[j][i] = []
                    ad1_to_ad2_heights[j][i] = []
        
                cas = SlabAdsorptionSites(ad2,facet,allow_6fold=False,composition_effect=False,
                                    label_sites=True,
                                    surrogate_metal=metal)
                sites2 = cas.get_sites()
                occ = get_occupied_sites(ad2,sites2,nslab)
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
    
    if not symmetric:
        for j in range(len(ad1s)):
            ad1_sites = ad1_to_ad1_sites[j]
            
            ad2_sites = []
            ad2_geoms = []
            heights = []
            for i,sites in ad1_to_ad2_sites[j].items():
                for k,site in enumerate(sites):
                    heights.append(ad1_to_ad2_heights[j][i][k])
                    ad2_sites.append(ad1_to_ad2_sites[j][i][k])
                    ad2_geoms.append(ad2s[i])


            inds = get_unique_site_inds(ad2_sites,slab,fixed_point=ad1_sites[0]["position"])

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
                adpairs.append(ad)
                pairmols.append(amol)
    else: #symmetric case, monodentate-monodentate since ad2 must be monodentate
        ad2_site_pairs = []
        ad2_geoms = []
        heights = []
        ad1_inds = []
        for j in range(len(ad1s)):
            ad1_sites = ad1_to_ad1_sites[j]
            
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
            adpairs.append(ad)
            pairmols.append(amol)
    

    return adpairs,pairmols

def get_unique_site_inds(sites,slab,fixed_point=None,tol=0.15):
    fingerprints = []
    for k,site in enumerate(sites):
        if fixed_point is None:
            fingerprints.append((site["morphology"],site["site"]))
        else:
            bd,d = get_distances([site["position"]], [fixed_point], cell=slab.cell, pbc=(True,True,False))
            xydist = np.linalg.norm(bd[0][0][:2])
            zdist = bd[0][0][2]
            fingerprints.append((site["morphology"],site["site"],xydist,zdist,))
    
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

def setup_pair_opts_for_rxns(targetdir,tsdirs,coadnames,metal,facet,max_dist=3.0,imag_freq_max=150.0):
    pairdir = os.path.join(targetdir,"pairs")
    addir = os.path.join(os.path.split(tsdirs[0])[0],"Adsorbates")
    slabpath = os.path.join(os.path.split(tsdirs[0])[0],"slab.xyz")
    if not os.path.exists(pairdir):
        os.makedirs(pairdir)
    
    ads = set()  
    
    for tsdir in tsdirs:
        with open(os.path.join(tsdir,"info.json"),"r") as f:
            info = json.load(f)
        for molname in info["species_names"]+info["reverse_names"]:
            with open(os.path.join(addir,molname,"info.json"),"r") as f:
                molinfo = json.load(f)
            m = Molecule().from_adjacency_list(molinfo["adjlist"])
            if m.contains_surface_site():
                ads.add(molname)
    
    combs = []
    for adname in ads:
        for coadname in coadnames:
            tp = (adname,coadname)
            revtp = (coadname,adname)
            if (revtp not in combs) and (tp not in combs):
                combs.append(tp)
                
    outdirs = []
    for s in combs:
        name = "_".join(s)
        namedir = os.path.join(pairdir,name)
        if not os.path.exists(namedir):
            os.makedirs(namedir)
            ds = [os.path.join(addir,x) for x in s]
            pairs,pairmols = generate_pair_geometries(ds[0],ds[1],slabpath,metal,facet,
                                 max_dist=max_dist,imag_freq_max=imag_freq_max)
            for i,pair in enumerate(pairs):
                os.makedirs(os.path.join(namedir,str(i)))
                write(os.path.join(namedir,str(i),"init.xyz"), pair)
                moldict = {"adjlist": pairmols[i].to_adjacency_list()}
                with open(os.path.join(namedir,str(i),"info.json"),'w') as f:
                    json.dump(moldict,f)
                outdirs.append(os.path.join(namedir,str(i),"init.xyz"))
    return outdirs
