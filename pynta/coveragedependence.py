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
from pynta.main import *
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
        ad1s = get_unique_adsorbate_geometries(adpath1,mol1,sites,site_adjacency,atom_to_molecule_surface_atom_map1,
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
        ad2s = get_unique_adsorbate_geometries(adpath2,mol2,sites,site_adjacency,atom_to_molecule_surface_atom_map2,
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
    print('pairdir',pairdir)
    addir = os.path.join(os.path.split(tsdirs[0])[0],"Adsorbates")
    print('addir',addir)
    slabpath = os.path.join(os.path.split(tsdirs[0])[0],"slab.xyz")
    if not os.path.exists(pairdir):
        os.makedirs(pairdir)
    #if pairdir exists, complain
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
        print('namedir',namedir)
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
                print('outdirs',outdirs)
    return outdirs

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

def get_coadsorbate_information(coadnames,ads_dir,neighbor_sites,sites,site_adjacency,nslab,slab):
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
                               nslab)
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


def add_coadsorbate_3D(geo,site,ad,coads,site_stable_parameters,
                        atom_to_molecule_surface_atom_map,infocoad,coad_occ_dict,coad_height_map,coad2D,
                       metal,facet,sites,site_adjacency,nslab):
    if (site["site"],site["morphology"]) not in site_stable_parameters:
        return None
    for i,coad in enumerate(coads):
        #occ = get_occupied_sites(coad,sites,nslab)
        occ = coad_occ_dict[i]
        occsite = None
        for s in occ:
            if (s['bonding_index'] - nslab) in atom_to_molecule_surface_atom_map.keys():
                occsite = s
                break
        else:
            raise ValueError
        if occsite["site"] == site["site"] and occsite["morphology"] == site["morphology"]:
            #assert len(mol2D.get_all_labeled_atoms()) == 0, mol2D.to_adjacency_list()
            surf_ind = list(atom_to_molecule_surface_atom_map.keys())[0]
            h = coad_height_map[i]
            add_adsorbate_to_site(geo, coad[nslab:], surf_ind, site, height=h)
            return geo,coad2D
    else:
        return None,None

def add_coadsorbate_2D(mol2D,site,coad2D,slab,neighbor_sites_2D,site_2D_inds):
    if site_2D_inds:
        ind2D = site_2D_inds[0]
    else:
        return None
    siteatom = mol2D.atoms[ind2D]
    assert siteatom.site == site["site"]
    for a in siteatom.edges.keys():
        if not a.is_surface_site():
            return None
    c = coad2D.get_desorbed_molecules()[0]
    mol2D = mol2D.merge(c)
    ldict = mol2D.get_all_labeled_atoms()
    label = list(ldict.keys())[0]
    catom = list(ldict.values())[0]
    catom.label = ''
    if label == "*1":
        bd = Bond(siteatom,catom,order=1)
        catom.radical_electrons -= 1
    elif label == "*2":
        bd = Bond(siteatom,catom,order=2)
        if catom.radical_electrons >= 2:
            catom.radical_electrons -= 2
        else:
            catom.lone_pairs -= 1
    elif label == "*3":
        bd = Bond(siteatom,catom,order=3)
        if catom.radical_electrons >= 3:
            catom.radical_electrons -= 3
        elif catom.radical_electrons == 1:
            catom.radical_electrons = 0
            catom.lone_pairs -= 1
        elif catom.radical_electrons == 2:
            catom.radical_electrons = 1
            catom.lone_pairs -= 1
    elif label == "*4":
        bd = Bond(siteatom,catom,order=4)
        if catom.radical_electrons >= 4:
            catom.radical_electrons -= 4
        elif catom.radical_electrons == 1:
            catom.lone_pairs -= 2
        elif catom.radical_electrons == 2:
            catom.radical_electrons = 0
            catom.lone_pairs -= 1
        elif catom.radical_electrons == 3:
            catom.radical_electrons = 1
            catom.lone_pairs -= 1
    else:
        raise ValueError
    mol2D.add_bond(bd)
    mol2D.multiplicity = mol2D.get_radical_count() + 1
    mol2D.update_atomtypes()
    mol2D.update_connectivity_values()
    return mol2D

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

def get_best_adsorbate_xyz(adsorbate_path,sites,nslab):
    """
    load the adsorbates associated with the reaction and find the unique optimized
    adsorbate structures for each species
    returns a dictionary mapping each adsorbate name to a list of ase.Atoms objects
    """
    prefixes = os.listdir(adsorbate_path)
    with open(os.path.join(adsorbate_path,"info.json"),"r") as f:
        info = json.load(f)
    ase_to_mol_surface_atom_map = {int(k):int(v) for k,v in info["gratom_to_molecule_surface_atom_map"].items()}
    mol = Molecule().from_adjacency_list(info["adjlist"])
    geoms = []
    for prefix in prefixes:
        path = os.path.join(adsorbate_path,prefix,prefix+".xyz")
        if os.path.exists(path):
            geoms.append(path)
    xyzs = get_unique_sym(geoms)
    adsorbate = None
    min_energy = np.inf
    for xyz in xyzs:
        geo = read(xyz)
        occ = get_occupied_sites(geo,sites,nslab)
        required_surface_inds = set([ind+nslab for ind in ase_to_mol_surface_atom_map.keys()])
        found_surface_inds = set([site["bonding_index"] for site in occ])
        if len(occ) >= len(mol.get_adatoms()) and required_surface_inds.issubset(found_surface_inds):
            if geo.get_potential_energy() < min_energy:
                adsorbate = xyz
                min_energy = geo.get_potential_energy()

    return adsorbate

class CoverageDependence:
    def __init__(self,path,metal,nslab,surrogate_metal,surface_type,repeats,pynta_run_directory,software,software_kwargs,label,sites,site_adjacency,adsorbates=[],transition_states=dict(),coadsorbates=[],
                max_dist=3.0,frozen_layers=2,fmaxopt=0.05,TS_opt_software_kwargs=None,launchpad_path=None,
                 fworker_path=None,queue=False,njobs_queue=0,reset_launchpad=False,queue_adapter_path=None,
                 num_jobs=25):
        self.path = path
        self.metal = metal
        self.repeats = repeats
        self.nslab = np.product(repeats)
        self.pynta_run_directory = pynta_run_directory
        self.pairs_directory = os.path.join(self.path,"pairs")
        self.slab_path = os.path.join(self.pynta_run_directory,"slab.xyz")
        self.adsorbates_path = os.path.join(self.pynta_run_directory,"Adsorbates")
        self.software = software
        self.software_kwargs = software_kwargs
        self.surface_type = surface_type
        self.facet = metal + surface_type
        self.sites = sites
        self.site_adjacency = site_adjacency
        self.surrogate_metal = surrogate_metal
        #self.nslab = nslab
        
        self.software_kwargs_TS = deepcopy(software_kwargs)
        if TS_opt_software_kwargs:
            for key,val in TS_opt_software_kwargs.items():
                self.software_kwargs_TS[key] = val
        
        if adsorbates != []:
            raise ValueError("Not implemented yet")
        self.adsorbates = adsorbates
        self.transition_states = transition_states
        self.coadsorbates = coadsorbates
        self.max_dist = max_dist
        self.frozen_layers = frozen_layers
        self.layers = self.repeats[2]
        self.freeze_ind = int((self.nslab/self.layers)*self.frozen_layers)
        self.fmaxopt = fmaxopt
        self.label = label
        
        if launchpad_path:
            launchpad = LaunchPad.from_file(launchpad_path)
        else:
            launchpad = LaunchPad()

        if reset_launchpad:
            launchpad.reset('', require_password=False)
        self.launchpad = launchpad
        
        self.fworker_path = fworker_path
        if fworker_path:
            self.fworker = FWorker.from_file(fworker_path)
        else:
            self.fworker = FWorker()

        self.queue = queue
        self.fworker = None
        self.qadapter = None
        if fworker_path:
            self.fworker = FWorker.from_file(fworker_path)
        else:
            self.fworker = FWorker()

        self.njobs_queue = njobs_queue
        self.num_jobs = num_jobs

        if queue:
            self.qadapter = load_object_from_file(queue_adapter_path)

        self.reset_launchpad = reset_launchpad  
        self.pairs_fws = []
        self.fws = []
        
    def setup_pairs_calculations(self):
        tsdirs = [os.path.join(self.pynta_run_directory,t) for t in self.transition_states.keys()] #transition state folder (e.g.TS0) and subfolder (e.g. 76)
        outdirs = setup_pair_opts_for_rxns(self.path,tsdirs,self.coadsorbates,self.surrogate_metal,self.surface_type,max_dist=self.max_dist)
        print(outdirs)

        for d in outdirs:

            fwopt = optimize_firework(d,
                            self.software,"weakopt",
                            opt_method="MDMin",opt_kwargs={'dt': 0.05,"trajectory": "weakopt.traj"},software_kwargs=self.software_kwargs,order=0,
                            run_kwargs={"fmax" : 0.5, "steps" : 30},parents=[],
                              constraints=["freeze up to {}".format(self.freeze_ind)],
                            ignore_errors=True, metal=self.metal, facet=self.surface_type, priority=-1)
            fwopt2 = optimize_firework(os.path.join(os.path.split(d)[0],"weakopt.xyz"),
                            self.software,"out",
                            opt_method="QuasiNewton",opt_kwargs={"trajectory": "out.traj"},software_kwargs=self.software_kwargs,order=0,
                            run_kwargs={"fmax" : self.fmaxopt, "steps" : 70},parents=[fwopt],
                              constraints=["freeze up to {}".format(self.freeze_ind)],
                            ignore_errors=True, metal=self.metal, facet=self.surface_type, priority=-1)
                                
            self.pairs_fws.append(fwopt)
            self.pairs_fws.append(fwopt2)

        self.fws.extend(self.pairs_fws)
    
    def setup_sample_calculations(self):
        ts_xyzs = [os.path.join(self.pynta_run_directory,k,v,"opt.xyz") for k,v in self.transition_states.items()]
        ads = self.adsorbates
        for ts in self.transition_states.keys():
            with open(os.path.join(self.pynta_run_directory,ts,"info.json"),"r") as f:
                info = json.load(f)
                for name in info["species_names"]+info["reverse_names"]:
                    if name not in ads:
                        ads.append(name)
        
        ad_xyzs = []
        for ad in ads:
            p = os.path.join(self.pynta_run_directory,"Adsorbates",ad)
            with open(os.path.join(p,"info.json")) as f:
                info = json.load(f)
            mol = Molecule().from_adjacency_list(info["adjlist"])
            if mol.contains_surface_site():
                ad_xyz = get_best_adsorbate_xyz(p,self.sites,self.nslab)
                ad_xyzs.append(ad_xyz)

        outatoms,outmol2Dsad,outmol2Dsts,outxyzsad,outxyzsts=generate_coadsorbed_xyzs(self.path,ad_xyzs,ts_xyzs,
                            self.slab_path,self.pairs_directory,self.adsorbates_path,
                             self.coadsorbates,self.metal,self.surface_type,self.sites,self.site_adjacency,max_dist=self.max_dist)
        
        for xyzad in outxyzsad:
            fwopt = optimize_firework(xyzad,
                            self.software,"weakopt",
                            opt_method="MDMin",opt_kwargs={'dt': 0.05,"trajectory": "weakopt.traj"},software_kwargs=self.software_kwargs,order=0,
                            run_kwargs={"fmax" : 0.5, "steps" : 30},parents=[],
                              constraints=["freeze up to {}".format(self.freeze_ind)],
                            ignore_errors=True, metal=self.metal, facet=self.surface_type, priority=-1)
            fwopt2 = optimize_firework(os.path.join(os.path.split(xyzad)[0],"weakopt.xyz"),
                            self.software,"out",
                            opt_method="QuasiNewton",opt_kwargs={"trajectory": "out.traj"},software_kwargs=self.software_kwargs,order=0,
                            run_kwargs={"fmax" : self.fmaxopt, "steps" : 70},parents=[fwopt],
                              constraints=["freeze up to {}".format(self.freeze_ind)],
                            ignore_errors=True, metal=self.metal, facet=self.surface_type, priority=-1)
            fwvib = vibrations_firework(os.path.join(os.path.split(xyzad)[0],"out.xyz"),self.software,"outvib",
                                software_kwargs=self.software_kwargs,
                                parents=[fwopt2],constraints=["freeze up to {}".format(self.nslab)],socket=False,ignore_errors=True)
            self.fws.append(fwopt)
            self.fws.append(fwopt2)
            self.fws.append(fwvib)
        for xyzts in outxyzsts:
            fwopt = optimize_firework(xyzts,
                            self.software,"out",
                            opt_kwargs={"trajectory": "out.traj"},software_kwargs=self.software_kwargs_TS,order=1,
                            run_kwargs={"fmax" : self.fmaxopt, "steps" : 70},parents=[],
                              constraints=["freeze up to {}".format(self.freeze_ind)],
                            ignore_errors=True, metal=self.metal, facet=self.surface_type, sella=True, priority=-1)
            fwvib = vibrations_firework(os.path.join(os.path.split(xyzad)[0],"out.xyz"),self.software,"outvib",
                                parents=[fwopt],constraints=["freeze up to {}".format(self.nslab)],socket=False,ignore_errors=True)
            self.fws.append(fwopt)
            self.fws.append(fwvib)
        
    def launch(self,single_job=False):
        """
        Call appropriate rapidfire function
        """
        if self.queue:
            rapidfirequeue(self.launchpad,self.fworker,self.qadapter,njobs_queue=self.njobs_queue,nlaunches="infinite")
        elif not self.queue and (self.num_jobs == 1 or single_job):
            rapidfire(self.launchpad,self.fworker,nlaunches="infinite")
        else:
            launch_multiprocess(self.launchpad,self.fworker,"INFO","infinite",self.num_jobs,5)

    def execute(self,run_pairs=True,run_samples=False):
        if run_pairs:
            self.setup_pairs_calculations()
        if run_samples:
            self.setup_sample_calculations()
        """
        generate and launch a Pynta Fireworks Workflow
        """
        wf = Workflow(self.fws, name=self.label)
        self.launchpad.add_wf(wf)

        self.launch()