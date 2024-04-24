from molecule.molecule import Molecule,Atom,Bond
from ase.io import read, write
from ase.geometry import get_distances
from ase.geometry.analysis import Analysis
from ase.io.trajectory import Trajectory
from ase.vibrations import VibrationsData
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from pynta.utils import get_unique_sym, get_occupied_sites, sites_match, SiteOccupationException
from pynta.mol import *
import numpy as np
import json
import shutil
import os
import logging

def generate_site_molecule(slab, target_sites, sites, site_adjacency, max_dist=None):
    
    #find overall neighboring sites
    if max_dist:
        neighbor_sites = []
        ninds = []
        for i,site in enumerate(sites):
            for target_site in target_sites:
                v,dist = get_distances([site["position"]], [target_site["position"]], cell=slab.cell, pbc=slab.pbc)

                if np.linalg.norm(v[:2]) < max_dist:
                    neighbor_sites.append(site)
                    ninds.append(i)
                    break
    else:
        neighbor_sites = sites[:]
        ninds = list(range(len(sites)))
    
    atoms = []
    
    for nsite in neighbor_sites:
        atoms.append(Atom(element="X", lone_pairs=0, site=nsite["site"], morphology=nsite["morphology"]))
    
    mol = Molecule(atoms=atoms,multiplicity=1)
    
    
    for i,nind in enumerate(ninds):
        adjinds = site_adjacency[nind]
        for adjind in adjinds:
            if adjind in ninds:
                j = ninds.index(adjind)
                if not mol.has_bond(mol.atoms[i],mol.atoms[j]):
                    mol.add_bond(Bond(mol.atoms[i],mol.atoms[j],1.0))
    
    return mol,neighbor_sites,ninds


def generate_adsorbate_molecule(adslab, sites, site_adjacency, nslab, max_dist=None, allowed_site_dict=dict()):
    #adsorbate
    ad = adslab[nslab:]
    adanalysis = Analysis(ad)
    adadj = adanalysis.adjacency_matrix[0]
    adatoms = []
    for at in ad:
        sym = at.symbol
        if sym in ['N','P']:
            adatoms.append(Atom(element=sym, lone_pairs=1))
        elif sym in ['O','S']:
            adatoms.append(Atom(element=sym, lone_pairs=2))
        elif sym in ["F","Cl","Br","I"]:
            adatoms.append(Atom(element=sym, lone_pairs=3))
        elif sym in ["Si","C","H","Li"]:
            adatoms.append(Atom(element=sym, lone_pairs=0))
        else:
            raise ValueError("Cannot handle atom {}".format(sym))
            
    admol = Molecule(atoms=adatoms)
    for i in range(len(adatoms)):
        for j in range(len(adatoms)):
            if i == j:
                continue
            if adadj[i,j] and not admol.has_bond(admol.atoms[i],admol.atoms[j]):
                admol.add_bond(Bond(admol.atoms[i],admol.atoms[j],1.0))
    
    #slab
    occ = get_occupied_sites(adslab,sites,nslab,allowed_site_dict=allowed_site_dict)
    target_sites = [site for site in sites if any((oc["position"] == site["position"]).all() for oc in occ)]
    
    #find overall neighboring sites
    if max_dist:
        neighbor_sites = []
        ninds = []
        for i,site in enumerate(sites):
            for target_site in target_sites:
                v,dist = get_distances([site["position"]], [target_site["position"]], cell=adslab.cell, pbc=adslab.pbc)
                if np.linalg.norm(v[:2]) < max_dist:
                    neighbor_sites.append(site)
                    ninds.append(i)
                    break
    else:
        neighbor_sites = sites[:]
        ninds = list(range(len(sites)))
        
    atoms = []
    
    for nsite in neighbor_sites:
        atoms.append(Atom(element="X", lone_pairs=0, site=nsite["site"], morphology=nsite["morphology"]))
    
    mol = Molecule(atoms=atoms,multiplicity=1)
    
    
    for i,nind in enumerate(ninds):
        adjinds = site_adjacency[nind]
        for adjind in adjinds:
            if adjind in ninds:
                j = ninds.index(adjind)
                if not mol.has_bond(mol.atoms[i],mol.atoms[j]):
                    mol.add_bond(Bond(mol.atoms[i],mol.atoms[j],1.0))
    
    #merge adsorbate and slab
    fullmol = Molecule(atoms=mol.atoms+admol.atoms)
    fullmol.multiplicity = 1
    
    for oc in occ:
        adind = oc["bonding_index"] - nslab + len(neighbor_sites)
        nind = [ninds[i] for i,nsite in enumerate(neighbor_sites) if (nsite["position"]==oc["position"]).all()][0]
        site_index = ninds.index(nind)
        fullmol.add_bond(Bond(fullmol.atoms[site_index],fullmol.atoms[adind],0))
    
    return fullmol,neighbor_sites,ninds

def fix_bond_orders(mol,allow_failure=False,cleanup_surface_bonds=True):
    atom_order = ["C","Si","N","P","O","S","F","Cl","Br","I","Li","H"]#["H","Li","I","Br","Cl","F","S","O","P","N","Si","C"]#["C","Si","N","P","O","S","F","Cl","Br","I","Li","H"]
    for target in atom_order:
        for a in mol.atoms:
            if a.symbol == target:
                fix_atom(mol,a,allow_failure=allow_failure,cleanup_surface_bonds=cleanup_surface_bonds)
    return None

def get_octet_deviation(atom):
    if atom.symbol in ["Li","H"]:
        return 2 - (atom.lone_pairs*2 + 2*sum(x.order for x in atom.bonds.values() if x.order >= 0.5))
    else:
        return 8 - (atom.lone_pairs*2 + 2*sum(x.order for x in atom.bonds.values() if x.order >= 0.5))
    

def choose_bond_to_fix(symbols,octets,orders,bonded_to_surface):
    atom_order = [["X"],["C","Si"],["N","P"],["O","S"]]
    priority = []
    for i in range(len(symbols)):
        sym = symbols[i]
        if (sym != "X" and octets[i] <= 0) or orders[i] == "R":
            priority.append(0)
        else:
            for j,atoms in enumerate(atom_order):
                if sym in atoms:
                    if not bonded_to_surface[i]:
                        priority.append(j+11)
                    else:
                        priority.append(j+1)
    return np.argmax(priority),np.max(priority)

class TooManyElectronsException(Exception):
    pass

def fix_atom(mol,atom,allow_failure=False,cleanup_surface_bonds=True):
    delta = get_octet_deviation(atom)
    
    if delta < 0: #too many electrons
        logging.error(mol.to_adjacency_list())
        raise TooManyElectronsException("Cannot solve problem of too many electrons not counting surface bonds and all ad bonds are 1")
    elif delta > 0: #too few electrons
        atoms = [k for k,v in atom.bonds.items()]
        symbols = [k.symbol for k in atoms]
        octets = [get_octet_deviation(k) for k in atoms]
        orders = [v.get_order_str() for k,v in atom.bonds.items()]
        bonded_to_surface = [a.is_bonded_to_surface() for a in atoms]
        deviation = True
        while delta > 0:
            ind,val = choose_bond_to_fix(symbols,octets,orders,bonded_to_surface)
            if val == 0:
                if allow_failure:
                    break
                else:
                    raise ValueError("Could not fix bonds")
            bd = atom.bonds[atoms[ind]]
            bd.increment_order()
            delta -= 2
            octets[ind] -= 2
        
    #at end clean up improper surface bonds
    if cleanup_surface_bonds:
        to_remove = []
        for k,v in atom.bonds.items():
            if k.symbol == 'X' and v.order == 0:    
                to_remove.append(v)
        
        for v in to_remove:
            mol.remove_bond(v)
        
def generate_adsorbate_2D(atoms, sites, site_adjacency, nslab, max_dist=3.0, cut_multidentate_off_num=None, allowed_structure_site_structures=None):
    admol,neighbor_sites,ninds = generate_adsorbate_molecule(atoms, sites, site_adjacency, nslab, max_dist=max_dist)
    
    if cut_multidentate_off_num:
        bds_to_remove = []
        for i in range(len(admol.atoms)-cut_multidentate_off_num,len(admol.atoms)):
            if admol.atoms[i].is_bonded_to_surface():         
                for a,b in admol.atoms[i].edges.items():
                    if not a.is_surface_site() and (a.is_bonded_to_surface() or admol.atoms.index(a) < len(admol.atoms)-cut_multidentate_off_num):
                        if b not in bds_to_remove:
                            bds_to_remove.append(b)
        for b in bds_to_remove:
            admol.remove_bond(b)
    
    fix_bond_orders(admol)
    
    if allowed_structure_site_structures: #we are going to analyze the initial occupational analysis and update it based on expectations
        allowed_site_dict = dict()
        slabless = remove_slab(admol)
        for site_structs in allowed_structure_site_structures:
            struct = generate_without_site_info(site_structs[0])
            grp_struct = struct.to_group()
            mappings = slabless.find_subgraph_isomorphisms(grp_struct,save_order=True)
            considered_sites = []
            for mapping in mappings:
                atouts = [at for at in mapping.keys() if at.is_surface_site()]
                if set(atouts) in considered_sites:
                    continue
                else:
                    considered_sites.append(set(atouts))
                
                subgraph,inds = pluck_subgraph(slabless,atouts[0])
                for site_struct in site_structs:
                    if subgraph.is_isomorphic(site_struct,save_order=True):
                        break
                else:
                    for atout in atouts:
                        adatom = [a for a in atout.bonds.keys() if not a.is_surface_site()][0]
                        ind = slabless.atoms.index(adatom) - len([a for a in slabless.atoms if a.is_surface_site()])
                        struct_ind = [grp_struct.atoms.index(a) for aout,a in mapping.items() if aout==atout][0]
                        sitetype = [(site_struct.atoms[struct_ind].site,site_struct.atoms[struct_ind].morphology) for site_struct in site_structs]
                        if ind+nslab in allowed_site_dict.keys():
                            allowed_site_dict[ind+nslab].extend(sitetype)
                        else:
                            allowed_site_dict[ind+nslab] = sitetype

        admol,neighbor_sites,ninds = generate_adsorbate_molecule(atoms, sites, site_adjacency, nslab, max_dist=max_dist, allowed_site_dict=allowed_site_dict)

        if cut_multidentate_off_num:
            bds_to_remove = []
            for i in range(len(admol.atoms)-cut_multidentate_off_num,len(admol.atoms)):
                if admol.atoms[i].is_bonded_to_surface():         
                    for a,b in admol.atoms[i].edges.items():
                        if not a.is_surface_site() and (a.is_bonded_to_surface() or admol.atoms.index(a) < len(admol.atoms)-cut_multidentate_off_num):
                            if b not in bds_to_remove:
                                bds_to_remove.append(b)
            for b in bds_to_remove:
                admol.remove_bond(b)
            
        fix_bond_orders(admol)
    
    admol.update_atomtypes()
    admol.update_connectivity_values()
    
    return admol,neighbor_sites,ninds

def generate_TS_2D(atoms, info_path,  metal, facet, sites, site_adjacency, nslab, imag_freq_path=None,
                     max_dist=3.0, cut_multidentate_off_num=None, allowed_structure_site_structures=None,
                     site_bond_cutoff=3.0):

    admol,neighbor_sites,ninds = generate_adsorbate_molecule(atoms, sites, site_adjacency, 
                                                             nslab, max_dist=max_dist) 

    if cut_multidentate_off_num:
        bds_to_remove = []
        for i in range(len(admol.atoms)-cut_multidentate_off_num,len(admol.atoms)):
            if admol.atoms[i].is_bonded_to_surface():         
                for a,b in admol.atoms[i].edges.items():
                    if not a.is_surface_site() and (a.is_bonded_to_surface() or admol.atoms.index(a) < len(admol.atoms)-cut_multidentate_off_num):
                        if b not in bds_to_remove:
                            bds_to_remove.append(b)
        for b in bds_to_remove:
            admol.remove_bond(b)
    
    with open(info_path) as f:
        info = json.load(f)
    
    occ = get_occupied_sites(atoms,sites,nslab,site_bond_cutoff=site_bond_cutoff)
    
    template_mol_map = [{int(k):v for k,v in x.items()}  for x in info["template_mol_map"]]
    molecule_to_atom_maps = [{int(k):v for k,v in x.items()}  for x in info["molecule_to_atom_maps"]]
    broken_bonds,formed_bonds = get_broken_formed_bonds(Molecule().from_adjacency_list(info["reactants"]),
                                                        Molecule().from_adjacency_list(info["products"]))
    rbonds = broken_bonds | formed_bonds
    if info["forward"]:
        template = Molecule().from_adjacency_list(info["reactants"])
    else:
        template = Molecule().from_adjacency_list(info["products"])
    
    for bd in rbonds:
        inds2D = []
        for label in bd:
            a = template.get_labeled_atoms(label)[0]
            aind = template.atoms.index(a)
            aseind = get_ase_index(aind,template_mol_map,molecule_to_atom_maps,nslab,info["ads_sizes"])
            if aseind:
                inds2D.append(aseind - nslab + len(neighbor_sites))

        if len(inds2D) == 1:
            i = inds2D[0]
            for batom in admol.atoms[i].bonds.keys(): #try to find a surface bond
                if batom.is_surface_site():
                    inds2D.append(admol.atoms.index(batom))
                    break
            else: #didn't find surface bond, need to identify the site
                pos = None
                for site in occ:
                    if site['bonding_index'] == i - len(neighbor_sites) + nslab:
                        pos = site["position"]
                        break
                else:
                    logging.error("couldn't find right site bond")
                    logging.error(admol.to_adjacency_list())
                    logging.error(i - len(neighbor_sites))
                    raise SiteOccupationException
                for j,s in enumerate(neighbor_sites):
                    if (s["position"] == pos).all():
                        inds2D.append(j)
                        break
                else:
                    raise IndexError("could not find matching site")

        if admol.has_bond(admol.atoms[inds2D[0]],admol.atoms[inds2D[1]]):
            b = admol.get_bond(admol.atoms[inds2D[0]],admol.atoms[inds2D[1]])
            if b.get_order_str() != "R":
                b.set_order_str("R")
            else: 
                b.set_order_str("S") #this ensures the the split_ts leaves this bond connected 
        else:
            admol.add_bond(Bond(admol.atoms[inds2D[0]],admol.atoms[inds2D[1]],"R"))

    fix_bond_orders(admol,allow_failure=True)

    
    if allowed_structure_site_structures: #we are going to analyze the initial occupational analysis and update it based on expectations
        allowed_site_dicts = []
        restricted_inds = set()
        rs = split_ts_to_reactants(admol) #get reactants/products
        for r in rs:
            allowed_site_dict = dict()
            slabless = remove_slab(r) #take the individual component
            for site_structs in allowed_structure_site_structures:
                struct = generate_without_site_info(site_structs[0])
                grp_struct = struct.to_group() #known stable adsorbates
                mappings = slabless.find_subgraph_isomorphisms(grp_struct,save_order=True) #try to find mappings to known stable structures 
                for mapping in mappings: #go through each mapping
                    atouts = [at for at in mapping.keys() if at.is_surface_site()]
                    subgraph,inds = pluck_subgraph(slabless,atouts[0]) #extract the matching subgraph
                    for site_struct in site_structs: #go through the known stable structures
                        if subgraph.is_isomorphic(site_struct,save_order=True):
                            break
                    else:
                        
                        for atout in atouts:
                            adatom = [a for a in atout.bonds.keys() if not a.is_surface_site()][0]
                            ind = slabless.atoms.index(adatom) - len([a for a in slabless.atoms if a.is_surface_site()])
                            struct_ind = [grp_struct.atoms.index(a) for aout,a in mapping.items() if aout==atout][0]
                            assert site_structs[0].atoms[struct_ind].is_surface_site()
                            sitetype = [(site_struct.atoms[struct_ind].site,site_struct.atoms[struct_ind].morphology) for site_struct in site_structs]
                            if ind+nslab in allowed_site_dict.keys():
                                allowed_site_dict[ind+nslab].extend(sitetype)
                            else:
                                allowed_site_dict[ind+nslab] = sitetype

            allowed_site_dicts.append(allowed_site_dict)
            restricted_inds |= set(allowed_site_dict.keys())
            
        allowed_site_dict = dict()
        for ind in restricted_inds:
            for d in allowed_site_dicts:
                if ind in d.keys():
                    if ind not in allowed_site_dict.keys():
                        allowed_site_dict[ind] = set(d[ind])
                    else:
                        allowed_site_dict[ind] &= set(d[ind])

        allowed_site_dict = {k: list(v) for k,v in allowed_site_dict.items()}
    else:
        allowed_site_dict = dict()

    admol,neighbor_sites,ninds = generate_adsorbate_molecule(atoms, sites, site_adjacency, 
                                                             nslab, max_dist=max_dist, allowed_site_dict=allowed_site_dict) 

    if cut_multidentate_off_num:
        bds_to_remove = []
        for i in range(len(admol.atoms)-cut_multidentate_off_num,len(admol.atoms)):
            if admol.atoms[i].is_bonded_to_surface():         
                for a,b in admol.atoms[i].edges.items():
                    if not a.is_surface_site() and (a.is_bonded_to_surface() or admol.atoms.index(a) < len(admol.atoms)-cut_multidentate_off_num):
                        if b not in bds_to_remove:
                            bds_to_remove.append(b)
        for b in bds_to_remove:
            admol.remove_bond(b)

    for bd in rbonds:
        inds2D = []
        for label in bd:
            a = template.get_labeled_atoms(label)[0]
            aind = template.atoms.index(a)
            aseind = get_ase_index(aind,template_mol_map,molecule_to_atom_maps,nslab,info["ads_sizes"])
            if aseind:
                inds2D.append(aseind - nslab + len(neighbor_sites))

        if len(inds2D) == 1:
            i = inds2D[0]
            for batom in admol.atoms[i].bonds.keys(): #try to find a surface bond
                if batom.is_surface_site():
                    inds2D.append(admol.atoms.index(batom))
                    break
            else: #didn't find surface bond, need to identify the site
                pos = None
                for site in occ:
                    if site['bonding_index'] == i - len(neighbor_sites) + nslab:
                        pos = site["position"]
                        break
                else:
                    raise ValueError("didn't find occupied site")
                for j,s in enumerate(neighbor_sites):
                    if (s["position"] == pos).all():
                        inds2D.append(j)
                        break
                else:
                    raise IndexError("could not find matching site")

        if admol.has_bond(admol.atoms[inds2D[0]],admol.atoms[inds2D[1]]):
            b = admol.get_bond(admol.atoms[inds2D[0]],admol.atoms[inds2D[1]])
            if b.get_order_str() != "R":
                b.set_order_str("R")
            else: #diffusion on surface results in two R bonds to different sites
                if imag_freq_path is None:
                    logging.warning("Could not resolve diffusion type R-bonds due to lack of imaginary frequency info")
                    continue
                if admol.atoms[inds2D[0]].is_surface_site():
                    site_atom_center = admol.atoms[inds2D[0]]
                    Ratom = admol.atoms[inds2D[1]]
                    if allowed_structure_site_structures:
                        valid_sites = allowed_site_dict[inds2D[1] - len(neighbor_sites) + nslab]
                else:
                    site_atom_center = admol.atoms[inds2D[1]]
                    Ratom = admol.atoms[inds2D[0]]
                    if allowed_structure_site_structures:
                        valid_sites = allowed_site_dict[inds2D[0] - len(neighbor_sites) + nslab]
                site_center_ind = admol.atoms.index(site_atom_center)
                
                #get imaginary frequency motion
                tr_ts = Trajectory(imag_freq_path)
                v_ts,_ = get_distances(tr_ts[15].positions[admol.atoms.index(Ratom) - len(neighbor_sites) + nslab], tr_ts[16].positions[admol.atoms.index(Ratom) - len(neighbor_sites) + nslab], cell=atoms.cell, pbc=atoms.pbc)
                v_ts = v_ts[0,0,:]
                v_ts /= np.linalg.norm(v_ts)
                #position of atom
                target_pos = atoms.positions[admol.atoms.index(Ratom) - len(neighbor_sites) + nslab]
                site_pair = None
                best_metric = None
                iters = 1
                for i,site1 in enumerate(neighbor_sites):
                    pos1 = site1["position"]
                    site1_id = (site1["site"],site1["morphology"])
                    for j in range(i):
                        site2 = neighbor_sites[j]
                        if i == j:
                            continue
                        site2_id = (site2["site"],site2["morphology"])
                        pos2 = site2["position"]
                        v,dist = get_distances([pos1],[pos2],cell=atoms.cell,pbc=atoms.pbc)
                        dist = dist[0][0]
                        v = v[0,0,:]
                        v /= np.linalg.norm(v)
                        motion_alignment = abs(np.dot(v_ts,v))
                        v1t,dist1t = get_distances([pos1],[target_pos],cell=atoms.cell,pbc=atoms.pbc)
                        dist1t = dist1t[0][0]
                        v1t = v1t[0,0,:]
                        v2t,dist2t = get_distances([pos2],[target_pos],cell=atoms.cell,pbc=atoms.pbc)
                        v2t = v2t[0,0,:]
                        dist2t = dist2t[0][0]
                        halfwayness = np.linalg.norm(v1t+v2t)/(dist2t+dist1t)
                        dist_site = max(dist1t,dist2t)
                        if not np.isnan(motion_alignment) and motion_alignment > 0.95 and dist < 3.0 and dist_site < 2.0 and (allowed_structure_site_structures and site1_id in valid_sites and site2_id in valid_sites):
                            iters += 1
                            metric = 1.0 / (halfwayness*dist_site)
                            if site_pair is None:
                                site_pair = [i,j]
                                best_metric = metric
                            elif metric > best_metric:
                                site_pair = [i,j]
                                best_metric = metric
                if site_pair is None:
                    raise ValueError("Imaginary frequency doesn't mean minimal constraints for proper TS motion")
                admol.remove_bond(b)
                admol.add_bond(Bond(admol.atoms[site_pair[0]],Ratom,"R"))
                admol.add_bond(Bond(admol.atoms[site_pair[1]],Ratom,"R")) 
        else:
            admol.add_bond(Bond(admol.atoms[inds2D[0]],admol.atoms[inds2D[1]],"R"))

    fix_bond_orders(admol,allow_failure=True) #Reaction bonds prevent proper octet completion
    
    #remove surface bonds that abused for octet completion in reactions
    for atom in admol.atoms:
        surfbds = [bd for at,bd in atom.bonds.items() if at.is_surface_site() and bd.get_order_str() != 'R']
        if len(surfbds) == 1:
            Nrxnbd = len([bd for at,bd in atom.bonds.items() if bd.get_order_str() == 'R' and not at.is_surface_site()])
            if Nrxnbd >= 2 and abs(Nrxnbd/2 - surfbds[0].order) < 1e-3:
                admol.remove_edge(surfbds[0])
        
    admol.update_atomtypes()
    admol.update_connectivity_values()
    for i in range(len(neighbor_sites)):
        assert neighbor_sites[i]["site"] == admol.atoms[i].site
        assert neighbor_sites[i]["morphology"] == admol.atoms[i].morphology
    return admol,neighbor_sites,ninds

def tagsites(atoms,sites):
    aview = deepcopy(atoms)
    anames = ['He','Ne','Ar','Kr','Xe','Rn']
    for i,site in enumerate(sites):
        add_adsorbate_to_site(aview,Atoms(anames[i], [(0, 0, 0)]), 0, sites[i], height=1.0)
    return aview

def split_ts_to_reactants(ts2d,tagatoms=False):
    rs = []
    rbondnum = 0
    for bd in ts2d.get_all_edges():
        if bd.is_reaction_bond():
            rbondnum += 1
    
    combs = list(itertools.product([0,1],repeat=rbondnum))

    for comb in combs:
        r = deepcopy(ts2d)
        iters = 0
        for i,bd in enumerate(r.get_all_edges()):
            if bd.is_reaction_bond():
                if tagatoms:
                    if not bd.atom1.is_surface_site() and bd.atom1.is_bonded_to_surface():
                        bd.atom1.label = "*"
                    if not bd.atom2.is_surface_site() and bd.atom2.is_bonded_to_surface():
                        bd.atom2.label = "*"
                    
                if comb[iters] == 0:
                    r.remove_bond(bd)
                else:
                    bd.set_order_str('S')
                iters += 1
            else:
                bd.set_order_str('S')
        try:
            fix_bond_orders(r)
            r.update(sort_atoms=False)
            r.update_connectivity_values()
            rs.append(r)
        except (TooManyElectronsException,ValueError):
            pass
    
    return rs
      
def generate_allowed_structure_site_structures(adsorbate_dir,sites,site_adjacency,nslab,max_dist=3.0,cut_multidentate_off_num=None):
    allowed_structure_site_structures = []
    for ad in os.listdir(adsorbate_dir):
        if not os.path.isdir(os.path.join(adsorbate_dir,ad)):
            continue
        with open(os.path.join(adsorbate_dir,ad,"info.json"),'r') as f:
            info = json.load(f)
            target_mol = Molecule().from_adjacency_list(info["adjlist"])
            atom_to_molecule_surface_atom_map = {int(k): int(v) for k,v in info["gratom_to_molecule_surface_atom_map"].items()}
        if not target_mol.contains_surface_site():
            continue
        geoms = get_adsorbate_geometries(os.path.join(adsorbate_dir,ad),target_mol,sites,atom_to_molecule_surface_atom_map,nslab)
        mols = []
        for geom in geoms:
            admol,neighbor_sites,ninds = generate_adsorbate_2D(geom, sites, site_adjacency, nslab, max_dist=max_dist, cut_multidentate_off_num=cut_multidentate_off_num)
            m = remove_slab(admol)
            if not generate_without_site_info(m).is_isomorphic(target_mol,save_order=True): #geom didn't optimize to target
                continue
            for mol in mols:
                if mol.is_isomorphic(m,save_order=True):
                    break
            else:
                m.update_atomtypes()
                m.update_connectivity_values()
                mols.append(m)
        allowed_structure_site_structures.append(mols)
        
    return allowed_structure_site_structures

def get_best_adsorbate_geometries(adsorbate_path,aseinds,siteinfo,sites,site_adjacency,imag_freq_max=150.0):
    """
    load the adsorbates geometries find the best matching unique optimized
    adsorbate structures for each species
    siteinfo = [("ontop","terrace"),("edge","terrace")]
    aseinds = [37,38]
    """
    info_path = os.path.join(adsorbate_path,"info.json")
    with open(info_path,"r") as f:
        info = json.load(f)
    nslab = info["nslab"]
    atom_to_molecule_surface_atom_map = {int(k):int(v) for k,v in info["gratom_to_molecule_surface_atom_map"].items()}
    mol = Molecule().from_adjacency_list(info["adjlist"])
    prefixes = os.listdir(adsorbate_path)
    if len(prefixes) == 2: #includes info.json
        return read(os.path.join(adsorbate_path,"0","0.xyz"))
    geoms = []
    bestxyz = None
    best = None
    best_fingerprint = None
    best_score = None
    best_energy = None
    target_fingerprint = [d for x in siteinfo for d in x]
    for prefix in prefixes:
        if prefix == "info.json":
            continue
        path = os.path.join(adsorbate_path,prefix,prefix+".xyz")
        freq_path = os.path.join(adsorbate_path,prefix,"vib.json_vib.json")
        if os.path.exists(path) and os.path.exists(freq_path):
            with open(freq_path,"r") as f:
                freq_dict = json.load(f)
            freq = np.complex128(freq_dict["frequencies"][0])
            if (np.isreal(freq) or freq.imag < imag_freq_max):
                geoms.append(path)
                
    xyzs = geoms
    
    adsorbates = []
    for xyz in xyzs:
        geo = read(xyz)
        occ = get_occupied_sites(geo,sites,nslab)
        required_surface_inds = set([ind+nslab for ind in atom_to_molecule_surface_atom_map.keys()])
        found_surface_inds = set([site["bonding_index"] for site in occ])
        if len(occ) >= len(mol.get_adatoms()) and required_surface_inds.issubset(found_surface_inds):
            if best is None:
                best = geo
                bestxyz = xyz
                best_energy = geo.get_potential_energy()
                assert len(siteinfo) <= len(occ)
                for slist in itertools.permutations(occ,r=len(siteinfo)):
                    sinfo = [(site["site"],site["morphology"]) for site in slist]
                    fingerprint = [d for x in sinfo for d in x]
                    if best_fingerprint is None:
                        best_fingerprint = fingerprint
                        best_score = np.equal(np.array(best_fingerprint),np.array(target_fingerprint)).sum()
                    else:
                        score = np.equal(np.array(fingerprint),np.array(target_fingerprint)).sum()
                        if score > best_score:
                            best_fingerprint = fingerprint
                            best_score = score
                                  
            else:
                for slist in itertools.permutations(occ,r=len(siteinfo)):
                    sinfo = [(site["site"],site["morphology"]) for site in slist]
                    fingerprint = [d for x in sinfo for d in x]
                    score = np.equal(np.array(fingerprint),np.array(target_fingerprint)).sum()
                    if score > best_score or score == best_score and geo.get_potential_energy() < best_energy:
                        best = geo
                        bestxyz = xyz
                        best_energy = geo.get_potential_energy()
                        best_fingerprint = fingerprint
                        best_score = score
    
    return best,bestxyz

def get_best_reaction_adsorbate_geometries(admol,admol_neighbors,nslab2D,adsorbates_path,info,sites,
                                           site_adjacency,imag_freq_max=150.0,
                                           assume_reactant_product_templates_same_order=False):

    assert assume_reactant_product_templates_same_order #otherwise tricky to get formed bonds
    
    #get reverse template map
    if info["forward"]:
        template = Molecule().from_adjacency_list(info["reactants"])
        rev_template = Molecule().from_adjacency_list(info["products"])
    else:
        template = Molecule().from_adjacency_list(info["products"])
        rev_template = Molecule().from_adjacency_list(info["reactants"])
    
    broken_bonds,formed_bonds = get_broken_formed_bonds(template,rev_template)
    
    rev_mols = []
    for n in info["reverse_names"]:
        molinfo_path = os.path.join(adsorbates_path,n,"info.json")
        with open(molinfo_path,"r") as f:
            molinfo = json.load(f)
        mol = Molecule().from_adjacency_list(molinfo["adjlist"])
        rev_mols.append(mol)
    
    
    rev_ads_sizes = [ads_size(m) for m in rev_mols]
    molecule_to_atom_rev_maps = []
    for m in rev_mols:
        ads,mol_to_atoms_map = get_adsorbate(mol)
        molecule_to_atom_rev_maps.append(mol_to_atoms_map)
    molecule_to_atom_rev_maps_invert = [{int(v):int(k) for k,v in x.items()} for x in molecule_to_atom_rev_maps]
    
    template_rev_mol_map = get_template_mol_map(rev_template,rev_mols)
    
    template_rev_mol_map = get_template_mol_map(rev_template,rev_mols)
    template_rev_mol_map_invert = [{int(v):int(k) for k,v in x.items()} for x in template_rev_mol_map]
    template_mol_map = [{int(k):int(v) for k,v in x.items()} for x in info["template_mol_map"]]
    template_mol_map_invert = [{int(v):int(k) for k,v in x.items()} for x in info["template_mol_map"]]
    molecule_to_atom_maps = [{int(k):int(v) for k,v in x.items()} for x in info["molecule_to_atom_maps"]]
    molecule_to_atom_maps_invert = [{int(v):int(k) for k,v in x.items()} for x in info["molecule_to_atom_maps"]]
    ads_sizes = info["ads_sizes"]
    nslab = info["nslab"]
    mols = [Molecule().from_adjacency_list(x) for x in info["mols"]]
    
    #map adsorbed atoms on admol to templates
    adatoms = [a for a in admol.get_adatoms() if not a.is_surface_site()]
    adinds = [admol.atoms.index(adatom) - nslab2D for adatom in adatoms] #aseinds
    template_adinds = [ase_to_template_index(a,template_mol_map_invert,molecule_to_atom_maps_invert,
                                             ads_sizes) for a in adinds]
    forward_geoms = []
    for i,ad in enumerate(info["species_names"]):
        tempinds = [k for k in template_adinds if k in template_mol_map[i].keys()] 
        aseinds = [get_ase_index(k,template_mol_map,molecule_to_atom_maps,info["nslab"],ads_sizes) for k in tempinds]
        siteinfo = []
        for j,a in enumerate(adatoms):
            for bd in admol.get_bonds(a).values():
                if bd.atom1.is_surface_site():
                    ind = admol.atoms.index(bd.atom1)
                    indad = admol.atoms.index(bd.atom2)
                    break
                elif bd.atom2.is_surface_site():
                    ind = admol.atoms.index(bd.atom2)
                    indad = admol.atoms.index(bd.atom1)
                    break
            indtmp = ase_to_template_index(indad-nslab2D,template_mol_map_invert, 
                                           molecule_to_atom_maps_invert, ads_sizes)
            molind,mind = get_mol_index(indtmp,template_mol_map)
            
            if molind == i and mols[i].atoms[mind].is_bonded_to_surface():
                site = admol_neighbors[ind]
                siteinfo.append((site["site"],site["morphology"]))
        
        best = get_best_adsorbate_geometries(os.path.join(adsorbates_path,ad),aseinds,siteinfo,sites,site_adjacency,
                                             imag_freq_max=imag_freq_max)[0]
        forward_geoms.append(best)
        
    reverse_geoms = []  
    for j,ad in enumerate(info["reverse_names"]):
        tempinds = [k for k in template_adinds if k in template_rev_mol_map[j].keys()] 
        aseinds = [get_ase_index(k,template_mol_map,molecule_to_atom_maps,info["nslab"],ads_sizes) for k in tempinds]
        siteinfo = []
        for i,a in enumerate(adatoms):
            for bd in admol.get_bonds(a).values():
                if bd.atom1.is_surface_site():
                    ind = admol.atoms.index(bd.atom1)
                    indad = admol.atoms.index(bd.atom2)
                    break
                elif bd.atom2.is_surface_site():
                    ind = admol.atoms.index(bd.atom2)
                    indad = admol.atoms.index(bd.atom1)
                    break
            #assuming rev and forward template indices are the same
            indtmp = ase_to_template_index(indad-nslab2D,template_mol_map_invert, 
                                           molecule_to_atom_maps_invert, ads_sizes)

            molind,mind = get_mol_index(indtmp,template_rev_mol_map)
            if molind == j and rev_mols[j].atoms[mind].is_bonded_to_surface():
                site = admol_neighbors[ind]
                siteinfo.append((site["site"],site["morphology"]))

        best = get_best_adsorbate_geometries(os.path.join(adsorbates_path,ad),aseinds,siteinfo,sites,site_adjacency,
                                             imag_freq_max=imag_freq_max)[0]
        reverse_geoms.append(best)
    
    return forward_geoms,reverse_geoms
