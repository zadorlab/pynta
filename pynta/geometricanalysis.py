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

def fix_bond_orders(mol,allow_failure=False,keep_binding_vdW_bonds=False,keep_vdW_surface_bonds=False):
    atom_order = ["C","Si","N","P","O","S","F","Cl","Br","I","Li","H"]#["H","Li","I","Br","Cl","F","S","O","P","N","Si","C"]#["C","Si","N","P","O","S","F","Cl","Br","I","Li","H"]
    finished = False
    cleanup_surface_bonds = not (keep_vdW_surface_bonds or keep_binding_vdW_bonds)
    atoms_to_fix = [a for a in mol.atoms if not a.is_surface_site()]
    while len(atoms_to_fix) > 0:
        atoms_to_remove = []
        atoms_modified = False
        
        for a in atoms_to_fix:
            boo = fix_atom_only_required(mol,a,allow_failure=allow_failure,cleanup_surface_bonds=cleanup_surface_bonds)
            if boo:
                atoms_modified = True
                atoms_to_remove.append(a)
                
        for a in atoms_to_remove:
            atoms_to_fix.remove(a)
        
        if not atoms_modified: 
            atoms_modified = False
            a_to_remove = None
            for target in atom_order:
                for a in atoms_to_fix:
                    if a.symbol == target:
                        boo = fix_atom(mol,a,allow_failure=allow_failure,cleanup_surface_bonds=cleanup_surface_bonds)
                        if boo:
                            atoms_modified = True
                            a_to_remove = a
                            break
                if atoms_modified:
                    break
            
            if a_to_remove is not None:
                atoms_to_fix.remove(a_to_remove)     
    
    if keep_vdW_surface_bonds:
        vdW_bonds = []
        for bd in mol.get_all_edges():
            if bd.order == 0:
                if not bd.atom1.is_surface_site() and not bd.atom2.is_surface_site():
                    vdW_bonds.append(bd)
        for bd in vdW_bonds:
            mol.remove_bond(bd)
            
    elif keep_binding_vdW_bonds:
        mtest = mol.copy(deep=True)
        
        vdW_bonds = []
        for bd in mtest.get_all_edges():
            if bd.order == 0:
                vdW_bonds.append(bd)
        for bd in vdW_bonds:
            mtest.remove_bond(bd)
            
        split_structs,atom_mapping = split_adsorbed_structures(mtest,clear_site_info=True,atom_mapping=True)
        
        gas_cluster_inds = []
        for sstruct in split_structs:
            if not any(a.is_surface_site() for a in sstruct.atoms): #nominal gas phase construct
                gas_cluster_inds.extend([atom_mapping[a] for a in sstruct.atoms])
        
        bd_to_remove = []
        for bd in mol.get_all_edges():
            if bd.order == 0:
                if (bd.atom1.is_surface_site() and mol.atoms.index(bd.atom2) in gas_cluster_inds) or (bd.atom2.is_surface_site() and mol.atoms.index(bd.atom1) in gas_cluster_inds):
                    pass
                else:
                    bd_to_remove.append(bd)
        
        for bd in bd_to_remove:
            mol.remove_bond(bd)
    
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

class FailedFixBondsException(Exception):
    pass 

def fix_atom_only_required(mol,atom,allow_failure=False,cleanup_surface_bonds=True):
    """
    Only increment bonds for which there is no other way to achieve octet
    """
    fixed = False
    if atom.is_bonded_to_surface(): #can always pull electrons from surface
        return False
    
    delta = get_octet_deviation(atom)
    
    if delta < 0: #too many electrons
        raise TooManyElectronsException("Cannot solve problem of too many electrons not counting surface bonds and all ad bonds are 1: {}".format(mol.to_adjacency_list()))
    elif delta > 0: #too few electrons
        atoms = [k for k,v in atom.bonds.items()] #none of these will be sites
        octets = [get_octet_deviation(k) for k in atoms]
        if len([oct for oct in octets if oct > 0]) != 1 and sum(octets) > delta: #if there is more than enough electrons and multiple ways to distribute them exit
            return False
        
        symbols = [k.symbol for k in atoms]
        
        orders = [v.get_order_str() for k,v in atom.bonds.items()]
        bonded_to_surface = [a.is_bonded_to_surface() for a in atoms]

        while delta > 0:
            ind,val = choose_bond_to_fix(symbols,octets,orders,bonded_to_surface)
            fixed = True
            if val == 0:
                if allow_failure:
                    break
                else:
                    logging.error(mol.to_adjacency_list())
                    logging.error(mol.atoms.index(atom))
                    raise FailedFixBondsException
            bd = atom.bonds[atoms[ind]]
            bd.increment_order()
            delta -= 2
            octets[ind] -= 2
    
    else:
        fixed = True #we're happy with this 
        
    #at end clean up improper surface bonds
    if cleanup_surface_bonds:
        to_remove = []
        for k,v in atom.bonds.items():
            if k.symbol == 'X' and v.order == 0:    
                to_remove.append(v)
        
        for v in to_remove:
            mol.remove_bond(v)
    
    return fixed
            
def fix_atom(mol,atom,allow_failure=False,cleanup_surface_bonds=True):
    
    fixed = False
    delta = get_octet_deviation(atom)
    
    if delta < 0: #too many electrons
        raise TooManyElectronsException("Cannot solve problem of too many electrons not counting surface bonds and all ad bonds are 1: {}".format(mol.to_adjacency_list()))
    elif delta > 0: #too few electrons
        atoms = [k for k,v in atom.bonds.items()]
        symbols = [k.symbol for k in atoms]
        octets = [get_octet_deviation(k) for k in atoms]
        orders = [v.get_order_str() for k,v in atom.bonds.items()]
        bonded_to_surface = [a.is_bonded_to_surface() for a in atoms]
        while delta > 0:
            ind,val = choose_bond_to_fix(symbols,octets,orders,bonded_to_surface)
            fixed = True
            if val == 0:
                if allow_failure:
                    break
                else:
                    raise FailedFixBondsException
            bd = atom.bonds[atoms[ind]]
            bd.increment_order()
            delta -= 2
            octets[ind] -= 2
    else:
        fixed = True
        
    #at end clean up improper surface bonds
    if cleanup_surface_bonds:
        to_remove = []
        for k,v in atom.bonds.items():
            if k.symbol == 'X' and v.order == 0:    
                to_remove.append(v)
        
        for v in to_remove:
            mol.remove_bond(v)
    
    return fixed

def generate_adsorbate_2D(atoms, sites, site_adjacency, nslab, max_dist=3.0, cut_off_num=None, allowed_structure_site_structures=None,
                          keep_binding_vdW_bonds=False, keep_vdW_surface_bonds=False):
    admol,neighbor_sites,ninds = generate_adsorbate_molecule(atoms, sites, site_adjacency, nslab, max_dist=max_dist)

    if cut_off_num:
        bds_to_remove = []
        for i in range(len(admol.atoms)-cut_off_num,len(admol.atoms)):         
            for a,b in admol.atoms[i].edges.items():
                if not a.is_surface_site() and (a.is_bonded_to_surface() or admol.atoms.index(a) < len(admol.atoms)-cut_off_num):
                    if b not in bds_to_remove:
                        bds_to_remove.append(b)
        for b in bds_to_remove:
            admol.remove_bond(b)
    
    fix_bond_orders(admol,keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
    
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

        if cut_off_num:
            bds_to_remove = []
            for i in range(len(admol.atoms)-cut_off_num,len(admol.atoms)):        
                for a,b in admol.atoms[i].edges.items():
                    if not a.is_surface_site() and (a.is_bonded_to_surface() or admol.atoms.index(a) < len(admol.atoms)-cut_off_num):
                        if b not in bds_to_remove:
                            bds_to_remove.append(b)
            for b in bds_to_remove:
                admol.remove_bond(b)
            
        fix_bond_orders(admol,keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
    
    admol.update_multiplicity()
    admol.update_atomtypes()
    admol.update_connectivity_values()
    admol.identify_ring_membership()
    
    return admol,neighbor_sites,ninds

def generate_TS_2D(atoms, info_path,  metal, facet, sites, site_adjacency, nslab, imag_freq_path=None,
                     max_dist=3.0, cut_off_num=None, allowed_structure_site_structures=None,
                     site_bond_cutoff=3.0,keep_binding_vdW_bonds=False,keep_vdW_surface_bonds=False,reaction_adsorbate_bond_cutoff=2.3):

    admol,neighbor_sites,ninds = generate_adsorbate_molecule(atoms, sites, site_adjacency, 
                                                             nslab, max_dist=max_dist)

    if cut_off_num:
        bds_to_remove = []
        for i in range(len(admol.atoms)-cut_off_num,len(admol.atoms)):
            for a,b in admol.atoms[i].edges.items():
                if not a.is_surface_site() and (a.is_bonded_to_surface() or admol.atoms.index(a) < len(admol.atoms)-cut_off_num):
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
        inds3D = []
        for label in bd:
            a = template.get_labeled_atoms(label)[0]
            aind = template.atoms.index(a)
            aseind = get_ase_index(aind,template_mol_map,molecule_to_atom_maps,nslab,info["ads_sizes"])
            if aseind:
                inds3D.append(aseind)
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
                    raise SiteOccupationException("couldn't find right site bond")
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
            if len(inds3D) == 2:
                v,dist = get_distances(atoms.positions[inds3D[0]],atoms.positions[inds3D[1]],cell=atoms.cell,pbc=atoms.pbc)
                if dist < reaction_adsorbate_bond_cutoff:
                    admol.add_bond(Bond(admol.atoms[inds2D[0]],admol.atoms[inds2D[1]],"R"))
            else:
                admol.add_bond(Bond(admol.atoms[inds2D[0]],admol.atoms[inds2D[1]],"R"))
            
    fix_bond_orders(admol,allow_failure=True,keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)

    
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

    if cut_off_num:
        bds_to_remove = []
        for i in range(len(admol.atoms)-cut_off_num,len(admol.atoms)):  
            for a,b in admol.atoms[i].edges.items():
                if not a.is_surface_site() and (a.is_bonded_to_surface() or admol.atoms.index(a) < len(admol.atoms)-cut_off_num):
                    if b not in bds_to_remove:
                        bds_to_remove.append(b)
        for b in bds_to_remove:
            admol.remove_bond(b)

    for bd in rbonds:
        inds2D = []
        inds3D = []
        for label in bd:
            a = template.get_labeled_atoms(label)[0]
            aind = template.atoms.index(a)
            aseind = get_ase_index(aind,template_mol_map,molecule_to_atom_maps,nslab,info["ads_sizes"])
            if aseind:
                inds2D.append(aseind - nslab + len(neighbor_sites))
                inds3D.append(aseind)

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
                        if not np.isnan(motion_alignment) and motion_alignment > 0.95 and dist < 3.0 and dist_site < 2.0 and (not allowed_structure_site_structures or (site1_id in valid_sites and site2_id in valid_sites)):
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
            if len(inds3D) == 2:
                v,dist = get_distances(atoms.positions[inds3D[0]],atoms.positions[inds3D[1]],cell=atoms.cell,pbc=atoms.pbc)
                if dist < reaction_adsorbate_bond_cutoff:
                    admol.add_bond(Bond(admol.atoms[inds2D[0]],admol.atoms[inds2D[1]],"R"))
            else:
                admol.add_bond(Bond(admol.atoms[inds2D[0]],admol.atoms[inds2D[1]],"R"))

    fix_bond_orders(admol,allow_failure=True,keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds) #Reaction bonds prevent proper octet completion
    
    #remove surface bonds that abused for octet completion in reactions
    for atom in admol.atoms:
        surfbds = [bd for at,bd in atom.bonds.items() if at.is_surface_site() and bd.get_order_str() != 'R']
        if len(surfbds) == 1:
            Nrxnbd = len([bd for at,bd in atom.bonds.items() if bd.get_order_str() == 'R' and not at.is_surface_site()])
            if Nrxnbd >= 2 and abs(Nrxnbd/2 - surfbds[0].order) < 1e-3:
                admol.remove_edge(surfbds[0])
        
    admol.update_multiplicity()
    admol.update_atomtypes()
    admol.update_connectivity_values()
    admol.identify_ring_membership()
    
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

def split_ts_to_reactants(ts2d,tagatoms=False,keep_binding_vdW_bonds=False,keep_vdW_surface_bonds=False):
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
            fix_bond_orders(r,keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
            r.update(sort_atoms=False)
            r.update_connectivity_values()
            r.identify_ring_membership()
            rs.append(r)
        except (TooManyElectronsException,ValueError,FailedFixBondsException):
            pass
    
    return rs
      
def generate_allowed_structure_site_structures(adsorbate_dir,sites,site_adjacency,nslab,max_dist=3.0,cut_off_num=None):
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
            admol,neighbor_sites,ninds = generate_adsorbate_2D(geom, sites, site_adjacency, nslab, max_dist=max_dist, cut_off_num=cut_off_num)
            m = remove_slab(admol)
            if not generate_without_site_info(m).is_isomorphic(target_mol,save_order=True): #geom didn't optimize to target
                continue
            for mol in mols:
                if mol.is_isomorphic(m,save_order=True):
                    break
            else:
                m.update_atomtypes()
                m.update_connectivity_values()
                m.identify_ring_membership()
                mols.append(m)
        if mols:
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
        return read(os.path.join(adsorbate_path,"0","0.xyz")),os.path.join(adsorbate_path,"0","0.xyz")
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
    adatoms = np.unique(np.array([a for a in admol.get_adatoms() if not a.is_surface_site()])).tolist()
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

def ase_to_template_index(aseind, template_mol_map_invert, molecule_to_atom_maps_invert, ads_sizes, nslab=None, ):
    if nslab:
        aseind = aseind - nslab
        
    #go to molind
    n = 0
    for i,d in enumerate(molecule_to_atom_maps_invert):
        if aseind > n + len(d) - 1:
            n += ads_sizes[i]
            continue
        else:
            molnum = i
            molind = molecule_to_atom_maps_invert[molnum][aseind-n]
            break
    else:
        raise ValueError
    
    return template_mol_map_invert[molnum][molind]

def get_bond_factors(ts_path,adsorbates_path,metal,facet,sites,site_adjacency,max_dist=3.0,imag_freq_max=150.0):
    
    info_path = os.path.join(os.path.split(os.path.split(ts_path)[0])[0],"info.json")
    with open(info_path,"r") as f:
        info = json.load(f)
    
    template_mol_map = [{int(k):v for k,v in x.items()}  for x in info["template_mol_map"]]
    molecule_to_atom_maps = [{int(k):v for k,v in x.items()}  for x in info["molecule_to_atom_maps"]]
    broken_bonds,formed_bonds = get_broken_formed_bonds(Molecule().from_adjacency_list(info["reactants"]),
                                                        Molecule().from_adjacency_list(info["products"]))
    rbonds = broken_bonds | formed_bonds
    nslab = info["nslab"]
    
    atms = read(ts_path)
    
    admol,neighbor_sites,ninds = generate_TS_2D(atms, info_path, metal, facet, sites, site_adjacency, 
                                                nslab, max_dist=max_dist)
        
    occ = get_occupied_sites(atms,sites,nslab)
    occ_bd_lengths = {site["bonding_index"]:site['bond_length'] for site in occ}
    
    reactants,products = get_best_reaction_adsorbate_geometries(admol,neighbor_sites,len([a for a in admol.atoms if a.is_surface_site()]),
                                                adsorbates_path,info,sites,site_adjacency,
                                                                imag_freq_max=imag_freq_max,
                                               assume_reactant_product_templates_same_order=True)

    if info["forward"]:
        template = Molecule().from_adjacency_list(info["reactants"])
        rev_template = Molecule().from_adjacency_list(info["products"])
    else:
        template = Molecule().from_adjacency_list(info["products"])
        rev_template = Molecule().from_adjacency_list(info["reactants"])

    rev_mols = []
    molecule_to_atom_rev_maps = []
    for n in info["reverse_names"]:
        molinfo_path = os.path.join(adsorbates_path,n,"info.json")
        with open(molinfo_path,"r") as f:
            molinfo = json.load(f)
        mol = Molecule().from_adjacency_list(molinfo["adjlist"])
        rev_mols.append(mol)
        molecule_to_atom_rev_maps.append({int(v):int(k) for k,v in molinfo["atom_to_molecule_atom_map"].items()})

    template_rev_mol_map = get_template_mol_map(rev_template,rev_mols)
    rev_ads_sizes = [ads_size(mol) for mol in rev_mols]

    template_mol_map_invert = [{int(v):int(k) for k,v in x.items()} for x in info["template_mol_map"]]
    molecule_to_atom_maps_invert = [{int(v):int(k) for k,v in x.items()} for x in info["molecule_to_atom_maps"]]
    ads_sizes = info["ads_sizes"]
    mols = [Molecule().from_adjacency_list(x) for x in info["mols"]]
    
    adjlist_to_length = dict()
    adjlist_to_well_length = dict()
    
    for bd in rbonds:
        aseinds = []
        inds2D = []
            
        for label in bd:
            a = template.get_labeled_atoms(label)[0]
            aind = template.atoms.index(a)
            aseind = get_ase_index(aind,template_mol_map,molecule_to_atom_maps,nslab,info["ads_sizes"])
            if aseind:
                aseinds.append(aseind)
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
                    raise SiteOccupationException(admol.to_adjacency_list())
                for j,s in enumerate(neighbor_sites):
                    if (s["position"] == pos).all():
                        inds2D.append(j)
                        break
        
        m = admol.copy(deep=True)
        for i in inds2D:
            m.atoms[i].label = "*"
        
        templateinds = [ase_to_template_index(aind, template_mol_map_invert, molecule_to_atom_maps_invert, ads_sizes, nslab) for aind in aseinds]
        if len(templateinds) == 2:
            forward = template.has_bond(template.atoms[templateinds[0]],template.atoms[templateinds[1]])
        else:
            forward = template.atoms[templateinds[0]].is_bonded_to_surface()

        if len(aseinds) == 2:
            L = atms.get_distance(aseinds[0],aseinds[1],mic=True)
            if forward:
                molind = get_mol_index(templateinds[0],template_mol_map)[0]
                aind1 = get_ase_index(templateinds[0],template_mol_map,molecule_to_atom_maps,nslab,ads_sizes)
                aind2 = get_ase_index(templateinds[1],template_mol_map,molecule_to_atom_maps,nslab,ads_sizes)
                if mols[molind].contains_surface_site():
                    Lwell = reactants[molind].get_distance(aind1,aind2,mic=True)
                else:
                    Lwell = reactants[molind].get_distance(aind1-nslab,aind2-nslab,mic=True)
            else:
                molind,mind = get_mol_index(templateinds[0],template_rev_mol_map)
                aind1 = get_ase_index(templateinds[0],template_rev_mol_map,molecule_to_atom_rev_maps,nslab,rev_ads_sizes)
                aind2 = get_ase_index(templateinds[1],template_rev_mol_map,molecule_to_atom_rev_maps,nslab,rev_ads_sizes)                
                if rev_mols[molind].contains_surface_site():
                    Lwell = products[molind].get_distance(aind1,aind2,mic=True)
                else:
                    Lwell = products[molind].get_distance(aind1-nslab,aind2-nslab,mic=True)

        elif len(aseinds) == 1:
            L = occ_bd_lengths[aseinds[0]]
            if forward:
                molind,mind = get_mol_index(templateinds[0],template_mol_map)
                aind1 = molecule_to_atom_maps[molind][mind]+nslab
                mocc = get_occupied_sites(reactants[molind],sites,nslab)
                mocc_bd_lengths = {site["bonding_index"]:site['bond_length'] for site in mocc}
                Lwell = mocc_bd_lengths[aind1]
            else:
                molind,mind = get_mol_index(templateinds[0],template_rev_mol_map)
                aind1 = molecule_to_atom_maps[molind][mind]+nslab
                mocc = get_occupied_sites(products[molind],sites,nslab)
                mocc_bd_lengths = {site["bonding_index"]:site['bond_length'] for site in mocc}
                Lwell = mocc_bd_lengths[aind1]
                
        if Lwell > 2.0:
            assert False
        adjlist_to_length[m.to_adjacency_list()] = L
        adjlist_to_well_length[m.to_adjacency_list()] = Lwell
    
    adjlist_to_bond_factor = dict()
    for k in adjlist_to_length.keys():
        adjlist_to_bond_factor[k] = adjlist_to_length[k]/adjlist_to_well_length[k]
    
    return adjlist_to_bond_factor,adjlist_to_length,adjlist_to_well_length

def get_unique_adsorbate_geometries(adsorbate_path,mol,sites,site_adjacency,atom_to_molecule_surface_atom_map,
                                    nslab,imag_freq_max=150.0,return_xyzs=False):
    """
    load the adsorbates associated with the reaction and find the unique optimized
    adsorbate structures for each species
    returns a dictionary mapping each adsorbate name to a list of ase.Atoms objects
    """
    prefixes = os.listdir(adsorbate_path)
    geoms = []
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
    
    ase_atoms = [read(geom) for geom in geoms]
    mols = [generate_adsorbate_2D(geom, sites, site_adjacency, nslab, max_dist=np.inf)[0] for geom in ase_atoms]
    filtered_mols = []
    filtered_geoms = []
    filtered_xyzs = []
    filtered_energies = []
    for j,m in enumerate(mols):
        for i,fm in enumerate(filtered_mols):
            if fm.is_isomorphic(m,save_order=True):
                if filtered_energies[i] > ase_atoms[j].get_potential_energy():
                    filtered_mols[i] = m 
                    filtered_geoms[i] = ase_atoms[j]
                    filtered_energies[i] = ase_atoms[j].get_potential_energy()
                    filtered_xyzs[i] = geoms[j]
        else:
            filtered_mols.append(m)
            filtered_geoms.append(ase_atoms[j])
            filtered_energies.append(ase_atoms[j].get_potential_energy())
            filtered_xyzs.append(geoms[j])
            
    adsorbates = []
    xyzs_out = []
    for j,geo in enumerate(filtered_geoms):
        occ = get_occupied_sites(geo,sites,nslab)
        required_surface_inds = set([ind+nslab for ind in atom_to_molecule_surface_atom_map.keys()])
        found_surface_inds = set([site["bonding_index"] for site in occ])
        if len(occ) >= len(mol.get_adatoms()) and required_surface_inds.issubset(found_surface_inds):
            adsorbates.append(geo)
            xyzs_out.append(filtered_xyzs[j])
    
    if not return_xyzs:
        return adsorbates
    else:
        return adsorbates,xyzs_out

def get_adsorbate_geometries(adsorbate_path,mol,sites,atom_to_molecule_surface_atom_map,
                                    nslab,imag_freq_max=150.0):
    """
    load the adsorbates associated with the reaction and find the unique optimized
    adsorbate structures for each species
    returns a dictionary mapping each adsorbate name to a list of ase.Atoms objects
    """
    prefixes = os.listdir(adsorbate_path)
    geoms = []
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
                
    xyzs = get_unique_sym(geoms)
    adsorbates = []
    for xyz in xyzs:
        geo = read(xyz)
        occ = get_occupied_sites(geo,sites,nslab)
        required_surface_inds = set([ind+nslab for ind in atom_to_molecule_surface_atom_map.keys()])
        found_surface_inds = set([site["bonding_index"] for site in occ])
        if len(occ) >= len(mol.get_adatoms()) and required_surface_inds.issubset(found_surface_inds):
            adsorbates.append(geo)

    return adsorbates

def get_lowest_adsorbate_energies(adsorbates_path):
    """
    load the adsorbates geometries find the lowest energy correct structure and the 2D representation
    """
    ad_energy_dict = dict()
    for ad in os.listdir(adsorbates_path):
        path = os.path.join(adsorbates_path,ad)
        if not os.path.isdir(path):
            continue
        with open(os.path.join(path,"info.json"),"r") as f:
            info = json.load(f)
            mol = Molecule().from_adjacency_list(info["adjlist"])
        Es = get_adsorbate_energies(path)[0]
        ad_energy_dict[mol] = np.inf
        for E in Es.values():
            if E < ad_energy_dict[mol]:
                ad_energy_dict[mol] = E
        
    return ad_energy_dict

def extract_pair_graph(atoms,sites,site_adjacency,nslab,max_dist,cut_off_num=None,allowed_structure_site_structures=None,
                       is_ts=False,info_path=None,metal=None,facet=None,imag_freq_path=None):
    if is_ts:
        admol,neighbor_sites,ninds = generate_TS_2D(atoms, info_path, metal, facet, sites, site_adjacency, nslab, imag_freq_path=imag_freq_path,
                        max_dist=max_dist, cut_off_num=cut_off_num, allowed_structure_site_structures=allowed_structure_site_structures,
                        )
    else:
        admol,neighbor_sites,ninds = generate_adsorbate_2D(atoms, sites, site_adjacency, nslab, 
                                                       max_dist=max_dist,cut_off_num=cut_off_num,allowed_structure_site_structures=allowed_structure_site_structures)
    return reduce_graph_to_pairs(admol)

def get_adsorbate_energies(ad_path,atom_corrections=None,include_zpe=True):
    """
    get ZPE corrected adsorbate energies
    """
    dirs = os.listdir(ad_path)
    slab = read(os.path.join(os.path.split(os.path.split(ad_path)[0])[0],"slab.xyz"))
    Eslab = slab.get_potential_energy()
    with open(os.path.join(ad_path,"info.json")) as f:
        info = json.load(f)

    m = Molecule().from_adjacency_list(info["adjlist"])
    if atom_corrections:
        AEC = 0.0
        for atom in m.atoms:
            s = str(atom.element)
            if not "X" in s:
                AEC += atom_corrections[s]
    else: #if no corrections don't add a correction
        AEC = 0.0
    
    gasphase = len(info["gratom_to_molecule_surface_atom_map"]) == 0
    spin = (Molecule().from_adjacency_list(info["adjlist"]).multiplicity - 1.0)/2.0
    Es = dict()
    thermos = dict()
    fs = dict()
    for d in dirs:
        if d == "info.json":
            continue
        optdir = os.path.join(ad_path,d,d+".xyz")
        freqdir = os.path.join(ad_path,d,"vib.json_vib.json")
        if not (os.path.exists(optdir) and os.path.exists(freqdir)):
            continue
        sp = read(optdir)
        E = sp.get_potential_energy()
        if gasphase:
            vibdata = get_vibdata(optdir,freqdir,0)
        else:
            vibdata = get_vibdata(optdir,freqdir,len(slab))

        freqs = vibdata.get_frequencies().tolist()
        fs[d] = freqs

        ZPE = vibdata.get_zero_point_energy()

        if not gasphase:

            if include_zpe:
                Es[d] = E + ZPE - Eslab - AEC
            else:
                Es[d] = E - Eslab - AEC
            thermos[d] = HarmonicThermo(np.array([x for x in np.real(vibdata.get_energies()) if x > 0.0]),potentialenergy=E-Eslab-AEC)
        else:
            if len(sp) == 1:
                geometry = 'monatomic'
            elif any([ I < 1.0e-3 for I in sp.get_moments_of_inertia()]):
                geometry = "linear"
            else:
                geometry = "nonlinear"

            if include_zpe:
                Es[d] = E + ZPE - AEC
            else:
                Es[d] = E - AEC

            thermos[d] = IdealGasThermo(np.real(vibdata.get_energies()),geometry,
                                        potentialenergy=E - AEC,atoms=sp,symmetrynumber=1,
                                        natoms=len(sp),spin=spin)

    return Es,thermos,fs

def get_vibdata(dopt,dvib,nslab,gas_phase=False):
    xyz = read(dopt)
    if not gas_phase:
        inds = range(nslab,len(xyz))
    else:
        inds = range(0,len(xyz))
    with open(dvib,'r') as f:
        out = json.load(f)
    sh = (len(inds),3,len(inds),3)
    H = np.asarray(out["hessian"])
    if H.shape != sh:
        H = H.reshape(sh)
    vib = VibrationsData(xyz,H,inds)
    return vib

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
    if list(coad2D.get_surface_sites()[0].bonds.values())[0].is_van_der_waals():
        c = coad2D.copy(deep=True)
        c.clear_labeled_atoms()
        bd = list(c.get_surface_sites()[0].bonds.values())[0]
        if bd.atom1.is_surface_site():
            bd.atom2.label = "*0"
            c.remove_atom(bd.atom1)
        else:
            bd.atom1.label = "*0"
            c.remove_atom(bd.atom2)
    else:
        c = coad2D.get_desorbed_molecules()[0]
    mol2D = mol2D.merge(c)
    ldict = mol2D.get_all_labeled_atoms()
    label = list(ldict.keys())[0]
    catom = list(ldict.values())[0]
    catom.label = ''
    if label == "*0":
        bd = Bond(siteatom,catom,order='vdW')
    elif label == "*1":
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
    mol2D.identify_ring_membership()
    return mol2D

def mol_to_atoms(admol,slab,sites,metal,partial_atoms=None,partial_admol=None):
    """Generate a 3D initial guess for a given 2D admol configuration

    Args:
        admol (_type_): 2D representation of the configuration
        slab (_type_): 3D Atoms object for the slab
        sites (_type_): list of sites
        metal (_type_): metal of the slab
        partial_atoms (_type_, optional): a 3D Atoms object representing part of the target configuraitons. Defaults to None.
        partial_admol (_type_, optional): 2D representation of the partial_atoms Atoms object. Defaults to None.

    Raises:
        ValueError: _description_
        ValueError: _description_

    Returns:
        _type_: Atoms object corresponding to the admol 2D configuration
    """
    if partial_atoms and partial_admol:
        atoms = deepcopy(partial_atoms)
        gpartial_admol = partial_admol.to_group()
        admol_atom_order = admol.atoms[:]
        gpartial_admol_order = gpartial_admol.atoms[:]
        initial_map = dict()
        for i,a in enumerate(admol.atoms):
            if a.is_surface_site():
                assert a.site == gpartial_admol.atoms[i].site[0]
                initial_map[a] = gpartial_admol.atoms[i]
        try:
            subisos = admol.find_subgraph_isomorphisms(gpartial_admol,initial_map=initial_map,save_order=True)
        except ValueError:
            subisos = admol.find_subgraph_isomorphisms(gpartial_admol,initial_map=initial_map)
            admol.atoms = admol_atom_order
            gpartial_admol.atoms = gpartial_admol_order
        if len(subisos) == 0:
            raise ValueError("partial_admol is not subgraph isomorphic to admol: {0}, {1}".format(gpartial_admol.to_adjacency_list(),admol.to_adjacency_list()))
        subiso = subisos[0]
        atoms_in_partial = [a for a in subiso.keys() if not a.is_surface_site()]
        split_structs,adsorbed_atom_dict = split_adsorbed_structures(admol,clear_site_info=False,adsorption_info=True,atoms_to_skip=atoms_in_partial)
    elif partial_atoms is None and partial_admol is None:
        atoms = deepcopy(slab)
        split_structs,adsorbed_atom_dict = split_adsorbed_structures(admol,clear_site_info=False,adsorption_info=True)
    else:
        raise ValueError("Must include both partial_atoms, partial_admol to start from partial")
    
    for st in split_structs:
        if len(st.atoms) == 0:
            continue
        ad,mol_to_atoms_map = get_adsorbate(st)
        
        if len(st.atoms) > 2: #not atomic adsorbate
            rot_vec = np.array([0.0,0.0,0.0])
            for atst in st.atoms:
                if atst.is_bonded_to_surface():
                    adatom_molind = st.atoms.index(atst)
                    for bdat in atst.bonds.keys():
                        if not bdat.is_surface_site():
                            rot_vec += ad.positions[mol_to_atoms_map[st.atoms.index(bdat)]] - ad.positions[mol_to_atoms_map[adatom_molind]]
            rot_vec /= np.linalg.norm(rot_vec)
            ad.rotate(rot_vec,[0.0,0.0,1.0]) #rotate bonds toward the +z-axis
        else:
            adatom_molind = st.atoms.index([a for a in st.atoms if not a.is_surface_site()][0])

        atom_surf_inds = []
        ad_sites = []
        for a,sind in adsorbed_atom_dict.items():
            if a in st.atoms:
                assert sites[sind]["site"] == a.site, (admol.to_adjacency_list(),partial_admol.to_adjacency_list())
                atom_surf_inds.append(mol_to_atoms_map[adatom_molind])
                ad_sites.append(sites[sind])

        atoms,_,_ = place_adsorbate_covdep(ad,atoms,atom_surf_inds,ad_sites,metal)

    return atoms

def get_best_adsorbate_xyz(adsorbate_path,sites,site_adjacency,nslab,allowed_structure_site_structures,keep_binding_vdW_bonds,keep_vdW_surface_bonds):
    """
    load the adsorbates associated with the reaction and find the unique optimized
    adsorbate structures for each species
    returns a dictionary mapping each adsorbate name to a list of ase.Atoms objects
    """
    with open(os.path.join(adsorbate_path,"info.json"),"r") as f:
        info = json.load(f)
        
    mol = Molecule().from_adjacency_list(info["adjlist"])
    
    adsorbate = None
    min_energy = np.inf
    Es,_,freq = get_adsorbate_energies(adsorbate_path)
    for prefix in Es.keys():
        xyz = os.path.join(adsorbate_path,prefix,prefix+".xyz")
        strind = os.path.split(xyz)[1].split(".")[0]
        if strind not in Es.keys():
            continue
        geo = read(xyz)
        try:
            admol,_,_ = generate_adsorbate_2D(geo, sites, site_adjacency, nslab, max_dist=np.inf, allowed_structure_site_structures=allowed_structure_site_structures,
                                                keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
        except (SiteOccupationException,TooManyElectronsException,ValueError,FailedFixBondsException) as e:
            continue
        molp = split_adsorbed_structures(admol,clear_site_info=True)[0]
        if molp.is_isomorphic(mol,save_order=True) and Es[strind] < min_energy: #if matches target and is lower in energy
            adsorbate = xyz
            min_energy = Es[strind]

    return adsorbate

def validate_TS_configs(path,sites,site_adjacency,nslab,irc_concern_len=8):
    """analyze all TSs run for a given TS target in a Pynta run to determine validity

    Args:
        path (_type_): path to the overall TS target directory
        sites (_type_): list of all sites
        site_adjacency (_type_): site adjacency information
        nslab (_type_): number of atoms in the slab
        irc_concern_len (int, optional): length below which IRCs are considered short. Defaults to 8.

    Returns:
        dictionary mapping each TS guess index to overall validity assessment, dictionary mapping each TS guess index to the analysis information used to decide that
    """
    valid_dict = dict()
    valid_info = dict()
    for p in os.listdir(path):
        if not os.path.isdir(os.path.join(path,p)):
            continue
        if os.path.exists(os.path.join(path,p,"vib.json_vib.json")):
            valid,vinfo = validate_TS(os.path.join(path,p,"opt.xyz"),sites,site_adjacency,nslab,irc_concern_len=irc_concern_len)
            valid_dict[p] = valid
            valid_info[p] = vinfo
    return valid_dict,valid_info
            
def validate_TS(ts_path,sites,site_adjacency,nslab,irc_path1=None,irc_path2=None,info_path=None,adsorbate_dir=None,imag_freq_path=None,irc_concern_len=8):
    """Assess the validity of a transition state. Combines information from the TS geometry, imaginary frequency and if available the IRCs. 
    

    Args:
        ts_path (_type_): path to the TS optimized geometry file
        sites (_type_): list of all sites
        site_adjacency (_type_): site adjacency information
        nslab (_type_): number of atoms in the slab
        irc_path1 (_type_): path to the first IRC trajectory file
        irc_path2 (_type_): path to the second IRC trajectory file
        info_path (_type_, optional): path to the TS info.json file. Defaults to None.
        adsorbate_dir (_type_, optional): adsorbates directory corresponding to this TS. Defaults to None.
        imag_freq_path (_type_, optional): path to the imaginary frequency trajectory file. Defaults to None.
        irc_concern_len (int, optional): length below which IRCs are considered short. Defaults to 8.

    Returns:
        Overall assessment of TS validity (bool), detailed assessment criteria information (dict)
    """
    if adsorbate_dir is None:
        adsorbate_dir = os.path.join(os.path.split(os.path.split(os.path.split(ts_path)[0])[0])[0],"Adsorbates")
    
    if irc_path1 is None and os.path.exists(os.path.join(os.path.split(ts_path)[0],"irc_forward.traj")) and os.path.exists(os.path.join(os.path.split(ts_path)[0],"irc_reverse.traj")):
        irc_path1 = os.path.join(os.path.split(ts_path)[0],"irc_forward.traj")
        irc_path2 = os.path.join(os.path.split(ts_path)[0],"irc_reverse.traj")
    
    if info_path is None:
        info_path = os.path.join(os.path.split(os.path.split(ts_path)[0])[0],"info.json")
        
    if imag_freq_path is None:
        imag_freq_path = os.path.join(os.path.split(ts_path)[0],"vib.0.traj")
    
    with open(info_path,'r') as f:
        info = json.load(f)
    
    reactants = Molecule().from_adjacency_list(info["reactants"])
    products = Molecule().from_adjacency_list(info["products"])
    
    vdW_surface_bonds = False
    binding_vdW_bonds = False
    for mol in [reactants,products]:
        for bd in mol.get_all_edges():
            if bd.order == 0:
                if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
                    binding_vdW_bonds = True
                    m = mol.copy(deep=True)
                    b = m.get_bond(m.atoms[mol.atoms.index(bd.atom1)],m.atoms[mol.atoms.index(bd.atom2)])
                    m.remove_bond(b)
                    out = m.split()
                    if len(out) == 1: #vdW bond is not only thing connecting adsorbate to surface
                        vdW_surface_bonds = True
                        
    ts_val,struct_match,bond_alignment_validation,reaction_bond_alignments,fixed_bond_alignments = validate_TS_geometry(ts_path,reactants,products,sites,site_adjacency,nslab,info_path=info_path,
                         imag_freq_path=imag_freq_path,allowed_structure_site_structures=None,keep_binding_vdW_bonds=binding_vdW_bonds,keep_vdW_surface_bonds=vdW_surface_bonds)
    
    if irc_path1 is not None:
        both_valid,one_endpoint_valid,endpoints_match,target_endpoints_match,short_irc = validate_TS_ircs(irc_path1,irc_path2,reactants,products,sites,site_adjacency,nslab,
                                                                allowed_structure_site_structures=None,
                        binding_vdW_bonds=binding_vdW_bonds, vdW_surface_bonds=vdW_surface_bonds, irc_concern_len=irc_concern_len)
        
        d = {"TS_direct_Validation": ts_val, "IRC_Validation": both_valid, "One_Endpoint_Valid": one_endpoint_valid, "IRC_Endpoints_Match": endpoints_match,
                                                    "Short_IRC": short_irc, "Struct Validation": struct_match,
                                                    "Freq Alignment Validation": bond_alignment_validation, "Reaction_Bond_Freq_Alignments": reaction_bond_alignments,
                                                    "Fixed_Bond_Freq_Alignments": fixed_bond_alignments}
        if one_endpoint_valid is not None:
            return both_valid or (((one_endpoint_valid and endpoints_match) or target_endpoints_match) and ts_val), d
        else:
            return ts_val,d
    else:
        d = {"TS_direct_Validation": ts_val, "IRC_Validation": None, "One_Endpoint_Valid": None,"IRC_Endpoints_Match": None,
                                                    "Short_IRC": None, "Struct Validation": struct_match, "Freq Alignment Validation": bond_alignment_validation, 
                                                    "Reaction_Bond_Freq_Alignments": reaction_bond_alignments,
                                                    "Fixed_Bond_Freq_Alignments": fixed_bond_alignments}
        return ts_val,d

def validate_TS_ircs(irc_path1,irc_path2,reactants,products,sites,site_adjacency,nslab,allowed_structure_site_structures=None,
                     binding_vdW_bonds=False, vdW_surface_bonds=False, irc_concern_len=8, vdW_site_bond_cutoff=2.0):
    """Assess the validity of a transition state based on the intrinsic reaction coordinate calculation (IRC)

    Args:
        irc_path1 (_type_): path to the first IRC trajectory file
        irc_path2 (_type_): path to the second IRC trajectory file
        reactants (_type_): Molecule representation of reactants
        products (_type_): Molecule representation of products
        sites (_type_): list of all sites
        site_adjacency (_type_): site adjacency information
        nslab (_type_): number of atoms in the slab
        allowed_structure_site_structures (_type_, optional): allowed site binding information. Defaults to None.
        binding_vdW_bonds (bool, optional): if the target reactants/products have binding surface vdW bonds. Defaults to False.
        vdW_surface_bonds (bool, optional): if the target reactants/products have vdW surface vdW bonds that are not binding. Defaults to False.
        irc_concern_len (int, optional): length below which IRCs are considered short. Defaults to 8.
        vdW_site_bond_cutoff (float, optional): bond length above which surface vdW bonds to non-ontop sites are removed. Defaults to 2.0.

    Returns:
        overall assessment of TS validity (bool), if at least one IRC endpoint matches reactants/products (bool), if both IRC endpoints are the same (bool), if reactants and products are the same (bool), if the IRC was short (bool)
    """
    tr1 = Trajectory(irc_path1)
    try:
        admol1,neighbor_sites1,ninds1 = generate_adsorbate_2D(tr1[-1], sites, site_adjacency, nslab, max_dist=np.inf, cut_off_num=None, allowed_structure_site_structures=allowed_structure_site_structures,
                            keep_binding_vdW_bonds=True, keep_vdW_surface_bonds=vdW_surface_bonds)
    except (SiteOccupationException,FailedFixBondsException,TooManyElectronsException):
        return None,None,None,None,None
    
    vdw_ase_inds = []
    vdw_bds = []
    for bd in admol1.get_all_edges():
        if bd.is_van_der_waals():
            if bd.atom1.is_surface_site() and not bd.atom2.is_surface_site():
                vdw_ase_inds.append(admol1.atoms.index(bd.atom2) - len(neighbor_sites1) + nslab)
                vdw_bds.append(bd)
            elif not bd.atom1.is_surface_site() and bd.atom2.is_surface_site():
                vdw_ase_inds.append(admol1.atoms.index(bd.atom1) - len(neighbor_sites1) + nslab)
                vdw_bds.append(bd)
                
    occs = get_occupied_sites(tr1[-1],sites,nslab)
    
    vdw_bd_to_remove = []
    for i,ind in enumerate(vdw_ase_inds):
        for occ in occs:
            if occ["bonding_index"] == ind:
                if not binding_vdW_bonds and (occ["site"] == "ontop" or occ["bond_length"] > vdW_site_bond_cutoff):
                    vdw_bd_to_remove.append(vdw_bds[i])
    
    for bd in vdw_bd_to_remove:
        admol1.remove_bond(bd)
    
    split_structs1 = split_adsorbed_structures(admol1,clear_site_info=True,adsorption_info=False,atoms_to_skip=None,
                              atom_mapping=False,split_sites_with_multiple_adsorbates=True)
    
    tr2 = Trajectory(irc_path2)
    try:
        admol2,neighbor_sites2,ninds2 = generate_adsorbate_2D(tr2[-1], sites, site_adjacency, nslab, max_dist=np.inf, cut_off_num=None, allowed_structure_site_structures=allowed_structure_site_structures,
                            keep_binding_vdW_bonds=True, keep_vdW_surface_bonds=vdW_surface_bonds)
    except (SiteOccupationException,FailedFixBondsException,TooManyElectronsException):
        return False,None,None,None,None
    
    vdw_ase_inds = []
    vdw_bds = []
    for bd in admol2.get_all_edges():
        if bd.is_van_der_waals():
            if bd.atom1.is_surface_site() and not bd.atom2.is_surface_site():
                vdw_ase_inds.append(admol2.atoms.index(bd.atom2) - len(neighbor_sites2) + nslab)
                vdw_bds.append(bd)
            elif not bd.atom1.is_surface_site() and bd.atom2.is_surface_site():
                vdw_ase_inds.append(admol2.atoms.index(bd.atom1) - len(neighbor_sites2) + nslab)
                vdw_bds.append(bd)
                
    occs = get_occupied_sites(tr2[-1],sites,nslab)
    
    vdw_bd_to_remove = []
    for i,ind in enumerate(vdw_ase_inds):
        for occ in occs:
            if occ["bonding_index"] == ind:
                if not binding_vdW_bonds and (occ["site"] == "ontop" or occ["bond_length"] > vdW_site_bond_cutoff):
                    vdw_bd_to_remove.append(vdw_bds[i])
    
    for bd in vdw_bd_to_remove:
        admol2.remove_bond(bd)
        
    split_structs2 = split_adsorbed_structures(admol2,clear_site_info=True,adsorption_info=False,atoms_to_skip=None,
                              atom_mapping=False,split_sites_with_multiple_adsorbates=True)
        
    rs = reactants.split()
    ps = products.split()
    for r in rs+ps:
        r.update_multiplicity()
        r.update_atomtypes()
        r.update_connectivity_values()
        r.identify_ring_membership()
        
    reactant_1_mapping = dict()
    reactant_1_match = True
    for i,st in enumerate(split_structs1):
        for j,r in enumerate(rs):
            if j in reactant_1_mapping.keys():
                continue
            if st.is_isomorphic(r,save_order=True,strict=False):
                reactant_1_mapping[j] = i
                break 
        else:
            reactant_1_match = False
    
    reactant_2_mapping = dict()
    reactant_2_match = True
    for i,st in enumerate(split_structs2):
        for j,r in enumerate(rs):
            if j in reactant_2_mapping.keys():
                continue
            if st.is_isomorphic(r,save_order=True,strict=False):
                reactant_2_mapping[j] = i
                break 
        else:
            reactant_2_match = False
    
    product_1_mapping = dict()
    product_1_match = True
    for i,st in enumerate(split_structs1):
        for j,r in enumerate(ps):
            if j in product_1_mapping.keys():
                continue
            if st.is_isomorphic(r,save_order=True,strict=False):
                product_1_mapping[j] = i
                break 
        else:
            product_1_match = False
    
    product_2_mapping = dict()
    product_2_match = True
    for i,st in enumerate(split_structs2):
        for j,r in enumerate(ps):
            if j in product_2_mapping.keys():
                continue
            if st.is_isomorphic(r,save_order=True,strict=False):
                product_2_mapping[j] = i
                break 
        else:
            product_2_match = False
    
    #Check if IRC reactants and products are the same
    endpoints_mapping = dict() 
    endpoints_match = True
    for i,st in enumerate(split_structs1):
        for j,r in enumerate(split_structs2):
            if j in endpoints_mapping.keys():
                continue
            if st.is_isomorphic(r,save_order=True,strict=False):
                endpoints_mapping[j] = i
                break 
        else:
            endpoints_match = False
    
    #Check if target reactants and products are the same
    target_endpoints_mapping = dict() 
    target_endpoints_match = True
    for i,st in enumerate(rs):
        for j,r in enumerate(ps):
            if j in target_endpoints_mapping.keys():
                continue
            if st.is_isomorphic(r,save_order=True,strict=False):
                target_endpoints_mapping[j] = i
                break 
        else:
            target_endpoints_match = False
    
    if (reactant_1_match and product_2_match) or (reactant_2_match and product_1_match):
        both_valid = True 
    else:
        both_valid = False
    
    if reactant_1_match or product_2_match or reactant_2_match or product_1_match:
        one_endpoint_valid = True 
    else:
        one_endpoint_valid = False
    
    if len(tr1) < irc_concern_len or len(tr2) < irc_concern_len:
        short_irc = True
    else:
        short_irc = False
    
    return both_valid,one_endpoint_valid,endpoints_match,target_endpoints_match,short_irc

def validate_TS_geometry(opt_path,reactants,products,sites,site_adjacency,nslab,info_path=None,
                         imag_freq_path=None,allowed_structure_site_structures=None,keep_binding_vdW_bonds=False,keep_vdW_surface_bonds=False,
                         min_reaction_bond_alignment=0.7,max_fixed_bond_alignment=0.4):
    """Evaluates the validity of a transition state based on its geometry and frequencies
    1) the geometry is converted into graph form and compared with the target transition state graph
    2) the imaginary frequency projection onto covalent bonds is examined to see if covalent bonds are being correctly stretched/not stretched

    Args:
        opt_path (_type_): path to the optimized TS geometry
        reactants (_type_): Molecule representation of reactants
        products (_type_): Molecule representation of products
        sites (_type_): list of all sites
        site_adjacency (_type_): site adjacency information
        nslab (_type_): number of atoms in the slab
        info_path (_type_, optional): path to the TS info.json file. Defaults to None.
        imag_freq_path (_type_, optional): path to the imaginary frequency trajectory file. Defaults to None.
        allowed_structure_site_structures (_type_, optional): allowed site binding information. Defaults to None.
        keep_binding_vdW_bonds (bool, optional): during 3D to 2D keep vdW bonds that are needed to bind to the surface. Defaults to False.
        keep_vdW_surface_bonds (bool, optional): during 3D to 2D keep vdW bonds in general. Defaults to False.
        min_reaction_bond_alignment (float, optional): the minimum allowed bond alignment with imaginary frequency for a covalent bond breaking/forming in the reaction. Defaults to 0.7.
        max_fixed_bond_alignment (float, optional): the maximum allowed bond alignment with imaginary frequency for a covalent bond unchanged in the reaction. Defaults to 0.3.

    Returns:
        overall assessment of validity (bool), if the TS geometry structurally matches the target (bool), list of breaking/forming bond alignments (list), list of fixed bond alignments (list)
    """
    atoms = read(opt_path)
    
    broken_bonds,formed_bonds = get_broken_formed_bonds(reactants,products)
    target_TS = reactants.copy(deep=True)
    
    for labels in list(broken_bonds)+list(formed_bonds):
        label1,label2 = list(labels)
        a1 = target_TS.get_labeled_atoms(label1)[0]
        a2 = target_TS.get_labeled_atoms(label2)[0]
        if target_TS.has_bond(a1,a2):
            bd = target_TS.get_bond(a1,a2)
            bd.set_order_str('R')
        else:
            bd = Bond(a1,a2,order='R')
            target_TS.add_bond(bd)
    
    target_TS.clear_labeled_atoms()
    
    target_split = target_TS.split()
    
    for st in target_split:
        st.update_multiplicity()
        st.update_atomtypes()
        st.update_connectivity_values()
        st.identify_ring_membership()
    try:
        admol,neighbor_sites,ninds = generate_TS_2D(atoms, info_path,  None, None, sites, site_adjacency, nslab, imag_freq_path=imag_freq_path,
                        max_dist=np.inf, cut_off_num=None, allowed_structure_site_structures=allowed_structure_site_structures,
                        keep_binding_vdW_bonds=keep_binding_vdW_bonds,keep_vdW_surface_bonds=keep_vdW_surface_bonds)
        
        
        split_mols = split_adsorbed_structures(admol,clear_site_info=True,adsorption_info=False,atoms_to_skip=None,
                                atom_mapping=False,split_sites_with_multiple_adsorbates=True)

        for st in split_mols:
            st.update_multiplicity()
            st.update_atomtypes()
            st.update_connectivity_values()
            st.identify_ring_membership()
            
        mapping = dict()
        match = True
        for i,st in enumerate(split_mols):
            for j,r in enumerate(target_split):
                if j in mapping.keys():
                    continue
                if st.is_isomorphic(r,save_order=True,strict=False):
                    mapping[j] = i
                    break 
            else:
                match = False
    except (SiteOccupationException,TooManyElectronsException,ValueError):
        return False,False,False,[],[]
    
    #get imaginary frequency motion
    tr_ts = Trajectory(imag_freq_path)
    v_ts = np.array([get_distances(tr_ts[15].positions[i], tr_ts[16].positions[i], cell=tr_ts[15].cell, pbc=tr_ts[15].pbc)[0][0][0] for i in range(nslab,len(tr_ts[15].positions))])
    v_ts = v_ts.flatten()
    v_ts /= np.linalg.norm(v_ts)
    
    Nrbonds = len([bd for bd in admol.get_all_edges() if bd.is_reaction_bond()])
    
    reaction_bond_alignments = []
    fixed_bond_alignments = []
    for bd in admol.get_all_edges():
        if bd.atom1.is_surface_site() or bd.atom2.is_surface_site():
            continue
        else:
            admol_ind1 = admol.atoms.index(bd.atom1)
            admol_ind2 = admol.atoms.index(bd.atom2)
            ase_ind1 = admol_ind1 - len(neighbor_sites) + nslab
            ase_ind2 = admol_ind2 - len(neighbor_sites) + nslab
            bvec = get_distances(atoms.positions[ase_ind1],atoms.positions[ase_ind2],cell=atoms.cell,pbc=atoms.pbc)[0][0][0]
            v_b = []
            for i in range(nslab,len(atoms.positions)):
                if i == ase_ind1:
                    v_b.append(bvec)
                elif i == ase_ind2:
                    v_b.append(-bvec)
                else:
                    v_b.append(np.zeros(3))
            v_b = np.array(v_b).flatten()
            v_b /= np.linalg.norm(v_b)
            bond_alignment = Nrbonds * np.abs(np.dot(v_ts,v_b))
            if bd.is_reaction_bond():
                reaction_bond_alignments.append(bond_alignment)
            else:
                fixed_bond_alignments.append(bond_alignment)
    
    if fixed_bond_alignments:
        min_r_bond_alignment = np.max(fixed_bond_alignments) < max_fixed_bond_alignment
    else:
        min_r_bond_alignment = True
    
    if reaction_bond_alignments:
        bond_alignment_validation = np.min(reaction_bond_alignments) > min_reaction_bond_alignment and min_r_bond_alignment
    else:
        bond_alignment_validation = min_r_bond_alignment
    
    return match and bond_alignment_validation,match,bond_alignment_validation,reaction_bond_alignments,fixed_bond_alignments