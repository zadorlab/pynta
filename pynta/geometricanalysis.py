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
