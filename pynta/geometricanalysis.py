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
