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
