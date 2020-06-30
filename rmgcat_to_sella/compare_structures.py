#!/usr/bin/env python3
import yaml
import os
import networkx as nx
from ase.io import read, write
from ase.data import covalent_radii

from catkit import Gratoms

import numpy as np
from spglib import get_symmetry
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from pathlib import Path

from .adjacency_to_3d import rmgcat_to_gratoms
from .graph_utils import node_test

def get_unique(adsorbates):
    unique = []
    comparator = SymmetryEquivalenceCheck()
    for atoms1 in adsorbates:
        tmp1 = atoms1.copy()
        tmp1.pbc = True
        for atoms2 in unique:
            tmp2 = atoms2.copy()
            tmp2.pbc = True
            if comparator.compare(tmp1, tmp2):
                break
        else:
            unique.append(atoms1)
    return unique


def get_unique_TS_estimate_before_xtb(facetpath, TSestimateDir):
    good_minima = []
    result_list = []
    unique_index = []
    gpath = os.path.join(facetpath, TSestimateDir)
    # for geom in sorted(os.listdir(gpath), key=str):
    #     geomDir = os.path.join(gpath, geom)
    #     if os.path.isdir(geomDir):
    trajlist = sorted(Path(gpath).glob('*xyz'), key = str)
    for traj in trajlist:
        print(traj)
        minima = read(traj)
        minima.pbc = True
        comparator = SymmetryEquivalenceCheck()
        result = comparator.compare(minima, good_minima)
        result_list.append(result)
        if result is False:
            good_minima.append(minima)
    for num, res in enumerate(result_list):
        if res is False:
            unique_index.append(str(num).zfill(3))
    return unique_index


def get_unique_minima_after_opt(facetpath, minimaDir):
    good_minima = []
    result_list = []
    unique_index = []
    gpath = os.path.join(facetpath, minimaDir)
    # for geom in sorted(os.listdir(gpath), key=str):
    #     geomDir = os.path.join(gpath, geom)
    #     if os.path.isdir(geomDir):
    trajlist = sorted(Path(gpath).glob('*final.xyz'), key = str)
    for traj in trajlist:
        print(traj)
        minima = read(traj)
        minima.pbc = True
        comparator = SymmetryEquivalenceCheck()
        result = comparator.compare(minima, good_minima)
        result_list.append(result)
        if result is False:
            good_minima.append(minima)
    for num, res in enumerate(result_list):
        if res is False:
            unique_index.append(str(num).zfill(3))
    return unique_index


def get_unique_minima_eh(yamlfile, facetpath, slab):
    with open(yamlfile, 'r') as f:
        yamltxt = f.read()
    rxn_spec = yaml.safe_load(yamltxt)

    all_species = []
    for rxn in rxn_spec:
        grtmp1, _ = rmgcat_to_gratoms(rxn['reactant'].split('\n'))
        grtmp2, _ = rmgcat_to_gratoms(rxn['product'].split('\n'))

        # for sp1 in grtmp1 + grtmp2:
        for sp1 in grtmp2:
            for sp2 in all_species:
                if nx.is_isomorphic(sp1.graph, sp2.graph, node_test):
                    break
            else:
                all_species.append(sp1)

    for species in all_species:
        path = os.path.join(facetpath, 'minima', str(species.symbols))
        fnames = os.listdir(path)
        all_confs = []
        for fname in sorted(fnames, key=str):
            if fname.endswith('.traj'):
                all_confs.append(read(os.path.join(path, fname), index=-1))
        unique_confs = find_all_unique(all_confs, slab, species)
        newpath = os.path.join(facetpath, 'minima_unique', str(species.symbols))
        os.makedirs(newpath, exist_ok=True)
        for i, conf in enumerate(unique_confs):
            fname_xyz = os.path.join(newpath, '{}.xyz'.format(str(i).zfill(2)))
            fname_png = os.path.join(newpath, '{}.png'.format(str(i).zfill(2)))
            write(fname_xyz, conf)
            write(fname_png, conf)

def get_unique_rxns(yamlfile, facetpath, slab):
    with open(yamlfile, 'r') as f:
        yamltxt = f.read()
    rxn_spec = yaml.safe_load(yamltxt)

    all_rxns = dict()
    for rxn in rxn_spec:
        # assigning first value, i.e, rmgcat_to_gratoms(rxn['reactant'].split('\n') to reactants and ignoring the second, i.e. None. In the end reatcants is a catkit gratom object
        reactants, _ = rmgcat_to_gratoms(rxn['reactant'].split('\n'))
        r_names = '+'.join([str(species.symbols) for species in reactants])
        products, _ = rmgcat_to_gratoms(rxn['product'].split('\n'))
        p_names = '+'.join([str(species.symbols) for species in products])
        rxn_name = r_names + '_' + p_names
        all_rxns[rxn_name] = (reactants, products)

    for name, rxn in all_rxns.items():
        # split name of rxns to reactant and product. Note, that loop has to be through keys and value
        r_name, p_name = name.split('_')
        rpath = os.path.join(facetpath, 'rxns', name, 'reactants')
        all_reactants = []
        for fname in os.listdir(rpath):
            if fname.endswith('.traj'):
                all_reactants.append(read(os.path.join(rpath, fname), index=-1))
        unique_reactants = find_all_unique(all_reactants, slab, rxn[0])

        ppath = os.path.join(facetpath, 'rxns', name, 'products')
        all_products = []
        for fname in os.listdir(ppath):
            if fname.endswith('.traj'):
                all_products.append(read(os.path.join(ppath, fname), index=-1))
        unique_products = find_all_unique(all_products, slab, rxn[1])

        newrpath = os.path.join(facetpath, 'rxns_unique', name, 'reactants')
        os.makedirs(newrpath, exist_ok=True)
        for i, conf in enumerate(unique_reactants):
            fname = os.path.join(newrpath, '{}.xyz'.format(str(i).zfill(2)))
            write(fname, conf)

        newppath = os.path.join(facetpath, 'rxns_unique', name, 'products')
        os.makedirs(newppath, exist_ok=True)
        for i, conf in enumerate(unique_products):
            fname = os.path.join(newppath, '{}.xyz'.format(str(i).zfill(2)))
            write(fname, conf)

def find_all_unique(adsorbates, slab, gr_ref):
    # adsorbates is a gratom object containig info about adsorbated specias + surface 
    # slab is a gratom object
    # g_ref is a gratom object. Is a set of connectivity (graph) representation of species. See species definition in get_unique_minima()
    if isinstance(gr_ref, Gratoms):
        gr_ref = [gr_ref]

    nslab = len(slab)
    valid_adsorbates = []
    for adsorbate in adsorbates:
        # Remove slab atoms, look only at adsorbate
        ads = adsorbate[nslab:]
        # Find the connectivity between atoms
        edges = []
        nads = len(ads)
        for i, ai in enumerate(ads):
            for j, aj in enumerate(ads):
                if j <= i:
                    continue
                rij = ads.get_distance(i, j, mic=True)
                rcut = 1.25 * (covalent_radii[ai.number] + covalent_radii[aj.number])
                if rij < rcut:
                    edges.append((i, j))
        #print(edges)
        grads = Gratoms(ads.symbols, edges=edges)
        grads_separate = []
        # Identify unique subgraphs, which are distinct but co-adsorbed species
        for i, subgraph in enumerate(nx.connected_component_subgraphs(grads.graph)):
            indices = list(subgraph.nodes)
            symbols = grads[indices].symbols
            new_indices = {old:new for new, old in enumerate(indices)}
            new_edges = []
            for edge in subgraph.edges:
                newa = new_indices[edge[0]]
                newb = new_indices[edge[1]]
                new_edges.append((newa, newb))
            new_gratoms = Gratoms(symbols, positions=grads[indices].positions, cell=grads.cell,
                                  pbc=grads.pbc, edges=new_edges)
            grads_separate.append(new_gratoms)
        # If the number of co-adsorbates doesn't match the reference, then
        # this structure isn't valid
        if len(gr_ref) != len(grads_separate):
            continue
        fail = False
        # Make sure each adsorbate has a matching species in the reference
        gr_ref_tmp = gr_ref.copy()
        #print(grads_separate, gr_ref)
        for ads in grads_separate:
            for ref in gr_ref_tmp:
                if nx.is_isomorphic(ads.graph, ref.graph, node_test):
                    gr_ref_tmp.remove(ref)
                    break
            else:
                fail = True
                break
        if fail:
            continue
        valid_adsorbates.append(adsorbate)
    unique_adsorbates = get_unique(valid_adsorbates)
    return unique_adsorbates
