from catkit import Gratoms
from catkit.gen.adsorption import AdsorptionSites, Builder
from catkit.build import molecule

from .adjacency_to_3d import get_edges, rmgcat_to_gratoms
from .find_all_nebs import get_all_species
from .graph_utils import node_test

from ase.io import read, write
from ase import Atoms

import numpy as np

import yaml

import os

import shutil

def estimate_struc_to_neb(slab, repeats, yamlfile, facetpath):
    with open(yamlfile, 'r') as f:
        yamltxt = f.read()
    reactions = yaml.safe_load(yamltxt)

    species_unique = dict()
    nslab = len(slab)

    speciesInd = []
    bonds = []

    for rxn in reactions:
        # transforming reactions data to gratom objects
        reactants, rbonds = rmgcat_to_gratoms(rxn['reactant'].split('\n'))
        products, pbonds = rmgcat_to_gratoms(rxn['product'].split('\n'))
        speciesInd += reactants + products
        bonds += rbonds + pbonds
        # print(pbonds)
        # print(products)

        r_unique = []
        p_unique = []

    symbols_list = []

    for rp, uniquelist in ((reactants, r_unique), (products, p_unique)):
        for species in rp:
            symbols = str(species.symbols)
            symbols_list.append(symbols)
            speciesdir = os.path.join(facetpath, 'minima_unique', symbols)
            if symbols not in species_unique:
                species_unique[symbols] = get_all_species(speciesdir)
            uniquelist.append(species_unique[symbols])

    r_name = '+'.join([str(species.symbols) for species in reactants])
    p_name = '+'.join([str(species.symbols) for species in products])

    rxn_name = r_name + '_' + p_name

    print(reactants)
    unique_species = []
    unique_bonds = []
    images = []

    # uniqueDir = os.path.join(facetpath, 'minima_unique')

    r_name_list = [str(species.symbols) for species in reactants]
    p_name_list = [str(species.symbols) for species in products]

    print(p_name_list)


    saveDir = f'./{facetpath}/neb_estimate/{rxn_name}/products'
    if os.path.exists(saveDir):
        shutil.rmtree(saveDir)
        os.makedirs(saveDir)
    else:
        os.makedirs(saveDir)
    # os.makedirs(saveDir)

    # ADSORBATES
    # struc1_tmp = molecule('CH2')[0]
    # struc2_tmp = molecule('H')[0]
    struc1_tmp = molecule(p_name_list[0])[0]
    struc2_tmp = molecule(p_name_list[1])[0]

    pos1 = struc1_tmp.get_positions() 
    pos2 = struc2_tmp.get_positions() 

    # Moving x coordinate of one of co-adsorbed species by defined value dist
    dist = 3.0

    pos2_updated = [[pos1[0][0] + dist, pos1[0][1], pos1[0][2]]]

    struc1 = Gratoms(p_name_list[0], positions = pos1)
    struc2 = Gratoms(p_name_list[1], positions = pos2_updated)

    combined = struc1 + struc2

    #slab transfromed to gratom object

    slabedges, tags = get_edges(slab, True)


    grslab = Gratoms(numbers=slab.numbers,
                     positions=slab.positions,
                     cell=slab.cell,
                     pbc=slab.pbc,
                     edges=slabedges)
    grslab.arrays['surface_atoms'] = tags

    #building adsorbtion structures

    ads_builder = Builder(grslab)

    structs = ads_builder.add_adsorbate(combined, [0], -1, auto_construct=False)

    big_slab = slab * repeats
    nslab = len(slab)

    for i, struc in enumerate(structs):
        big_slab_ads = big_slab + struc[nslab:]
        write(os.path.join(saveDir, '{}'.format(str(i).zfill(2)) + '_' + p_name + '.png'), big_slab_ads)
        write(os.path.join(saveDir, '{}'.format(str(i).zfill(2)) + '_' + p_name + '.xyz'), big_slab_ads)


