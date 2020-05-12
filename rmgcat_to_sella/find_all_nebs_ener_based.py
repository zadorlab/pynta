import os

import yaml

import numpy as np

from ase.io import read, write
from ase.utils.structure_comparator import SymmetryEquivalenceCheck

from rmgcat_to_sella.find_all_nebs import find_all_nebs_ener_based, create_relax_jobs, generate_unique_combinations
from rmgcat_to_sella.adjacency_to_3d import rmgcat_to_gratoms

from spglib import get_symmetry


def get_all_species_ener_based(path):
    species = []
    for file in os.listdir(path):
        if file.endswith('.xyz'):
            fpath = os.path.join(path, file)
            species.append(read(fpath))
    return species


def generate_unique_combinations_ener_based(slab, species):
    symm = get_symmetry(slab)
    symops = list(zip(symm['rotations'], symm['translations']))
    nslab = len(slab)

    # print(species[0])
    if len(species) == 1:
        combos = species[0]
    elif len(species) == 2:
        combos = []
        for s1full in species[0]:
            s1 = s1full[nslab:].copy()
            # For each unique geometry of the first reactant, generate
            # all rotations and translations allowed by the crystal
            # symmetry of the slab.
            s1scpos = s1.get_scaled_positions()
            for s2full in species[1]:
                s2 = s2full[nslab:].copy()
                for rot, trans in symops:
                    newscpos = s1scpos @ rot.T + trans
                    newscpos %= 1.
                    newscpos %= 1.
                    s1.set_scaled_positions(newscpos)
                    combos.append(slab + s1 + s2)
    else:
        raise RuntimeError(
            "Only know how to handle at most 2 reactants at a time")

    print(f'Found {len(combos)} total combos')

    good_combos = []
    for combo in combos:
        ads = combo[nslab:]
        nads = len(ads)
        # If any of the atoms are overlapping, then forget
        # this structure.
        clashing = False
        if np.all(ads.get_all_distances(mic=True)[np.triu_indices(nads, k=1)] > 0.7) and np.all(ads.get_all_distances(mic=True)[np.triu_indices(nads, k=1)] < 3.0):
            good_combos.append(combo)

    print(f'Found {len(good_combos)} good (non-clashing) combos')

    good_unique_combos = []
    compare = []
    comparator = SymmetryEquivalenceCheck()
    for s1 in good_combos:
        tmp1 = s1.copy()
        tmp1.pbc = True
        if not comparator.compare(tmp1, compare):
            good_unique_combos.append(s1)
            compare.append(tmp1)
    print(f'Found {len(good_unique_combos)} combos that are both good and unique')

    return good_unique_combos


def find_all_nebs_ener_based(slab, yamlfile, facetpath):
    with open(yamlfile, 'r') as f:
        text = f.read()
    reactions = yaml.safe_load(text)

    # Symmetry operations for the bare slab
    symm = get_symmetry(slab)
    symops = list(zip(symm['rotations'], symm['translations']))

    species_unique = dict()
    nslab = len(slab)

    for rxn in reactions:
        reactants, _ = rmgcat_to_gratoms(rxn['reactant'].split('\n'))
        products, _ = rmgcat_to_gratoms(rxn['product'].split('\n'))
        if len(reactants) > 2 or len(products) > 2:
            raise RuntimeError(
                "Only know how to handle at most 2 reactants at a time")

        r_unique = []
        p_unique = []
        for rp, uniquelist in ((reactants, r_unique), (products, p_unique)):
            for species in rp:
                symbols = str(species.symbols)
                speciesdir = os.path.join(
                    facetpath, 'minima_unique_ener_based', symbols)
                if symbols not in species_unique:
                    species_unique[symbols] = get_all_species_ener_based(
                        speciesdir)
                uniquelist.append(species_unique[symbols])

        r_name = '+'.join([str(species.symbols) for species in reactants])
        p_name = '+'.join([str(species.symbols) for species in products])

        rxn_name = r_name + '_' + p_name

        unique_reactants = generate_unique_combinations_ener_based(
            slab, r_unique)
        unique_products = generate_unique_combinations_ener_based(
            slab, p_unique)

        rpath = os.path.join(
            facetpath, 'rxns_ener_based_const', rxn_name, 'reactants')
        ppath = os.path.join(
            facetpath, 'rxns_ener_based_const', rxn_name, 'products')

        os.makedirs(rpath, exist_ok=True)
        os.makedirs(ppath, exist_ok=True)

        for i, species in enumerate(unique_reactants):
            fname = os.path.join(rpath, '{}.xyz'.format(str(i).zfill(3)))
            write(fname, species)
            fname = os.path.join(rpath, '{}.png'.format(str(i).zfill(3)))
            write(fname, species)
        for i, species in enumerate(unique_products):
            fname = os.path.join(ppath, '{}.xyz'.format(str(i).zfill(3)))
            write(fname, species)
            fname = os.path.join(ppath, '{}.png'.format(str(i).zfill(3)))
            write(fname, species)
