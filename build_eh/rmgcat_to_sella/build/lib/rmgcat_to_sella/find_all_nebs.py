import os
import itertools

import numpy as np

from ase.geometry import find_mic
from ase.io import read, write
from ase.utils.structure_comparator import SymmetryEquivalenceCheck

import yaml

from spglib import get_symmetry
from .adjacency_to_3d import rmgcat_to_gratoms

def get_all_species(path):
    idx = 0
    species = []
    while True:
        fname = '{}.xyz'.format(str(idx).zfill(2))
        fpath = os.path.join(path, fname)
        if os.path.isfile(fpath):
            species.append(read(fpath))
        else:
            return species
        idx += 1

def generate_unique_combinations(slab, species):
    symm = get_symmetry(slab)
    symops = list(zip(symm['rotations'], symm['translations']))
    nslab = len(slab)

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
        raise RuntimeError("Only know how to handle at most 2 reactants at a time")

    print(f'Found {len(combos)} total combos')

    good_combos = []
    for combo in combos:
        ads = combo[nslab:]
        nads = len(ads)
        # If any of the atoms are overlapping, then forget
        # this structure.
        clashing = False
        if np.all(ads.get_all_distances(mic=True)[np.triu_indices(nads, k=1)] > 0.5):
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


def find_all_nebs(slab, reactionlist, facetpath):
    with open(reactionlist, 'r') as f:
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
            raise RuntimeError("Only know how to handle at most 2 reactants at a time")

        r_unique = []
        p_unique = []
        for rp, uniquelist in ((reactants, r_unique), (products, p_unique)):
            for species in rp:
                symbols = str(species.symbols)
                speciesdir = os.path.join(facetpath, 'minima_unique', symbols)
                if symbols not in species_unique:
                    species_unique[symbols] = get_all_species(speciesdir)
                uniquelist.append(species_unique[symbols])

        r_name = '+'.join([str(species.symbols) for species in reactants])
        p_name = '+'.join([str(species.symbols) for species in products])

        rxn_name = r_name + '_' + p_name

        unique_reactants = generate_unique_combinations(slab, r_unique)
        unique_products = generate_unique_combinations(slab, p_unique)
        rpath = os.path.join(facetpath, 'rxns', rxn_name, 'reactants')
        ppath = os.path.join(facetpath, 'rxns', rxn_name, 'products')
        os.makedirs(rpath, exist_ok=True)
        os.makedirs(ppath, exist_ok=True)
        for i, species in enumerate(unique_reactants):
            fname = os.path.join(rpath, '{}.xyz'.format(str(i).zfill(3)))
            write(fname, species)
        for i, species in enumerate(unique_products):
            fname = os.path.join(ppath, '{}.xyz'.format(str(i).zfill(3)))
            write(fname, species)


def create_relax_jobs(basepath, pytemplate, shtemplate=None):
    with open(pytemplate, 'r') as f:
        pytemplate = f.read()

    if shtemplate is not None:
        with open(shtemplate, 'r') as f:
            shtemplate = f.read()

    for rxn in os.listdir(basepath):
        rxnpath = os.path.join(basepath, rxn)
        if not os.path.isdir(rxnpath):
            continue
        for rp in os.listdir(rxnpath):
            rppath = os.path.join(rxnpath, rp)
            if not os.path.isdir(rppath):
                continue
            for structure in os.listdir(rppath):
                if structure.endswith('xyz'):
                    prefix = os.path.splitext(structure)[0]
                    fname = os.path.join(basepath, f'{rxn}_{rp}_{prefix}_relax.py')
                    with open(fname, 'w') as f:
                        f.write(pytemplate.format(rxn=rxn,
                                                  rp=rp,
                                                  prefix=prefix))
                    if shtemplate is None:
                        continue

                    fname = os.path.join(facetpath, f'{rxn}_{rp}_{prefix}_relax.sh')
                    with open(fname, 'w') as f:
                        f.write(shtemplate.format(rxn=rxn,
                                                  rp=rp,
                                                  prefix=prefix))


def set_up_nebs(yamlfile, facetpath, slab):
    with open(yamlfile, 'r') as f:
        yamltxt = f.read()
    rxn_spec = yaml.safe_load(yamltxt)

    nslab = len(slab)
    symm = get_symmetry(slab)
    symops = list(zip(symm['translations'], symm['rotations']))

    all_rxns = dict()
    for rxn in rxn_spec:
        reactants, _ = rmgcat_to_gratoms(rxn['reactant'].split('\n'))
        r_names = '+'.join([str(species.symbols) for species in reactants])
        products, _ = rmgcat_to_gratoms(rxn['product'].split('\n'))
        p_names = '+'.join([str(species.symbols) for species in products])
        rxn_name = r_names + '_' + p_names
        all_rxns[rxn_name] = (reactants, products)

    for name, rxn in all_rxns.items():
        rxnpath = os.path.join(facetpath, 'rxns_unique', name)
        rpath = os.path.join(rxnpath, 'reactants')

        reactants = []
        for fname in os.listdir(rpath):
            if fname.endswith('.xyz'):
                reactants.append(read(os.path.join(rpath, fname)))

        ppath = os.path.join(rxnpath, 'products')
        products = []
        for fname in os.listdir(ppath):
            if fname.endswith('.xyz'):
                products.append(read(os.path.join(ppath, fname)))

        tmpatoms = reactants[0][nslab:].copy()
        atomic_numbers = sorted(list(set(tmpatoms.numbers)))
        all_atomic_numbers = list(tmpatoms.numbers)
        n_element = []
        for num in atomic_numbers:
            n_element.append(all_atomic_numbers.count(num))

        # Get all reactants and products into a consistent atom order

        new_reactants = []
        for react in reactants:
            rads = react[nslab:]
            el_indices = {num: [] for num in atomic_numbers}
            for i, atom in enumerate(rads):
                el_indices[atom.number].append(i)
            new_indices = []
            for num in atomic_numbers:
                new_indices += el_indices[num]
            #react[nslab:] = rads[new_indices]
            new_reactants.append(react[:nslab] + rads[new_indices])
        reactants = new_reactants

        new_products = []
        for prod in products:
            pads = prod[nslab:]
            el_indices = {num: [] for num in atomic_numbers}
            for i, atom in enumerate(pads):
                el_indices[atom.number].append(i)
            new_indices = []
            for num in atomic_numbers:
                new_indices += el_indices[num]
            #prod[nslab:] = pads[new_indices]
            new_products.append(prod[:nslab] + pads[new_indices])
        products = new_products

        permute = []
        start = 0
        for num, nel in zip(atomic_numbers, n_element):
            permute.append(list(itertools.permutations(range(start, start+nel))))
            start += nel

        all_reactions = []
        for react, prod in itertools.product(reactants, products):
            rads = react[nslab:]
            pads = prod[nslab:]
            shortest_dist = np.linalg.norm(pads.positions - rads.positions)
            shortest = (list(range(len(rads))), np.zeros(3, dtype=int), np.eye(3, dtype=int))
            for trans, rot in symops:
                tmp = pads.copy()
                scpos = tmp.get_scaled_positions()
                scpos = scpos @ rot.T + trans
                scpos %= 1.
                scpos %= 1.
                tmp.set_scaled_positions(scpos)
                for permutation in itertools.product(*permute):
                    new_indices = np.concatenate(permutation)
                    tmp2 = tmp[new_indices]
                    dxs, _ = find_mic(tmp2.positions - rads.positions, slab.cell, pbc=True)
                    dist = np.linalg.norm(dxs)
                    if dist < shortest_dist:
                        shortest_dist = dist
                        shortest = (new_indices, trans, rot)
            new_indices, trans, rot = shortest
            newprod = prod[:nslab] + pads[new_indices]
            scpos = newprod.get_scaled_positions()
            scpos = scpos @ rot.T + trans
            scpos %= 1.
            scpos %= 1.
            newprod.set_scaled_positions(scpos)

            rslab = react[:nslab]
            pslab = newprod[:nslab]

            slab_indices = -np.ones(nslab, dtype=int)
            for i, atom in enumerate(pslab):
                _, dxs = find_mic(atom.position - rslab.positions, rslab.cell, pbc=True)
                new_index = np.argmin(dxs)
                assert slab_indices[new_index] == -1
                slab_indices[new_index] = i

            newprod = pslab[slab_indices] + newprod[nslab:]
            newprod.pbc = react.pbc
            # now, reposition atoms in the product so they are close to the corresponding atoms
            # in the reactant
            newpos = np.zeros_like(newprod.positions)
            for i, (ratom, patom) in enumerate(zip(react, newprod)):
                dx, _ = find_mic((patom.position - ratom.position).reshape((1, 3)), slab.cell, pbc=True)
                newpos[i] = ratom.position + dx
            newprod.set_positions(newpos)

            all_reactions.append((react, newprod))

        for i, (react, prod) in enumerate(all_reactions):
            prefix = str(i).zfill(3)
            nebpath = os.path.join(facetpath, 'nebs', name, prefix)
            os.makedirs(nebpath, exist_ok=True)
            write(os.path.join(nebpath, 'reactants.xyz'), react)
            write(os.path.join(nebpath, 'products.xyz'), prod)

