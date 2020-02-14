#!/usr/bin/env python

import os
import yaml

import numpy as np
import networkx as nx

from ase import Atoms
from ase.io import read, write
from ase.data import covalent_radii
from ase.build import bulk, fcc111

from catkit.gen.utils import connectivity_to_edges, get_cutoff_neighbors, to_gratoms
from catkit.gen.molecules import get_3D_positions
from catkit.gen.adsorption import AdsorptionSites, Builder
from catkit.gen.surface import SlabGenerator
from catkit import Gratoms

from .graph_utils import node_test

# Instead of using CatKit's built in slab generator routines, we want to
# use pre-relaxed slab geometries to save computer time. In order to use
# our own geometries, we need to generate connectivities in a way that
# CatKit understands. CatKit has some built-in routines for doing this,
# but they are not great, so I wrote my own method. This can be
# used to convert any ASE Atoms object into a CatKit Gratoms object,
# though we are only currently using it for the catalyst geometry.
def get_edges(atoms, find_surface=False):
    # If the Atoms object is periodic, we need to check connectivity
    # across the unit cell boundary as well.
    tvecs = np.array([[0., 0., 0.]])
    if np.any(atoms.pbc):
        cell = atoms.cell
        # We are looking for atoms that are at most 2.5 times the
        # greatest covalent radius of all atoms in the system, so we
        # limit the number of periodic images we consider in each
        # direction to those that are within this cutoff radius of any
        # point on the face of the unit cell
        cutoff = max(covalent_radii[atoms.numbers]) * 2.5

        # If you want to understand the code below, read the `find_mic`
        # code in ASE, located at ase/geometry/geometry.py
        latt_len = np.sqrt((cell**2).sum(1))
        V = atoms.get_volume()
        padding = atoms.pbc * np.array(np.ceil(cutoff * np.prod(latt_len) /
                                               (V * latt_len)), dtype=int)
        offsets = np.mgrid[-padding[0]:padding[0]+1,
                           -padding[1]:padding[1]+1,
                           -padding[2]:padding[2]+1].T
        tvecs = np.dot(offsets, cell).reshape(-1, 3)

    natoms = len(atoms)
    edges = []
    if find_surface:
        pairvecs = np.zeros_like(atoms.positions)
    for atomi in atoms:
        for atomj in atoms:
            i = atomi.index
            j = atomj.index
            # Like above, consider only bonds where j >= i. Note, we do
            # need to consider bonds where j == i because of periodic
            # boundary conditions.
            if j < i:
                continue
            # 1.25 times the sum of the covalent radii was chosen based
            # on trial and error. Too small, and you miss neighbors.
            # Too big, and you start including next-nearest-neighbors.
            cutoff = 1.25 * (covalent_radii[atomi.number] + covalent_radii[atomj.number])
            # xij is the direct displacement vector in the central unit
            # cell.
            xij = atomj.position - atomi.position
            nbonds = 0
            # Loop over all neighboring unit cells
            for tvec in tvecs:
                # ...including the central unit cell. If i == j, then
                # explicitly skip the central unit cell.
                if i == j and np.all(tvec == [0., 0., 0.]):
                    continue
                # Count up the number of times i bonds to j via pbc
                if np.linalg.norm(xij + tvec) < cutoff:
                    if find_surface:
                        dx = xij + tvec
                        dx /= np.linalg.norm(dx)
                        pairvecs[i] -= dx
                        pairvecs[j] += dx
                    nbonds += 1
            # CatKit uses NetworkX "MultiGraph"s for periodic systems.
            # I'm not entirely sure how to interpret the edges for a
            # multigraph, but this is how CatKit wants multiple edges
            # between the same two atoms to be specified.
            for k in range(nbonds):
                edges.append((i, j, k))
    if not find_surface:
        return edges

    surface = np.zeros(len(atoms), dtype=int)
    for i, pairvec in enumerate(pairvecs):
        # Big pairvec means highly asymmetrical atoms, which
        # implies that it is near or at a surface
        if np.linalg.norm(pairvec) > 1.:
            # Use sign to determine whether it is at the top or
            # the bottom of the slab
            surface[i] = int(np.sign(pairvec[2]))
    return edges, surface

def rmgcat_to_gratoms(adjtxt):
    symbols = []
    edges = []
    tags = []
    bond_index = None
    for i, line in enumerate(adjtxt):
        if i == 0:
            continue
        if not line:
            break

        line = line.split()
        inc = 0
        if line[1][0] == '*':
            inc = 1
            tags.append(int(line[1][1]))
        else:
            tags.append(0)

        symbols.append(line[1 + inc])
        conn = line[5 + inc:]

        for bond in conn:
            j = int(bond.strip('{}').split(',')[0])
            if j > i:
                edges.append((i - 1, j - 1))

    gratoms = Gratoms(symbols, edges=edges)

    del_indices = []

    for i, atom in enumerate(gratoms):
        if atom.symbol == 'X':
            for j in gratoms.graph.neighbors(i):
                tags[j] *= -1
            del_indices.append(i)

    gratoms.set_tags(tags)
    del gratoms[del_indices]

    gratoms_list = []
    bonds = []
    for i, subgraph in enumerate(nx.connected_component_subgraphs(gratoms.graph)):
        indices = list(subgraph.nodes)
        symbols = gratoms[indices].symbols
        #new_gratoms = gratoms[indices].copy()
        new_indices = {old:new for new, old in enumerate(indices)}
        new_edges = []
        for edge in subgraph.edges:
            newa = new_indices[edge[0]]
            newb = new_indices[edge[1]]
            new_edges.append((newa, newb))
        new_gratoms = Gratoms(symbols, edges=new_edges)

        bond = None
        tags = new_gratoms.get_tags()
        for i, tag in enumerate(tags):
            if tag < 0:
                if bond is None:
                    bond = [i]
                elif len(bond) == 1:
                    bond.append(i)
                else:
                    raise RuntimeError('At most two bonds to the metal are allowed per adsorbate!')
                tags[i] = abs(tags[i])
        new_gratoms.set_tags(tags)
        bonds.append(bond)
        gratoms_list.append(new_gratoms)

    return gratoms_list, bonds

def adjacency_to_3d(reactionlist, slab, repeats, slabname):
    with open(reactionlist, 'r') as f:
        text = f.read()
    reactions = yaml.safe_load(text)

    species = []
    bonds = []
    for rxn in reactions:
        reactants, rbonds = rmgcat_to_gratoms(rxn['reactant'].split('\n'))
        products, pbonds = rmgcat_to_gratoms(rxn['product'].split('\n'))
        species += reactants + products
        bonds += rbonds + pbonds

    unique_species = []
    unique_bonds = []
    images = []
    # check if any products are the same as any reactants
    for species1, bond in zip(species, bonds):
        for species2 in unique_species:
            if nx.is_isomorphic(species1.graph, species2.graph, node_test):
                break
        else:
            images.append(get_3D_positions(species1))
            unique_species.append(species1)
            unique_bonds.append(bond)

    slabedges, tags = get_edges(slab, True)
    grslab = Gratoms(numbers=slab.numbers,
                     positions=slab.positions,
                     cell=slab.cell,
                     pbc=slab.pbc,
                     edges=slabedges)
    grslab.arrays['surface_atoms'] = tags

    ads_builder = Builder(grslab)

    structures = dict()

    for adsorbate, bond in zip(images, unique_bonds):
        if len(adsorbate) == 0:
            continue

        if bond is None:
            bond = [0]

        key = adsorbate.get_chemical_formula()
        try:
            structs = ads_builder.add_adsorbate(adsorbate, bonds=bond, index=-1)
            structures[key] = structs
        except IndexError:
            print(adsorbate, adsorbate.edges, adsorbate.get_tags())

    big_slab = slab * repeats
    nslab = len(slab)

    for species_name, adsorbate in structures.items():
        savedir = f'./{slabname}/{species_name}'
        os.makedirs(savedir, exist_ok=True)
        for j, structure in enumerate(adsorbate):
            big_slab_ads = big_slab + structure[nslab:]
            write(os.path.join(savedir, '{}.xyz'.format(str(j).zfill(2))), big_slab_ads)
            write(os.path.join(savedir, '{}.png'.format(str(j).zfill(2))), big_slab_ads)
            # write(os.path.join(savedir, '{}.png'.format(str(j).zfill(2))), big_slab_ads, rotation='10z,-80x')
