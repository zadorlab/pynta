#!/usr/bin/env python
import os
from pathlib import PosixPath
from typing import Tuple, List, Dict

import numpy as np
import networkx as nx

from ase.io import read, write
from ase.data import covalent_radii

from pynta.utils import edge_cases_bonded_dict, edge_cases_topology_dict
from pynta.excatkit.molecule import Molecule
from pynta.excatkit.adsorption import Builder
from pynta.excatkit.gratoms import Gratoms
from pynta.io import IO
from pynta.graph_utils import node_test

# Instead of using CatKit's built in slab generator routines, we want to
# use pre-relaxed slab geometries to save computer time. In order to use
# our own geometries, we need to generate connectivities in a way that
# CatKit understands. CatKit has some built-in routines for doing this,
# but they are not great, so I wrote my own method. This can be
# used to convert any ASE Atoms object into a CatKit Gratoms object,
# though we are only currently using it for the catalyst geometry.


class Adsorbates:
    ''' This class handles adsorbates and puts them on the surface'''

    def __init__(
            self,
            facetpath: str,
            slab: str,
            repeats: Tuple[int, int, int],
            yamlfile: str,
            creation_dir: PosixPath) -> None:
        ''' Initializing

        Parameters:
        ___________
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        slab : str
            a '.xyz' file name with the optimized slab
            e.g.
            'Cu_111_slab_opt.xyz'
        repeats: tuple
            specify reapeats in (x, y, z) direction,
            eg. (3, 3, 1)
        yamlfile : str
            a name of the .yaml file with reaction list

        '''
        self.facetpath = facetpath
        self.slab = slab
        self.repeats = repeats
        self.creation_dir = creation_dir
        self.yamlfile = os.path.join(self.creation_dir, yamlfile)

    def get_edges(
            self,
            find_surface: bool = False) -> Tuple[List[str], np.ndarray]:
        ''' Get adsorption edges

        Parameters:
        ___________
        find_surface : bool
            specify whether to include surface or not
            default = False

        Returns:
        ________
        edges : list[tuple]
            adsobrtion edges
        surface : numpy.ndarray
            an array with tags describing:
            top surface atoms (1)
            bottom surface (-1)
            and bulk atoms (0)
            Atoms with tags '1' are considered as the possible binding spots

        '''
        # read slab as an Atom object
        slab_atom = read(os.path.join(self.creation_dir, self.slab))

        # If the Atoms object is periodic, we need to check connectivity
        # across the unit cell boundary as well.
        tvecs = np.array([[0., 0., 0.]])
        if np.any(slab_atom.pbc):
            cell = slab_atom.cell
            # We are looking for atoms that are at most 2.5 times the
            # greatest covalent radius of all atoms in the system, so we
            # limit the number of periodic images we consider in each
            # direction to those that are within this cutoff radius of any
            # point on the face of the unit cell
            cutoff = max(covalent_radii[slab_atom.numbers]) * 2.5

            # If you want to understand the code below, read the `find_mic`
            # code in ASE, located at ase/geometry/geometry.py
            latt_len = np.sqrt((cell**2).sum(1))
            V = slab_atom.get_volume()
            padding = slab_atom.pbc * np.array(np.ceil(
                cutoff * np.prod(latt_len) / (V * latt_len)),
                dtype=int
            )
            offsets = np.mgrid[-padding[0]:padding[0] + 1,
                               -padding[1]:padding[1] + 1,
                               -padding[2]:padding[2] + 1].T
            tvecs = np.dot(offsets, cell).reshape(-1, 3)

        edges = []
        if find_surface:
            pairvecs = np.zeros_like(slab_atom.positions)
        for atomi in slab_atom:
            for atomj in slab_atom:
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
                cutoff = 1.25 * (
                    covalent_radii[atomi.number] + covalent_radii[atomj.number]
                )
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

        surface = np.zeros(len(slab_atom), dtype=int)
        for i, pairvec in enumerate(pairvecs):
            # Big pairvec means highly asymmetrical atoms, which
            # implies that it is near or at a surface
            if np.linalg.norm(pairvec) > 1.:
                # Use sign to determine whether it is at the top or
                # the bottom of the slab
                surface[i] = int(np.sign(pairvec[2]))
        return edges, surface

    # def get_how_many_atoms_before_adsorbate(self, adjtxt):
    #     n_surf_at_befor_ads = 0

    #     for line in adjtxt:
    #         if 'multiplicity' in line:
    #             continue
    #         if 'X' in line:
    #             n_surf_at_befor_ads += 1
    #         else:
    #             if n_surf_at_befor_ads != 0:
    #                 n_surf_at_befor_ads -= 1
    #             break
    #     return n_surf_at_befor_ads

    # def get_proper_atomic_indicies(self, adjtxt):
    #     correct_atomic_indicies = []

    #     remove_one_more = 0
    #     for i, line in enumerate(adjtxt):
    #         if 'multiplicity' in line:
    #             remove_one_more = 1
    #             continue
    #         if 'X' not in line:
    #             n = self.get_how_many_atoms_before_adsorbate(adjtxt)
    #             adsorbate_atom_idx = i - n - remove_one_more
    #             correct_atomic_indicies.append(adsorbate_atom_idx)
    #     return correct_atomic_indicies

    def rmgcat_to_gratoms(self, adjtxt):
        ''' Convert a slice of .yaml file to Catkit's Gratoms object

        Parameters:
        ___________

        adjtxt : list
            a list with a connectivity info for reactant or product
            as from the .yaml file.
            e.g. for given reaction (reactant or product)

            In .yaml file we have something like that:

                    multiplicity -187
                1 *1 C u0 p0 c0 { 2,S} {4,S}
                2    O u0 p0 c0 {1,S}
                3 *2 H u0 p0 c0 {5,S}
                4 *3 X u0 p0 c0 {1,S}
                5 *4 X u0 p0 c0 {3,S}

            but we need here a list like that:

            ['multiplicity -187', '1 *1 C u0 p0 c0 {2,S} {4,S}',
            '2    O u0 p0 c0 {1,S}', '3 *2 H u0 p0 c0 {5,S}',
            '4 *3 X u0 p0 c0 {1,S}', '5 *4 X u0 p0 c0 {3,S}', '']

            So it can be simply converted using the following:

            yamlfile = 'reactions.yaml'
            with open(yamlfile, 'r') as f:
                text = f.read()
            reactions = yaml.safe_load(text)
            for rxn in reactions:
                adjtxt = rxn['reactant'].split('\n')

        Returns:
        ________
        gratoms_list : list
            a Gratom like object
        bonds : list
            a list of bonds to the metal

        '''
        symbols = []
        edges = []
        tags = []
        # bond_index = None
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
        print(gratoms)

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
            print(indices)
            symbols = gratoms[indices].symbols
            # new_gratoms = gratoms[indices].copy()
            new_indices = {old: new for new, old in enumerate(indices)}
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
                        raise RuntimeError(
                            'At most two bonds to the metal are allowed per adsorbate!')
                    tags[i] = abs(tags[i])
            # tags = self.get_proper_atomic_indicies(adjtxt)
            # print(tags)
            tags = indices
            new_gratoms.set_tags(tags)
            bonds.append(bond)
            gratoms_list.append(new_gratoms)

        return gratoms_list, bonds

    def adjacency_to_3d(self) -> None:
        ''' Place adsorbates on the surface '''

        # # get all species
        # all_species_symbols = IO.get_all_unique_species(self.yamlfile)

        # # convert str to Gratoms and deal with edge cases
        # images = []
        # for sp_symbol in all_species_symbols:
        #     index = 0
        #     # deal with edge cases - which topology to use
        #     if sp_symbol in edge_cases_topology_dict.keys():
        #         index = edge_cases_topology_dict[sp_symbol]
        #     images.append(Molecule().molecule(sp_symbol)[index])

        # slabedges, tags = Adsorbates.get_edges(self, True)

        # # convert .xyz file to Atoms object
        # slab_atom = read(os.path.join(self.creation_dir, self.slab))
        # big_slab = slab_atom * self.repeats
        # nslab = len(slab_atom)

        reactions = IO.open_yaml_file(self.yamlfile)
        all_species_symbols = IO.get_all_unique_species(self.yamlfile)

        species = []
        bonds = []
        for rxn in reactions:
            reactants, rbonds = Adsorbates.rmgcat_to_gratoms(
                self, rxn['reactant'].strip().split('\n'))
            products, pbonds = Adsorbates.rmgcat_to_gratoms(
                self, rxn['product'].strip().split('\n'))
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
                images.append(Molecule().get_3D_positions(species1))
                unique_species.append(species1)
                unique_bonds.append(bond)

        slabedges, tags = Adsorbates.get_edges(self, True)
        slab_atom = read(os.path.join(self.creation_dir, self.slab))
        big_slab = slab_atom * self.repeats
        nslab = len(slab_atom)

        # transform Atoms to Gratoms object
        grslab = Gratoms(numbers=slab_atom.numbers,
                         positions=slab_atom.positions,
                         cell=slab_atom.cell,
                         pbc=slab_atom.pbc,
                         edges=slabedges)
        grslab.arrays['surface_atoms'] = tags

        # prepare surface for placing adsorbates
        ads_builder = Builder(grslab)

        # build adsorbates
        structures = dict()
        for sp_symbol, sp_gratoms in zip(all_species_symbols, images):
            if len(sp_gratoms) == 0:
                continue
            print('---')
            for atom in sp_gratoms:
                print(atom)
            # deal with edge cases - which atom connects to the surface
            bonded = [0]
            if sp_symbol in edge_cases_bonded_dict.keys():
                bonded = [edge_cases_bonded_dict[sp_symbol]]

            try:
                # put adsorbates on the surface
                structs = ads_builder.add_adsorbate(
                    sp_gratoms, index=-1, bonds=bonded)
                structures[str(sp_symbol)] = structs
            except IndexError:
                # TODO an idea to put -1 in adsorbate.get_tags() for
                # surface bonded atom
                print(sp_gratoms, sp_gratoms.edges,
                      sp_gratoms.get_tags())

        for sp_symbol, adsorbate in structures.items():
            # create dir
            savedir = os.path.join(
                self.creation_dir, self.facetpath, 'minima', sp_symbol)
            os.makedirs(savedir, exist_ok=True)

            for prefix, structure in enumerate(adsorbate):
                big_slab_ads = big_slab + structure[nslab:]
                # save adsorbates as .xyz and .png
                write(os.path.join(savedir, '{}.xyz'.format(
                    str(prefix).zfill(2))), big_slab_ads)
                write(os.path.join(savedir, '{}.png'.format(
                    str(prefix).zfill(2))), big_slab_ads)

    def create_relax_jobs(
            self,
            pytemplate: str,
            pseudopotentials: Dict[str, str],
            pseudo_dir: str,
            balsam_exe_settings: Dict[str, int],
            calc_keywords: Dict[str, str],
            shtemplate: str = None) -> None:
        ''' Create a submit scripts

        Parameters:
        __________
        pytemplate: python file
            a template to prepare submission scripts
            for adsorbate+surface minimization
        pseudopotentials: dict(str: str)
            a dictionary with QE pseudopotentials for all species.
            e.g.
            dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                H='H.pbe-kjpaw_psl.1.0.0.UPF',
                O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                C='C.pbe-n-kjpaw_psl.1.0.0.UPF',
                )
        pseudo_dir: str
            a path to the QE's pseudopotentials main directory
            e.g.
            '/home/mgierad/espresso/pseudo'
        balsam_exe_settings: dict{str: int}
            a dictionary with balsam execute parameters(cores, nodes, etc.),
            e.g.
            balsam_exe_settings={'num_nodes': 1,
                                   'ranks_per_node': 48,
                                   'threads_per_rank': 1}
        calc_keywords: dict{str: str}
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            calc_keywords={'kpts': (3, 3, 1), 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}
        shtemplate: .sh file
            optional, a .sh template(not required by the workflow)

        '''
        n_kpts = IO().get_kpoints(self.repeats)
        minimapath = os.path.join(self.creation_dir, self.facetpath, 'minima')
        with open(pytemplate, 'r') as f:
            pytemplate = f.read()

        if shtemplate is not None:
            with open(shtemplate, 'r') as f:
                shtemplate = f.read()

        for species in os.listdir(minimapath):
            speciespath = os.path.join(minimapath, species)
            if not os.path.isdir(speciespath):
                continue
            for structure in os.listdir(speciespath):
                if structure.endswith('xyz'):
                    prefix = os.path.splitext(structure)[0]
                    fname = os.path.join(
                        minimapath, self.facetpath + '_' + species + '_'
                        + prefix + '_relax.py')
                    with open(fname, 'w') as f:
                        f.write(pytemplate.format(
                            adsorbate=species,
                            prefix=prefix,
                            pseudopotentials=pseudopotentials,
                            pseudo_dir=pseudo_dir,
                            balsam_exe_settings=balsam_exe_settings,
                            calc_keywords=calc_keywords,
                            creation_dir=self.creation_dir,
                            repeats=self.repeats,
                            n_kpts=n_kpts))
                    if shtemplate is None:
                        continue
                    # optional
                    fname = os.path.join(
                        minimapath, f'{species}_{prefix}_relax.sh')
                    with open(fname, 'w') as f:
                        f.write(shtemplate.format(adsorbate=species,
                                                  prefix=prefix))
