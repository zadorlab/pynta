#!/usr/bin/env python
import os
from pathlib import PosixPath
import yaml
from typing import Tuple, List, Dict

import numpy as np

from ase.io import read, write
from ase.data import covalent_radii

from pynta.excatkit.molecule import Molecule
from pynta.excatkit.adsorption import Builder
from pynta.excatkit.gratoms import Gratoms
from pynta.graph_utils import node_test
from pynta.io import IO

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

    def adjacency_to_3d(self) -> None:
        ''' Place adsorbates on the surface '''

        all_species_symbols = IO.get_all_unique_species(self.yamlfile)
        print(all_species_symbols)

        images = []
        # convert str to Gratoms and deal with edge cases
        for species_symbol in all_species_symbols:
            index = 0
            if species_symbol == 'COOH':
                index = 1
            images.append(Molecule().molecule(species_symbol)[index])

        slabedges, tags = Adsorbates.get_edges(self, True)

        # convert .xyz file to Atoms object
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

        structures = dict()
        for sp_symbol, sp_gratoms in zip(all_species_symbols, images):
            if len(sp_gratoms) == 0:
                continue
            try:
                # put adsorbates on the surface
                structs = ads_builder.add_adsorbate(
                    sp_gratoms, index=-1, bonds=[0])
                structures[str(sp_symbol)] = structs
            except IndexError:
                # TODO an idea to put -1 in adsorbate.get_tags() for
                # surface bonded atom
                print(sp_gratoms, sp_gratoms.edges,
                      sp_gratoms.get_tags())

        for species_symbol, adsorbate in structures.items():
            # create dir
            savedir = os.path.join(
                self.creation_dir, self.facetpath, 'minima', species_symbol)
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
        pytemplate : python file
            a template to prepare submission scripts
            for adsorbate+surface minimization
        pseudopotentials : dict(str: str)
            a dictionary with QE pseudopotentials for all species.
            e.g.
            dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                H='H.pbe-kjpaw_psl.1.0.0.UPF',
                O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                C='C.pbe-n-kjpaw_psl.1.0.0.UPF',
                )
        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            '/home/mgierad/espresso/pseudo'
        balsam_exe_settings : dict{str:int}
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.
            balsam_exe_settings = {'num_nodes': 1,
                                   'ranks_per_node': 48,
                                   'threads_per_rank': 1}
        calc_keywords : dict{str:str}
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}
        shtemplate : .sh file
            optional, a .sh template (not required by the workflow)

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
