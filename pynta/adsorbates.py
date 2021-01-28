#!/usr/bin/env python
import os
from pathlib import PosixPath
from typing import Tuple, List, Dict

import numpy as np

from ase.io import read, write
from ase.data import covalent_radii

from pynta.excatkit.adsorption import Builder
from pynta.excatkit.gratoms import Gratoms
from pynta.io import IO


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

        Parameters
        ___________
        facetpath : str
            a path to the workflow's main dir
            e.g. ``'Cu_111'``
        slab : str
            a '.xyz' file name with the optimized slab
            e.g.
            ``'Cu_111_slab_opt.xyz'``
        repeats: tuple
            specify reapeats in (x, y, z) direction,
            eg.

            >>> repeats = (3, 3, 1)
        yamlfile : str
            a name of the .yaml file with reaction list
        creation_dir : posix
            a posix path to the working directory

        '''
        self.facetpath = facetpath
        self.slab = slab
        self.repeats = repeats
        self.creation_dir = creation_dir
        self.yamlfile = os.path.join(self.creation_dir, yamlfile)
        self.slab_atom = read(os.path.join(self.creation_dir, self.slab))
        self.big_slab = self.slab_atom * self.repeats
        self.nslab = len(self.slab_atom)

    def get_edges(
            self,
            find_surface: bool = False) -> Tuple[List[str], np.ndarray]:
        ''' Get adsorption edges

        Parameters
        ___________
        find_surface : bool
            specify whether to include surface or not
            default = False

        Returns
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

        # pairvecs = np.ndarray(0)
        pairvecs = np.zeros_like(slab_atom.positions)

        # if find_surface:
        #     pairvecs = np.zeros_like(slab_atom.positions)
        # else:
        #     pairvecs = np.ndarray(0)

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

    def get_grslab(self) -> Gratoms:
        ''' Convert surface slab Atoms object into Gratoms object

        Returns
        -------
        grslab : Gratoms
            Gratoms representation of the surface slab - ready to place
            adsorbates

        '''
        slabedges, tags = Adsorbates.get_edges(self, True)

        grslab = Gratoms(numbers=self.slab_atom.numbers,
                         positions=self.slab_atom.positions,
                         cell=self.slab_atom.cell,
                         pbc=self.slab_atom.pbc,
                         edges=slabedges)
        grslab.arrays['surface_atoms'] = tags

        return grslab

    def adjacency_to_3d(self) -> None:
        ''' Place adsorbates on the surface

        .. todo:: Add support for a bidentate adsorption

        '''
        all_species_symbols = IO.get_all_unique_species_symbols(self.yamlfile)
        all_images_with_bonds = IO.get_all_unique_images_with_bonds(
            self.yamlfile)

        # prepare surface for placing adsorbates
        grslab = self.get_grslab()
        ads_builder = Builder(grslab)

        # build adsorbates
        structures = dict()
        for sp_symbol, unique_images_with_bonds in \
                zip(all_species_symbols, all_images_with_bonds.values()):
            for bond, sp_gratoms in unique_images_with_bonds.items():

                if len(sp_gratoms) == 0:
                    continue

                # which atom connects to the surface
                bonded = [bond]

                try:
                    # put adsorbates on the surface
                    structs = ads_builder.add_adsorbate(
                        sp_gratoms, index=-1, bonds=bonded)
                    structures[str(sp_symbol)] = structs
                except IndexError:
                    print('sp_gratoms, sp_gratoms.edges, sp.gratoms.tags')
                    print(sp_gratoms, sp_gratoms.edges,
                          sp_gratoms.get_tags())

        for sp_symbol, adsorbate in structures.items():
            # create directory where all adsorbates are stored
            savedir = os.path.join(
                self.creation_dir, self.facetpath, 'minima', sp_symbol)
            os.makedirs(savedir, exist_ok=True)

            for prefix, structure in enumerate(adsorbate):
                big_slab_ads = self.big_slab + structure[self.nslab:]
                # save adsorbates as .xyz and .png
                write(os.path.join(savedir, '{}.xyz'.format(
                    str(prefix).zfill(2))), big_slab_ads)
                write(os.path.join(savedir, '{}.png'.format(
                    str(prefix).zfill(2))), big_slab_ads)

    def create_relax_jobs(
            self,
            socket_calculator: str,
            pytemplate: str,
            pseudopotentials: Dict[str, str],
            pseudo_dir: str,
            balsam_exe_settings: Dict[str, int],
            calc_keywords: Dict[str, str],
            shtemplate: str = None) -> None:
        ''' Create a submit scripts

        Parameters
        __________
        pytemplate: str
            a template to prepare submission scripts
            for adsorbate+surface minimization
        pseudopotentials: Dict[str, str]
            a dictionary with QE pseudopotentials for all species.
            e.g.

            >>> dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                    H='H.pbe-kjpaw_psl.1.0.0.UPF',
                    O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                    C='C.pbe-n-kjpaw_psl.1.0.0.UPF')

        pseudo_dir: str
            a path to the QE's pseudopotentials main directory
            e.g.
            ``'/home/mgierad/espresso/pseudo'``
        balsam_exe_settings: Dict[str, int]
            a dictionary with balsam execute parameters(cores, nodes, etc.),
            e.g.

            >>> balsam_exe_settings={'num_nodes': 1,
                'ranks_per_node': 48,
                'threads_per_rank': 1}

        calc_keywords: Dict[str, str]
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            >>> calc_keywords={'kpts': (3, 3, 1),
                    'occupations': 'smearing',
                    'smearing': 'marzari-vanderbilt',
                    'degauss': 0.01,
                    'ecutwfc': 40,
                    'nosym': True,
                    'conv_thr': 1e-11,
                    'mixing_mode': 'local-TF'}

        shtemplate: str
            optional, a path to :literal:`*.sh` template(not required by the
            workflow but possible to specified for special cases)

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
                            socket_calculator=socket_calculator,
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
