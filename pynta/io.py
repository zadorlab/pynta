import os
import shutil
import yaml
from pathlib import Path, PosixPath
from typing import List, Tuple, Optional, Dict

import numpy as np
import networkx as nx

from pynta.excatkit.molecule import Molecule
from pynta.excatkit.gratoms import Gratoms
from pynta.graph_utils import node_test

from ase.io import read, write
from ase.dft.kpoints import monkhorst_pack
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from ase.io.formats import UnknownFileTypeError


class IO():
    ''' Class for handling Input/Output and transforming it to more usefull
    format for the pynta '''

    @staticmethod
    def get_facetpath(
            symbol: str,
            surface_type: str) -> None:
        ''' Get a facetpath for a given surface defined by a
        symbol and a surface_type

        Parameters
        ----------
        symbol : str
            atomic symbol of the studied metal surface
            e.g. ``'Cu'``
        surface_type : str
            type of the surface, i.e. facet.
            e.g. ``'fcc111'``

        Returns
        -------
        facetpath : str
            a name of the facetpath,
            eg. ``'Cu_111'``

        '''
        nums = []
        for num in surface_type:
            try:
                int(num)
            except ValueError:
                continue
            nums.append(num)
        facet = ''.join(nums)
        facetpath = symbol + '_' + facet
        return facetpath

    @staticmethod
    def get_facetpaths(
            symbol: str,
            surface_types: List[str]) -> List[str]:
        ''' Generate a list with all facetpaths for a given surface defined by
        a symbol and a surface_type

        Parameters
        ----------
        symbol : str
            atomic symbol of the studied metal surface
            e.g. ``'Cu'``
        surface_types : List[str]
            a list with all surface types, i.e. facets.
            e.g.

            >>> surface_types = ['fcc111', 'fcc100']

        Returns
        -------
        facetpaths : List[str]
            a list with all facetpath names,
            e.g.

            >>> facetpaths = ['Cu_111', 'Cu_100']

        '''
        facetpaths = []
        for stype in surface_types:
            nums = []
            for num in stype:
                try:
                    int(num)
                except ValueError:
                    continue
                nums.append(num)
            facet = ''.join(nums)
            facetpath = symbol + '_' + facet
            facetpaths.append(facetpath)
        return facetpaths

    @staticmethod
    def get_kpoints(
            size: Tuple[int, int, int],
            get_uniq_kpts: bool = False) -> Tuple[int, Optional[np.ndarray]]:
        ''' Returns number of unique k-points for a given size of the slab

        Parameters
        ----------
        size : Tuple[int, int, int]
            a size or repeats of the slab,
            e.g.

            >>> size = (3, 3, 1)

        get_uniq_kpts : bool, optional
            If True, return size and an ndarray of unique kpoints
            Otherwise False. By default False

        Returns
        -------
        m_uniq_kpts : int
            a number of unique k-points
        uniq : ndarray
            an array with unique k-points, optional

        '''
        kpts = monkhorst_pack(size)
        half_kpts = len(kpts) // 2
        uniq = kpts[half_kpts:, ]
        m_uniq_kpts = len(uniq)
        return (m_uniq_kpts, uniq) if get_uniq_kpts else m_uniq_kpts

    def get_species_dict(
            self,
            yamlfile: str) -> Dict[str, List[str]]:
        ''' For a given reaction get a dictionary with all species that takes
        part in the reaction.

        Those species will be considered as a reacting species by the
        TS esitmate constructor

        Parameters
        ----------
        yamlfile : str
            a name of the :literal:`*.yaml` file with a reaction list

        Returns
        -------
        species_dict : Dict[str, List[str]]
            a dictionary where keys are reactions (in a rxn{#} format)
            and values are species considered to moved in that reaction
            e.g.

            >>> species_dict = {'rxn1': ['O', 'H'], 'rxn2': ['C', 'H']}

        '''
        species_dict = {}
        reactions = self.open_yaml_file(yamlfile)
        for num, rxn in enumerate(reactions):
            r_name_list, p_name_list = IO.get_reactants_and_products(
                rxn)
            if len(r_name_list) >= len(p_name_list):
                species_dict['rxn{}'.format(num)] = r_name_list
            else:
                species_dict['rxn{}'.format(num)] = p_name_list
        return species_dict

    @staticmethod
    def open_yaml_file(
            yamlfile: str) -> List[Dict[str, str]]:
        ''' Open yaml file with list of reactions

        Parameters
        ----------
        yamlfile : str
            a name of the .yaml file with a reaction list

        Returns
        -------__
        reactions : List[Dict[str,str]]
            a list with each reaction details stored as a dictionary

        '''
        with open(yamlfile, 'r') as f:
            yamltxt = f.read()
        reactions = yaml.safe_load(yamltxt)
        return reactions

    @staticmethod
    def get_all_unique_species_symbols(
            yamlfile: str) -> List[str]:
        ''' Generate a list with all unique species names

        Parameters
        ----------
        yamlfile : str
            a name of the :literal:`*.yaml` file with a reaction list

        Returns
        -------
        all_species : List[str]
            a list with all unique species names in all reactions

        '''
        reactions = IO().open_yaml_file(yamlfile)

        all_sp_tmp = []
        for rxn in reactions:
            reactants_rxn, products_rxn = IO.get_reactants_and_products(
                rxn)
            all_sp_tmp.append(reactants_rxn)
            all_sp_tmp.append(products_rxn)

        all_species = [
            species for sublist in all_sp_tmp for species in sublist]

        # remove empty elements, such as ''
        all_species = [species for species in all_species if species]

        # remove duplicates keeping order - dictionary would do the job
        all_unique_species = {species: True for species in all_species}

        return list(all_unique_species.keys())

    @staticmethod
    def get_reactants_and_products(
            rxn: Dict[str, str]) -> Tuple[List[str], List[str]]:
        ''' For a given rxn, get lists with all reactants and products

        Parameters
        ----------
        rxn : Dict[str,str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction :literal:`*.yaml` file to a
            single reaction :literal:`*.yaml` file

        Returns
        -------
        Tuple[List[str], List[str]]
            lists with species names for reactants and products

        '''
        raw_rxn_name = rxn['reaction']
        raw_products, raw_reactants = raw_rxn_name.split('<=>')

        tmp_reactants = [react[:react.find('*')].strip()
                         if '*' in react else
                         react[:react.find('(')].strip()
                         for react in raw_reactants.split('+')]

        tmp_products = [prod[:prod.find('*')].strip()
                        if '*' in prod else
                        prod[:prod.find('(')].strip()
                        for prod in raw_products.split('+')]

        # remove all 'X' species
        r_x_removed = list(filter(('X').__ne__, tmp_reactants))
        p_x_removed = list(filter(('X').__ne__, tmp_products))

        # remove all empty '' elements
        reactants = [react for react in p_x_removed if react]
        products = [prod for prod in r_x_removed if prod]
        return reactants, products

    @staticmethod
    def get_rxn_name(
            rxn: Dict[str, str]) -> str:
        ''' Get a reaction name for a given rxn

        Parameters
        ----------
        rxn : Dict[str,str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction :literal:`*.yaml` file to a
            single reaction :literal:`*.yaml` file

        Returns
        -------
        rxn_name: str
            a reaction name,
            e.g. ``'CO_C+O'``

        '''
        reactants, products = IO.get_reactants_and_products(rxn)
        rxn_name = '+'.join(reactants) + '_' + '+'.join(products)
        return rxn_name

    @staticmethod
    def get_xyz_from_traj(
            path_to_species: str) -> None:
        ''' Convert all ASE's traj files to .xyz files for a given species

        Parameters
        ----------
        path_to_minima : str
            a path to minima
            e.g. ``'Cu_111/minima'``
        species : str
            a species symbol
            e.g. ``'H'`` or ``'CO'``

        '''
        for traj in sorted(os.listdir(path_to_species), key=str):
            if traj.endswith('.traj'):
                src_traj_path = os.path.join(path_to_species, traj)
                des_traj_path = os.path.join(
                    path_to_species, traj[:-5] + '_final.xyz')
                write(des_traj_path, read(src_traj_path))

    def depends_on(
            self,
            facetpath: str,
            yamlfile: str,
            creation_dir: PosixPath) -> Dict[str, List[str]]:
        ''' Returns a dictionary of adsorbate + surface calculations
        (step 01; :literal:`*.py` files) that has to be finished before
        starting step 02 for a particular reaction

        Parameters
        ----------

        facetpath : str
            a path to the workflow's main dir
            e.g. ``'Cu_111'``
        yamlfile : str
            a name of the .yaml file with a reaction list
        creation_dir : str
            a path to the main working directory

        Returns
        -------

        dependancy_dict : [str:List[str]]
            a dictionary with keys being reaction names and values are lists
            of :literal:`*.py` files for step 01 that have to be finished to
            start 02 step for a given reaction
            e.g.

        '''
        path_to_minima = os.path.join(creation_dir, facetpath, 'minima')
        path_to_yamlfile = os.path.join(creation_dir, yamlfile)

        # get reactions from. .yaml file
        reactions = self.open_yaml_file(path_to_yamlfile)

        dependancy_dict = {}

        # loop through all reactions
        for rxn in reactions:
            # get list of reactant and product
            reactants, products = IO.get_reactants_and_products(
                rxn)
            # get reaction name
            rxn_name = IO.get_rxn_name(rxn)
            minima_py_list = []
            # loop through all reactants
            for reactant in reactants:
                # I have no idea why OH and HO is getting reverse
                # a workaround
                lookup_phrase = '{}_{}_*relax.py'.format(facetpath, reactant)
                # find matching reatants
                minima_py_files = Path(path_to_minima).glob(lookup_phrase)
                # append a list with minima that have to be calculated during
                # run_02 step
                for minima_py_file in minima_py_files:
                    minima_py_list.append(
                        os.path.split((str(minima_py_file)))[1])
            # loop through all products and do the same as for reactants
            for product in products:
                lookup_phrase = '{}_{}_*relax.py'.format(facetpath, product)
                minima_py_files = Path(path_to_minima).glob(lookup_phrase)
                for minima_py_file in minima_py_files:
                    minima_py_list.append(
                        os.path.split((str(minima_py_file)))[1])

            # create a dictionary with dependencies
            dependancy_dict[rxn_name] = minima_py_list
        return dependancy_dict

    @staticmethod
    def clean_finished_subjobs() -> None:
        ''' Move finished subjob files to finished_tmp_scripts directory '''
        dir_name = 'finished_tmp_scripts'
        os.makedirs(dir_name, exist_ok=True)
        for prefix in range(0, 6):
            prefix = str(prefix).zfill(2)
            keyphrase = prefix + '*out'
            files = Path(os.getcwd()).glob(keyphrase)
            for file in files:
                file = str(file)
                if os.path.getsize(file) != 0:
                    # move all not empty .out files
                    shutil.move(file, dir_name)
                    # and corresponding .py.out files
                    shutil.move(file[:-4], dir_name)

    def get_unique_adsorbates_prefixes(
            self,
            facetpath: str,
            yamlfile: str,
            creation_dir: PosixPath) -> Dict[str, List[str]]:
        ''' Get a dictionary with a list with prefixes of symmetry distinct
        conformers for a given adsorbate

        Parameters
        ----------
        facetpath : str
            a name of the facetpath,
            eg. ``'Cu_111'``
        yamlfile : str
            a name of the :literal:`*.yaml` file with a reaction list

        Returns
        -------
        unique_adsorbates_prefixes: Dict[str, List[str]]
            a dictionary with a list with prefixes of symmetry distinct
            conformers for a given adsorbate

        '''
        unique_adsorbates_prefixes = {}
        path_to_minima = os.path.join(creation_dir, facetpath, 'minima')
        all_species = IO.get_all_unique_species_symbols(yamlfile)
        for species in all_species:
            path_to_species = os.path.join(path_to_minima, species)
            uq_prefixes = IO.get_unique_prefixes(
                path_to_species)
            unique_adsorbates_prefixes[species] = uq_prefixes
        return unique_adsorbates_prefixes

    @staticmethod
    def get_unique_prefixes(
            path_to_species: str) -> List[str]:
        ''' Compare each conformers for a given adsorbate and returns a list
        with prefixes of a symmetrty dictinct structures

        Parameters
        ----------
        path_to_species : str
            a path to species
            e.g. ``'Cu_111/minima/CO'``

        Returns
        -------
        unique_minima_prefixes : List[str]
            a list with prefixes of symmetry distinct structures for a given
            adsorbate

        '''
        good_minima = []
        result_dict = {}
        unique_minima_prefixes = []
        trajlist = sorted(Path(path_to_species).glob('*traj'), key=str)
        for traj in trajlist:
            minima = read(traj)
            comparator = SymmetryEquivalenceCheck(to_primitive=True)
            result = comparator.compare(minima, good_minima)
            result_dict[str(os.path.basename(traj).split('.')[0])] = result
            if result is False:
                good_minima.append(minima)
        for prefix, result in result_dict.items():
            if result is False:
                unique_minima_prefixes.append(prefix)
        return unique_minima_prefixes

    @staticmethod
    def get_unique_final_ts_prefixes(
            path_to_ts: str) -> List[str]:
        ''' Compare each conformers for a given adsorbate and returns a list
        with prefixes of a symmetrty dictinct structures

        Parameters
        ----------
        path_to_species : str
            a path to species
            e.g. ``'Cu_111/minima/CO'``

        Returns
        -------
        unique_minima_prefixes : List[str]
            a list with prefixes of symmetry distinct structures for a given
            adsorbate

        '''
        unique_ts = []
        result_dict = {}
        unique_ts_prefixes = []
        trajlist = sorted(Path(path_to_ts).glob('**/*traj'), key=str)
        for traj in trajlist:
            try:
                ts = read(traj)
                comparator = SymmetryEquivalenceCheck(to_primitive=True)
                result = comparator.compare(ts, unique_ts)
                result_dict[str(os.path.basename(traj)[:2])] = result
                if result is False:
                    unique_ts.append(ts)
            # error handling while calculations are still running/unfinished
            except UnknownFileTypeError:
                pass
        for prefix, result in result_dict.items():
            if result is False:
                unique_ts_prefixes.append(prefix)
        return unique_ts_prefixes

    def get_list_all_rxns_names(
            self,
            yamlfile: str) -> List[str]:
        ''' Get a list with all reactions names

        Parameters
        ----------
        yamlfile : str
            a name of the :literal:`*.yaml` file with a reaction list

        Returns
        -------
        all_rxns : List[str]
            a list with all reactions names

        '''
        # open .yaml file
        reactions = self.open_yaml_file(yamlfile)

        all_rxns = []
        for rxn in reactions:
            rxn_name = self.get_rxn_name(rxn)
            all_rxns.append(rxn_name)
        return all_rxns

    @staticmethod
    def get_all_reacting_atoms(yamlfile: str) -> Dict[str, Dict[str, float]]:
        ''' Read a :literal:`*.yaml` file with all reactions and extract
        reacting atoms symbols and indicies - the one with asterisk in the
        :literal:`*.yaml` file

        Parameters
        ----------
        yamlfile : str
            a name of the :literal:`*.yaml` file with a reaction list

        Returns
        -------
        Dict[str, Dict[str, float]]
            a dictionary with keys are 'rxn_1, rxn_2...'' and values are
            another dicts, where keys are atomic symbols, values are indicies,
            as they appear in the :literal:`*.yaml` file,

            e.g. ``'OH'`` for ``'H + O -> OH'``

            >>> all_reacting_atoms = {'O': 0, 'H': 1}

        '''
        all_reacting_atoms = {}
        reactions = IO.open_yaml_file(yamlfile)

        for num, rxn in enumerate(reactions):
            r_name_list, p_name_list = IO.get_reactants_and_products(rxn)
            if len(r_name_list) <= len(p_name_list):
                easier_to_build = 'reactant'
            else:
                easier_to_build = 'product'

            reacting_species_connectivity = rxn[easier_to_build].strip().split(
                '\n')

            reacting_atoms_and_idxs = IO.get_reacting_atoms_idx_dict(
                reacting_species_connectivity)
            all_reacting_atoms['rxn_{}'.format(num)] = reacting_atoms_and_idxs
        return all_reacting_atoms

    @staticmethod
    def get_reacting_atoms_idx_dict(
            reacting_species_connectivity: List[str]) -> Dict[str, int]:
        ''' Get a dict with atomic symbols and indicies of reacting atoms -
        the one with asterisk in the :literal:`*.yaml` file

        Parameters
        ----------
        reacting_species_connectivity : List[str]
            conectivity info for the species that is easier to use as a
            TS guess skeleton

            e.g. ``'OH'`` for ``'H + O -> OH'``

            >>> reacting_species_connectivity = ['multiplicity -187',
                                                '1 *1 O u0 p0 c0 {2,S} {4,S}',
                                                '2 *2 H u0 p0 c0 {1,S}',
                                                '3 *3 X u0 p0 c0',
                                                '4    X u0 p0 c0 {1,S}',
                                                '']

        Returns
        -------
        reacting_idxs = Dict[str, int]
            keys are atomic symbols, values are indicies, as they appear in the
            :literal:`*.yaml` file,
            e.g. ``'OH'`` for ``'H + O -> OH'``

            >>> reacting_idxs = {'O': 0, 'H': 1}

        '''
        reacting_idxs = {}
        n_surf_at_befor_ads = 0
        remove_one_more = 0
        for line in reacting_species_connectivity:
            if 'multiplicity' in line:
                continue
            if 'X' in line:
                n_surf_at_befor_ads += 1
            else:
                break
        for num, line in enumerate(reacting_species_connectivity):
            if 'multiplicity' in line:
                remove_one_more = 1
            if '*' in line and 'X' not in line:
                atom_symbol = line.split()[2]
                reacting_idxs[atom_symbol] = (
                    num - n_surf_at_befor_ads - remove_one_more)
        return reacting_idxs

    @staticmethod
    def get_TS_guess_image(
            rxn: Dict[str, str],
            easier_to_build: str) -> Gratoms:
        ''' Convert RMGCat representation of species to Gratom object
        - the case of TS_guess

        .. todo:: There is only one element for every reaction tested. |br| There will be a problem for AX + BX -> CX + DX

        Parameters
        ----------
        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            :literal:`*.yaml` file
        easier_to_build : str
            a list of species that are considerd as the reactiong one

        Returns
        -------
        ts_guess_image : Gratoms
            Skeleton of an TS that is used to generate TS_guesses in
            :meth:`pynta.general_ts_guesses.decide`

        '''
        ts_guess = IO.rmgcat_to_gratoms(
            rxn[easier_to_build].strip().split('\n'))

        ts_guess_image = Molecule().get_3D_positions(ts_guess[0])
        return ts_guess_image

    @staticmethod
    def get_all_unique_images_with_bonds(
            yamlfile: str) -> Dict[int, Dict[int, Gratoms]]:
        ''' Return a list with all unique images (Gratoms) for all reactions

        Parameters
        ----------
        yamlfile : str a name of the :literal:`*.yaml` file with a reaction
            list

        Returns
        -------
        all_images_with_bonds : Dict[int, Dict[int, Gratoms]] a Dict with all
            unique Gratoms objects for all reactions together with information
            which atoms connects Gratoms object to the surface


        >>> all_images_with_bonds =
            {
             0: {0: Gratoms(symbols='OH', pbc=False, tags=...)},
             1: {0: Gratoms(symbols='O', pbc=False, tags=...)},
             2: {0: Gratoms(symbols='H', pbc=False, tags=...)},
             3: {0: Gratoms(symbols='OCH2OCH2', pbc=False, tags=...)},
             4: {0: Gratoms(symbols='OCH3', pbc=False, tags=...)},
             5: {1: Gratoms(symbols='OCH', pbc=False, tags=...)},
             6: {0: Gratoms(symbols='OHCH', pbc=False, tags=...)}
            }

        '''
        reactions = IO.open_yaml_file(yamlfile)

        all_images_unique = []
        all_images = {}
        for rxn in reactions:
            images = IO.get_images(rxn)
            all_images.update(images)

        # check if any species in one reaction is the same as any species
        # in the other reaction
        unique_bonds = []
        for bonded_dict in all_images.values():
            for bond, gratoms1 in bonded_dict.items():
                for gratoms2 in all_images_unique:
                    if nx.is_isomorphic(
                            gratoms1.graph, gratoms2.graph, node_test):
                        break
                else:
                    unique_bonds.append(bond)
                    all_images_unique.append(gratoms1)

        all_images_with_bonds = {}

        i = 0
        for uq_image, bond in zip(all_images_unique, unique_bonds):
            all_images_with_bonds[i] = {bond: uq_image}
            i += 1
        return all_images_with_bonds

    @staticmethod
    def get_images(rxn: Dict[str, str]) -> Dict[str, Gratoms]:
        ''' Convert RMGCat representation of species for a given rxn
        to a dict of Gratoms objects - both reactants and products combined.

        .. note:: Reactants or products can be symmetrically identical at this step.


        Parameters
        ----------
        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction :literal:`*.yaml` file to a
            single reaction :literal:`*.yaml` file

        Returns
        -------
        images : List[Gratoms]
            a dict with all unique Gratoms object representing each species
            taking part in reaction rxn

        '''
        species = []
        rxn_bonded_dict = {}
        for reag in ['reactant', 'product']:
            adjtxt_react = rxn[reag].strip().split('\n')
            reagent = IO.rmgcat_to_gratoms(adjtxt_react)
            reactant_bonded_idx = IO.get_surface_bonded(adjtxt_react)
            species += reagent
            bonded_dict_reactants = IO.get_bonded_dict(
                reagent, reactant_bonded_idx)
            rxn_bonded_dict.update(bonded_dict_reactants)

        unique_species = []
        images = {}
        # check if any products are the same as any reactants
        for species1, bond in zip(species, rxn_bonded_dict.values()):
            for species2 in unique_species:
                if nx.is_isomorphic(
                        species1.graph, species2.graph, node_test):
                    break
            else:
                image = Molecule().get_3D_positions(species1)
                images[str(species1.symbols)] = {bond: image}
                unique_species.append(species1)
        return images

    @staticmethod
    def get_bonded_dict(
            reactants: List[Gratoms],
            reactant_bonded_idx: List[int]) -> Dict[str, int]:
        ''' Create a dict with bonding information for a specifc reaction and
        specific type of reagents (reactants or products)

        Parameters
        ----------
        reactants : List[Gratoms]
            a list gratoms objects
        reactant_bonded_idx : List[int]
            a list with indicies of atoms that are connected to the surface

        Returns
        -------
        bonded_dict : Dict[str, int]
            a dictionary with keys being atomic symbols and values are indicies
            of atoms that connects to the surface. Indicies are folowing the
            order of atoms as in Gratoms objec

        '''
        bonded_dict = {}
        for r in reactants:
            if not reactant_bonded_idx:
                bonded_dict[str(r.symbols)] = 0
            for atom in r:
                if atom.tag in reactant_bonded_idx:
                    bonded_dict[str(r.symbols)] = atom.index
        return bonded_dict

    @staticmethod
    def get_surface_bonded(adjtxt) -> List[int]:
        ''' For each reactant in adjtxt from :literal:`*.yaml` file, get info
        which atoms are connected to surface

        Parameters
        ----------
        adjtxt : List[str]
            a list with a connectivity info for reactant or product
            as from the .yaml file

            e.g. for given reaction (reactant or product)

            In :literal:`*.yaml` file we have something like that::

                    multiplicity -187
                1 *1 C u0 p0 c0 {2,S} {4,S}
                2    O u0 p0 c0 {1,S}
                3 *2 H u0 p0 c0 {5,S}
                4 *3 X u0 p0 c0 {1,S}
                5 *4 X u0 p0 c0 {3,S}

            but we need here a list like that:

            >>> tmp_list = ['multiplicity -187', '1 *1 C u0 p0 c0 {2,S} {4,S}',
                            '2    O u0 p0 c0 {1,S}', '3 *2 H u0 p0 c0 {5,S}',
                            '4 *3 X u0 p0 c0 {1,S}', '5 *4 X u0 p0 c0 {3,S}',
                            '']

            So it can be simply converted using the following:

            >>> yamlfile = 'reactions.yaml'
            >>> with open(yamlfile, 'r') as f:
                    text = f.read()
            >>> reactions = yaml.safe_load(text)
            >>> for rxn in reactions:
                    adjtxt = rxn['reactant'].split('\\n')

        Returns
        -------
        idx_surface_bonded_atoms : List[int]
            a list with indicies of atoms connected to the surface. Each index
            refers to a position of a atom in :literal:`*.yaml` file

        '''
        symbols = []
        edges = []
        tags = []

        n_surf_at_befor_ads = 0

        start_idx = 1
        if 'multiplicity' in adjtxt[0]:
            start_idx -= 1

        for i, line in enumerate(adjtxt, start_idx):
            if i == 0:
                continue
            if 'X' in line and int(line.split()[0]) <= 2:
                n_surf_at_befor_ads += 1
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
        idx_surface_bonded_atoms = np.argwhere(
            gratoms.get_tags() < 0).flatten().tolist()

        idx_surface_bonded_atoms = [
            idx - n_surf_at_befor_ads for idx in idx_surface_bonded_atoms]
        return idx_surface_bonded_atoms

    @staticmethod
    def rmgcat_to_gratoms(
            adjtxt: List[str]) -> List[Gratoms]:
        ''' Convert a slice of :literal:`*.yaml` file to Catkit's Gratoms object

        Parameters
        ----------
        adjtxt : List[str]
            a list with a connectivity info for reactant or product
            as from the .yaml file

            e.g. for given reaction (reactant or product)

            In :literal:`*.yaml` file we have something like that::

                    multiplicity -187
                1 *1 C u0 p0 c0 {2,S} {4,S}
                2    O u0 p0 c0 {1,S}
                3 *2 H u0 p0 c0 {5,S}
                4 *3 X u0 p0 c0 {1,S}
                5 *4 X u0 p0 c0 {3,S}

            but we need here a list like that:

            >>> tmp_list = ['multiplicity -187', '1 *1 C u0 p0 c0 {2,S} {4,S}',
                            '2    O u0 p0 c0 {1,S}', '3 *2 H u0 p0 c0 {5,S}',
                            '4 *3 X u0 p0 c0 {1,S}', '5 *4 X u0 p0 c0 {3,S}',
                            '']

            So it can be simply converted using the following:

            >>> yamlfile = 'reactions.yaml'
            >>> with open(yamlfile, 'r') as f:
                    text = f.read()
            >>> reactions = yaml.safe_load(text)
            >>> for rxn in reactions:
                    adjtxt = rxn['reactant'].split('\\n')

        Returns
        -------
        gratoms_list : List[Gratoms]
            a list with all Gratoms objects for a give rxn

        '''
        symbols = []
        edges = []
        tags = []

        start_idx = 1
        if 'multiplicity' in adjtxt[0]:
            start_idx -= 1

        for i, line in enumerate(adjtxt, start_idx):
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
        graphs = [gratoms.graph.subgraph(
            c) for c in nx.connected_components(gratoms.graph)]
        for i, subgraph in enumerate(graphs):
            indices = list(subgraph.nodes)
            symbols = gratoms[indices].symbols
            new_indices = {old: new for new, old in enumerate(indices)}
            new_edges = []

            for edge in subgraph.edges:
                newa = new_indices[edge[0]]
                newb = new_indices[edge[1]]
                new_edges.append((newa, newb))

            new_gratoms = Gratoms(symbols, edges=new_edges)

            tags = indices
            new_gratoms.set_tags(tags)
            gratoms_list.append(new_gratoms)
        return gratoms_list

    @staticmethod
    def get_calculators(quantum_chemistry: str) -> Tuple[str, str]:
        ''' Get a proper names of quantum_chemistry calculators
        (socket and balsamsocket)

        Parameters
        ----------
        quantum_chemistry : str
            a keyword describing which quantum chemistry package to use

        Returns
        -------
        calculator, socket_calculator : Tuple[str, str]
            well formated calculator and socket_calculator that refers to
            spefic method in :mod: balsamcalc

        Raises
        ------
        Exception
            If not supported calculator was requested

        '''
        if quantum_chemistry == 'espresso':
            calculator = 'EspressoBalsam'
        elif quantum_chemistry == 'nwchem':
            calculator = 'NWChemBalsam'
        else:
            raise Exception('Pynta only supports the following '
                            'quantum chemistry packages (keywords): '
                            '   espresso'
                            '   nwchem')
        socket_calculator = calculator + 'SocketIO'

        return calculator, socket_calculator

    @staticmethod
    def set_calculators(
            executable: str,
            calculator: str,
            socket_calculator: str) -> None:
        ''' Create balsam socket and quantum chemistry calculator for requested
        quantum chemistry code

        Parameters
        ----------
        executable : str
            a path to system executable to python3

            >>> executable = '/Users/mgierad/opt/anaconda3/bin/python3'

        calculator : str
            name of quantum chemistry calculator
        socket_calculator : str
            name of quantum chemistry calculator bundle with balsam\

        '''
        balsamcalc_module = __import__('pynta.balsamcalc', fromlist=[
            socket_calculator])

        calc = getattr(balsamcalc_module, calculator)
        sock_calc = getattr(balsamcalc_module, socket_calculator)

        calc.exe = executable
        sock_calc.exe = executable
        calc.create_application()
        sock_calc.create_application()
