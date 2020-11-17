import os
import shutil
import yaml
import networkx as nx
from pathlib import Path, PosixPath
from typing import List, Tuple, Optional, Dict
import numpy as np

from pynta.excatkit.gratoms import Gratoms
from pynta.graph_utils import node_test

from ase.io import read, write
from ase.dft.kpoints import monkhorst_pack
from ase.utils.structure_comparator import SymmetryEquivalenceCheck


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
            e.g. 'Cu'
        surface_type : str
            type of the surface, i.e. facet.
            e.g. 'fcc111'

        Returns
        -------
        facetpath : str
            a name of the facetpath,
            eg. 'Cu_111'

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
        ''' Generate a list with all facetpaths for a
        given surface defined by a symbol and a surface_type

        Parameters
        ----------
        symbol : str
            atomic symbol of the studied metal surface
            e.g. 'Cu'
        surface_types : list(str)
            a list with all surface types, i.e. facets.
            e.g. ['fcc111', 'fcc100']

        Returns
        -------
        facetpaths : list(str)
            a list with all facetpath names,
            e.g. ['Cu_111', 'Cu_100']

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

        Parameters:
        ___________
        size : tuple(int, int, int):
            a size or repeats of the slab,
            e.g. (3, 3, 1)
        get_uniq_kpts : bool, optional
            If True, return size and an ndarray of unique kpoints
            Otherwise False.

        Returns:
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
            a name of the .yaml file with a reaction list

        Returns
        -------
        species_dict
            a dictionary where keys are reactions (in a rxn{#} format)
            and values are species considered to moved in that reaction
            e.g.
            species_dict = {'rxn1': ['O', 'H'], 'rxn2': ['C', 'H']}

        '''
        species_dict = {}
        reactions = self.open_yaml_file(yamlfile)
        for num, rxn in enumerate(reactions):
            r_name_list, p_name_list, _ = self.prepare_react_list(rxn)
            if len(r_name_list) >= len(p_name_list):
                species_dict['rxn{}'.format(num)] = r_name_list
            else:
                species_dict['rxn{}'.format(num)] = p_name_list
        return species_dict

    @staticmethod
    def open_yaml_file(
            yamlfile: str) -> List[Dict[str, str]]:
        ''' Open yaml file with list of reactions

        Parameters:
        ___________
        yamlfile : str
            a name of the .yaml file with a reaction list

        Returns:
        __________
        reactions : list[dict{str:str}]
            a list with each reaction details stored as a dictionary

        '''
        with open(yamlfile, 'r') as f:
            yamltxt = f.read()
        reactions = yaml.safe_load(yamltxt)
        return reactions

    def get_all_species(
            self,
            yamlfile: str) -> List[str]:
        ''' Generate a list with all unique species for all reactions
            combined

        Parameters:
        ___________
        yamlfile : str
            a name of the .yaml file with a reaction list

        Returns:
        ________

        all_species_unique : list[str]
            a list with all unique species

        '''
        reactions = self.open_yaml_file(yamlfile)
        all_species = []
        for rxn in reactions:
            r_name_list, p_name_list, _ = self.prepare_react_list(rxn)
            all_species.append(r_name_list)
            all_species.append(p_name_list)
        all_species_unique = list(
            set([sp for sublist in all_species for sp in sublist]))
        return all_species_unique

    def get_all_species_given_rxn(
            self,
            rxn: Dict[str, str]) -> List[str]:
        ''' Get the reaction name

        Paremeters:
        ___________

        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file into a single reaction
            .yaml file

        Returns:
        _______
        rxn_name : str
            a name of the reaction in the following format:
            'OH_H+O'
        '''
        r_name_list, p_name_list, _ = self.prepare_react_list(rxn)
        all_species_rxn = r_name_list + p_name_list
        return all_species_rxn

    def prepare_react_list(
            self,
            rxn: Dict[str, str]) -> Tuple[List[str], List[str], List[Gratoms]]:
        '''Convert yaml file to more useful format

        Paremeters:
        ___________

        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file

        Returns:
        _______
        r_name_list : list(str)
            a list with all reactants for the given reaction
        p_name_list : list(str)
            a list with all products for the given reaction
        images : list(Gratoms)
            a list of CatKit's Gratom object (both reactants and products)

        '''

        species_ind = []
        bonds = []
        unique_species = []
        unique_bonds = []
        images = []

        # transforming reactions data to gratom objects
        reactants, rbonds = self.rmgcat_to_gratoms(
            rxn['reactant'].split('\n'))
        products, pbonds = self.rmgcat_to_gratoms(
            rxn['product'].split('\n'))
        species_ind += reactants + products
        bonds += rbonds + pbonds
        # check if any products are the same as anßy reactants
        for species1, bond in zip(species_ind, bonds):
            for species2 in unique_species:
                if nx.is_isomorphic(species1.graph, species2.graph, node_test):
                    break
            else:
                # images.append(Molecule().get_3D_positions(species1))
                images.append(species1)
                unique_species.append(species1)
                unique_bonds.append(bond)

        r_name_list = [str(species.symbols) for species in reactants]
        p_name_list = [str(species.symbols) for species in products]

        return r_name_list, p_name_list, images

    def get_rxn_name(
            self,
            rxn: Dict[str, str]) -> str:
        ''' Get the reaction name

        Paremeters:
        ___________

        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file into a single reaction
            .yaml file

        Returns:
        _______
        rxn_name : str
            a name of the reaction in the following format:
            'OH_H+O'
        '''
        r_name_list, p_name_list, _ = self.prepare_react_list(rxn)

        r_name = '+'.join([species for species in r_name_list])
        p_name = '+'.join([species for species in p_name_list])

        rxn_name = r_name + '_' + p_name
        return rxn_name

    def get_list_all_rxns_names(
            self,
            yamlfile):
        ''' Get a list with all reactions names '''

        # open .yaml file
        reactions = self.open_yaml_file(yamlfile)

        all_rxns = []
        for rxn in reactions:
            rxn_name = self.get_rxn_name(rxn)
            all_rxns.append(rxn_name)
        return all_rxns

    @ staticmethod
    def rmgcat_to_gratoms(
            adjtxt: str) -> Tuple[List[Gratoms], List[int]]:
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
            if 'multiplicity' in line:
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
                    edges.append((i, j - 1))
                    # TODO bug -> change to the below for MG made reactions
                    # edges.append((i, j - 1))
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
        for i, subgraph in enumerate(
            nx.connected_component_subgraphs(gratoms.graph)
        ):
            # TODO bug -> remove [::-1] to have a default order of atoms
            indices = list(subgraph.nodes)[::-1]
            symbols = gratoms[indices].symbols
            print(symbols)
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
                            'At most two bonds to the metal are allowed '
                            'per adsorbate!'
                        )
                    tags[i] = abs(tags[i])
            new_gratoms.set_tags(tags)
            bonds.append(bond)
            gratoms_list.append(new_gratoms)
        print(gratoms_list)
        return gratoms_list, bonds

    @staticmethod
    def get_xyz_from_traj(
            path_to_species: str) -> None:
        ''' Convert all ASE's traj files to .xyz files for a given species

        Parameters:
        ___________
        path_to_minima : str
            a path to minima
            e.g. 'Cu_111/minima'
        species : str
            a species symbol
            e.g. 'H' or 'CO'

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
        (step 01; .py files) that has to be finished before starting step 02
        for a particular reaction

        Parameters:
        ___________

        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        yamlfile : str
            a name of the .yaml file with a reaction list
        creation_dir : str
            a path to the main working directory

        Returns:
        ________

        dependancy_dict : [str:list(str)]
            a dictionary with keys being reaction names and values are lists
            of .py files for step 01 that have to be finished to start 02 step
            for a given reaction
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
            r_name_list, p_name_list, _ = self.prepare_react_list(rxn)
            # get reaction name
            rxn_name = self.get_rxn_name(rxn)
            minima_py_list = []
            # loop through all reactants
            for reactant in r_name_list:
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
            for product in p_name_list:
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
        ''' Move finished subjob files to finised_tmp_scripts directory '''
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
            eg. 'Cu_111'
        yamlfile : str
            a name of the .yaml file with a reaction list

        Returns
        -------
        unique_adsorbates_prefixes: Dict[str, List[str]]
            a dictionary with a list with prefixes of symmetry distinct
            conformers for a given adsorbate

        '''
        unique_adsorbates_prefixes = {}
        path_to_minima = os.path.join(creation_dir, facetpath, 'minima')
        all_species = self.get_all_species(yamlfile)
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
            e.g. 'Cu_111/minima/CO'

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
