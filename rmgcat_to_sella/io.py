import os
import shutil
import yaml
import networkx as nx
from pathlib import Path

from catkit import Gratoms
from catkit.gen.molecules import get_3D_positions

from rmgcat_to_sella.graph_utils import node_test

from ase.io import read, write


class IO():
    """ Class for handling Input/Output and transforming it to more usefull
        format for the rmgcat_to_sella """

    def get_facetpaths(self, symbol, surface_types):
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

    def open_yaml_file(self, yamlfile):
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

    def get_all_species(self, yamlfile):
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

    def prepare_react_list(self, rxn):
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
                images.append(get_3D_positions(species1))
                unique_species.append(species1)
                unique_bonds.append(bond)

        r_name_list = [str(species.symbols) for species in reactants]
        p_name_list = [str(species.symbols) for species in products]

        return r_name_list, p_name_list, images

    def get_rxn_name(self, rxn):
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

    def get_list_all_rxns_names(self, yamlfile):
        ''' Get a list with all reactions names '''

        # open .yaml file
        reactions = self.open_yaml_file(yamlfile)

        all_rxns = []
        for rxn in reactions:
            rxn_name = self.get_rxn_name(rxn)
            all_rxns.append(rxn_name)
        return all_rxns

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
            indices = list(subgraph.nodes)
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
                            'At most two bonds to the metal are allowed '
                            'per adsorbate!'
                        )
                    tags[i] = abs(tags[i])
            new_gratoms.set_tags(tags)
            bonds.append(bond)
            gratoms_list.append(new_gratoms)

        return gratoms_list, bonds

    def get_xyz_from_traj(self, path_to_minima, species):
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
        species_path = os.path.join(path_to_minima, species)
        for traj in sorted(os.listdir(species_path), key=str):
            if traj.endswith('.traj'):
                src_traj_path = os.path.join(species_path, traj)
                des_traj_path = os.path.join(
                    species_path, traj[:-5] + '_final.xyz')
                write(des_traj_path, read(src_traj_path))

    def depends_on(self, facetpath, yamlfile):
        ''' Returns a dictionary of adsorbate + surface calculations
        (step 01; .py files) that has to be finished before starting step 02
        for a particular reaction

        Returns:
        ________

        dependancy_dict : [str:list(str)]
            a dictionary with keys being reaction names and values are lists
            of .py files for step 01 that have to be finished to start 02 step
            for a given reaction
            e.g.

        '''
        path_to_minima = os.path.join(facetpath, 'minima')

        # get reactions from. .yaml file
        reactions = self.open_yaml_file(yamlfile)

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
                # ?? is needed to prevent returning e.g.
                # OH_00_relax.py and O_00_relax.py
                # while
                # reactant or product = O
                # I have no idea why OH and HO is getting reverse
                # a workaround
                if reactant == 'OH':
                    reactant = 'HO'
                lookup_phrase = reactant + '_??_relax.py'
                # find matching reatants
                minima_py_files = Path(path_to_minima).glob(lookup_phrase)
                # append a list with minima that have to be calculated to
                # run 02 step
                for minima_py_file in minima_py_files:
                    minima_py_list.append(
                        os.path.split((str(minima_py_file)))[1])
            # loop through all products and do the same as for reactants
            for product in p_name_list:
                lookup_phrase = product + '_??_relax.py'
                minima_py_files = Path(path_to_minima).glob(lookup_phrase)
                for minima_py_file in minima_py_files:
                    minima_py_list.append(
                        os.path.split((str(minima_py_file)))[1])

            # create a dictionary with dependencies
            # {'reaction_name':[list_with_py_files_have_to_be_calculated]}
            dependancy_dict[rxn_name] = minima_py_list
        return dependancy_dict

    def clean_finished_subjobs(self):
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
