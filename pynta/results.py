from pynta.io import IO

import os
import numpy as np
import json
import sys
from ase.io import read
from pathlib import Path
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


class Results():
    def __init__(self) -> None:
        ''' A class to analyze the results of calculations '''

        input_file = 'input.json'
        try:
            with open(input_file, 'r') as f:
                input_json = json.load(f)

            surface_types_and_repeats = input_json['surface_types_and_repeats']
            metal_atom = input_json['metal_atom']
            self.yamlfile = input_json['yamlfile']
            self.facetpaths = IO().get_facetpaths(
                metal_atom, surface_types_and_repeats.keys())
            self.reactions = IO().open_yaml_file(self.yamlfile)
            self.ev_to_kjmol = 23.06035 * 4.184
            self.slab_paths = ['{}_big_slab_opt.xyz'.format(
                facetpath) for facetpath in self.facetpaths]
            self.current_dir = os.getcwd()

        except FileNotFoundError:
            print('!    input.json not found. \n'
                  '     Make sure {} matches exactly "input.json"\n'
                  '\n'
                  '     You can use this module only from your main working \n'
                  '     directory - the one with all your input files and \n'
                  '     "facetpath" dir and job_files dir.\n'
                  '\n'
                  '     Your current directory is:\n'
                  '     {}'.format(input_file, os.getcwd()))
            sys.exit()

    def get_reaction_energies_all(self) -> None:
        ''' Get reactiom energy (kj/mol) for all facetpaths and reactions

        Returns
        -------
        reaction_energies : Dict[str:float]
            a dictionary with reaction energies for a given facetpath and
            rxn_name, e.g.

            >>> reaction_energies = {'Cu_111_OH_O+H': 70.81}

        '''
        reaction_energies = {}
        for facetpath, slab_path in zip(self.facetpaths, self.slab_paths):
            minima_path = os.path.join(facetpath, 'minima')
            for rxn in self.reactions:
                r_name_list, p_name_list = IO.get_reactants_and_products(rxn)
                rxn_name = IO.get_rxn_name(rxn)
                key = facetpath + '_' + rxn_name
                reaction_energies[key] = float(self.get_reaction_energy(
                    minima_path, facetpath, r_name_list, p_name_list,
                    slab_path))
        return reaction_energies

    def get_reaction_energy(
            self,
            minima_path: str,
            facetpath: str,
            r_name_list: List[str],
            p_name_list: List[str],
            slab_path: str) -> str:
        ''' Calclate reaction energy as a difference between the most stable
        product and the most stable reactant

        Parameters
        ----------
        minima_path : str
            a path to main minima directory
            e.g.

            >>> minima_path - 'Cu_111/minima'

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        r_name_list : List[str]
            a list with all reactants for the given reaction
            e.g.

            >>> r_name_list = ['OH']

        p_name_list : List[str]
            a list with all products for the given reaction
            e.g.

            >>> p_name_list = ['O', 'H']

        slab_path : str
            a path to the slab
            e.g.

            >>> slab_path = 'Cu_111_slab_opt.xyz'

        Returns
        -------
        reaction_energy : str
            an energy of reaction calculated a difference between
            the most stable product and the most stable reactant

        '''
        r_ener_list, p_ener_list, slab_ener, nslabs = Results.get_data(
            minima_path, facetpath, r_name_list, p_name_list, slab_path)
        # Depending how the reactants and products are defined,
        # there are three options here:
        # e.g. AB --> A + B
        if len(p_ener_list) > len(r_ener_list):
            reaction_energy = (sum(p_ener_list) - slab_ener * nslabs -
                               r_ener_list) * self.ev_to_kjmol
        # e.g. A + B --> AB
        elif len(p_ener_list) < len(r_ener_list):
            reaction_energy = (p_ener_list + slab_ener * nslabs -
                               sum(r_ener_list)) * self.ev_to_kjmol
        # e.g. A + B --> C + D
        else:
            reaction_energy = (sum(p_ener_list) -
                               sum(r_ener_list)) * self.ev_to_kjmol
        reaction_energy = '{:.2f}'.format(round(reaction_energy[0], 3))
        return reaction_energy

    def get_barrier_all(self) -> Dict[str, Dict[str, str]]:
        ''' Get barrier heights for all rxn_names and facetpaths

        Returns
        -------
        ts_ener : Dict[str,Dict[str,str]]
            a dictionary with all barrier heights (kj/mol)
            in a format like below

            >>> ts_ener = {'Cu_111_OH_O+H':
                    {'TS_00': '155.27', 'TS_01': '157.97'}}

        '''
        ts_ener = {}
        for facetpath, slab_path in zip(self.facetpaths, self.slab_paths):
            minima_path = os.path.join(facetpath, 'minima')
            for rxn in self.reactions:
                r_name_list, p_name_list = IO.get_reactants_and_products(
                    rxn)
                rxn_name = IO.get_rxn_name(rxn)
                ts_path = os.path.join(self.current_dir,
                                       facetpath,
                                       rxn_name,
                                       'TS_estimate_unique')
                ts_ener[facetpath + '_' + rxn_name] = self.get_barrier(
                    minima_path, ts_path, facetpath, r_name_list, p_name_list,
                    slab_path)
        return ts_ener

    def get_barrier(
            self,
            minima_path: str,
            ts_path: str,
            facetpath: str,
            r_name_list: List[str],
            p_name_list: List[str],
            slab_path: str) -> Dict[str, float]:
        ''' Calculate reaction energy relatively to the most stable reactant

        Parameters
        ----------
        minima_path : str
            a path to main minima directory, e.g.

            >>> minima_path = 'Cu_111/minima'

        ts_path : str
            a path to main TS directory, e.g.

            >>> ts_path = 'Cu_111/TS_estimate_unique'

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        r_name_list : List[str]
            a list with all reactants for the given reaction
            e.g.

            >>> r_name_list = ['OH']

        p_name_list : List[str]
            a list with all products for the given reaction
            e.g.

            >>> p_name_list = ['O', 'H']

        slab_path : str
            a path to the slab
            e.g.

            >>> slab_path = 'Cu_111_slab_opt.xyz'


        Returns
        -------

        activation_barriers : Dict[str, float]
            a dictonary with all barrier heights relatively
            to the most stable reactant (in kJ/mol)

        '''

        r_ener_list, _, slab_ener, _ = Results.get_data(
            minima_path, facetpath, r_name_list, p_name_list, slab_path)
        tss_ener = Results.get_ts_ener(ts_path)
        tss_name = Results.format_TS_name(ts_path)

        activation_barriers = {}
        for ts_ener, ts_name in zip(tss_ener, tss_name):
            # Should be valud for any type of TS, no matter how many reactants
            # take part in the reaction, i.e. A + B -> AB or AB -> A + B
            barrier = (ts_ener + slab_ener * (len(r_ener_list) - 1) -
                       sum(r_ener_list)) * self.ev_to_kjmol
            activation_barriers['TS_' +
                                ts_name] = '{:.2f}'.format(barrier)
        return activation_barriers

    @staticmethod
    def get_data(
            minima_path: str,
            facetpath: str,
            r_name_list: List[str],
            p_name_list: List[str],
            slab_path: str) -> Tuple[List[float], List[float], float, int]:
        ''' Returns the lowest energies lists for reactants and products.

        Parameters
        ----------
        minima_path : str
            a path to main minima directory, e.g.

            >>> minima_path = 'Cu_111/minima'

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        r_name_list : List[str]
            a list with all reactants for the given reaction
            e.g.

            >>> r_name_list = ['OH']

        p_name_list : List[str]
            a list with all products for the given reaction
            e.g.

            >>> p_name_list = ['O', 'H']

        slab_path : str
            a path to the slab
            e.g.

            >>> slab_path = 'Cu_111_slab_opt.xyz'

        Returns
        -------
        r_ener_list : List[float]
            a list with the lowest energy conformer for each reactant
        p_ener_list : List[float]
            a list with the lowest energy conformer for each products
        slab_ener : float
            an energy of the representative slab
            of the size the same as reactants, TS, or product
        nslabs : int
            a number specifying how many additional slabs have to be
            considered in a stoichiometric reaction. Its defined as the
            absoulute valueof the difference between amount of products and
            reactants
            e.g. for
            O + H --> OH

            >>> nslabs = 1

            or, for C + O + H --> COH

            >>> nslabs = 2

        Raises
        ------
        TypeError
            If :literal:`*.out` files are missing either for reactants
            or products

        '''
        r_ener_list = []
        p_ener_list = []

        # get the lowest energy for all reactants
        for reactant in r_name_list:
            lowest_reactant_ener = Results.get_lowest_species_ener(
                minima_path, reactant, facetpath)
            r_ener_list.append(lowest_reactant_ener)
        # check if .out files for reactants are copied
        if None in r_ener_list:
            print('----')
            print(
                'Found None in reactants energy list. Missing .out files '
                'for reactants.')
            print('----')
            raise TypeError
        # get the lowest energy for all products
        for product in p_name_list:
            lowest_product_ener = Results.get_lowest_species_ener(
                minima_path, product, facetpath)
            p_ener_list.append(lowest_product_ener)
        # check if .out files for products are copied
        if None in p_ener_list:
            print('----')
            print('Found None in products energy list. Missing .out files'
                  'for products.')
            print('----')
            raise TypeError

        slab_ener = Results.get_slab_ener(slab_path)
        nslabs = abs(len(p_ener_list) - len(r_ener_list))

        return r_ener_list, p_ener_list, slab_ener, nslabs

    @staticmethod
    def get_slab_ener(
            slab_path: str) -> float:
        ''' Get energy of the slab

        Parameters
        ----------
        slab_path : str
            a path to the slab, e.g.

            >>> slab_path = 'Cu_111_slab_opt.xyz'

        Returns
        -------
        slab_ener : float
            an energy of the slab in eV

        '''
        slab = read(slab_path)
        slab_ener = slab.get_potential_energy()
        return slab_ener

    @staticmethod
    def get_ts_ener(
            ts_path: str) -> List[float]:
        ''' Get energy of all TSs

        Parameters
        ----------
        ts_path : str
            a path to main TS directory, e.g.

            >>> ts_path = 'Cu_111/TS_estimate_unique'

        Returns
        -------
        ts_ener_list : List[float]
            a sorted list with energy value for each TS

        '''
        ts_ener_dict = {}
        tss = Results.get_ts_out_files(ts_path)
        for ts in tss:
            with open(ts, 'r') as f:
                # error handiling while calculations are in progress and
                # .out file is empty or it is missing
                try:
                    data = f.readlines()
                    enerLine = data[-1]
                    enerVal = enerLine.split()
                    ts_ener_dict[ts] = float(enerVal[3])
                except IndexError:
                    pass
        ts_ener_list = list(ts_ener_dict.values())
        return ts_ener_list

    @staticmethod
    def get_lowest_species_ener(
            minima_path: str,
            species: str,
            facetpath: str) -> float:
        ''' Get the lowest energy of the most stable species

        Parameters
        ----------
        minima_path : str
            a path to main minima directory, e.g.

            >>> minima_path = 'Cu_111/minima'

        species : str
            a metal_atom of the species
            e.g.

            >>> species = 'OH'
            >>> species = 'H'
            >>> species = 'O'

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Returns
        -------
        lowest_species_ener: float
            conformer of the lowest energy among all
            calculated for the given species

        '''
        species_ener_dict = {}
        try:
            species_out_file_path_list = Results.get_species_out_files(
                minima_path, species, facetpath)
            for spiecies_out_file_path in species_out_file_path_list:
                with open(spiecies_out_file_path, 'r') as f:
                    data = f.readlines()
                    enerLine = data[-1]
                    enerVal = enerLine.split()
                    species_ener_dict[spiecies_out_file_path] = float(
                        enerVal[4])
            lowest_species_ener = min(species_ener_dict.values())
            return lowest_species_ener
        except ValueError:
            print('Minima .out files probably not copied successfully.')
            print('Check minima .out files. If missing, copy it.')

    @staticmethod
    def get_ts_out_files(
            ts_path: str) -> List[str]:
        ''' Get TS :literal:`*.out` files

        Parameters
        ts_path : str
            a path to main TS directory, e.g.

            >>> ts_path = 'Cu_111/TS_estimate_unique'

        Returns
        -------
        ts_out_file_list : List[str]
            a sorted list with paths to Sella's .out files for each TSs

        '''
        ts_out_file_list = []
        unique_ts_prefixes = IO.get_unique_final_ts_prefixes(ts_path)

        try:
            for outfile in os.listdir(ts_path):
                for uq_ts_prefix in unique_ts_prefixes:
                    if outfile.startswith(uq_ts_prefix) and outfile.endswith('out'):
                        uq_ts_outfile_path = os.path.join(ts_path, outfile)
                        ts_out_file_list.append(str(uq_ts_outfile_path))

        except FileNotFoundError:
            pass

        return sorted(ts_out_file_list)

    @staticmethod
    def get_species_out_files(
            minima_path: str,
            species: str,
            facetpath: str) -> List[str]:
        ''' Get .out files for each reactants

        Parameters
        ----------
        minima_path : str
            a path to main minima directory, e.g.

            >>> minima_path = 'Cu_111/minima'

        species : str
            a metal_atom of the species
            e.g.

            >>> species = 'OH'
            >>> species = 'H'
            >>> species = 'O'

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Returns
        -------
        species_out_file_path_list : List[str]
            a list with paths to all minima Sella's * out files for given
            species
            e.g.

            >>> species_out_file_path_list = ['Cu_111/minima/OH_01_relax.out',
                    'Cu_111/minima/OH_00_relax.out',
                    'Cu_111/minima/OH_03_relax.out',
                    'Cu_111/minima/OH_02_relax.out']

        '''
        species_out_file_path_list = []
        species = species + '_'
        outfile = '{}_{}*out'.format(facetpath, species)
        reactant_out_list = Path(minima_path).glob(outfile)
        for reactant_out_file in reactant_out_list:
            species_out_file_path_list.append(str(reactant_out_file))
        return species_out_file_path_list

    @staticmethod
    def format_TS_name(
            ts_path: str) -> List[str]:
        ''' Function to get prefixes of TSs

        Parameters
        ----------
        ts_path : str
            a path to main TS directory, e.g.

            >>> ts_path = 'Cu_111/TS_estimate_unique'

        Returns
        -------
        prefix_list : List[str]
            a list with all prefixes for TSs

        '''
        prefix_list = []
        ts_out_file_list = Results.get_ts_out_files(ts_path)
        for ts_out_file in ts_out_file_list:
            prefix = os.path.split(ts_out_file)[1].split('_')[0]
            prefix_list.append(prefix)
        return prefix_list

    @staticmethod
    def rxn_title(
            rxn_name: str) -> str:
        ''' Return rxn name with arrow between reactants and products

        Parameters
        ----------
        rxn_name : str
            a reaction formatted as:

            >>> rxn_name = 'OH_O+H'

        Returns
        ----------
        rxn_name_title : str
            a reaction name having a format

            >>> rxn_name_title = 'OH --> O+H'

        '''
        reactants, products = rxn_name.split('_')
        rxn_name_title = reactants + ' --> ' + products
        return rxn_name_title

    def plot(
            self,
            plot_filename: str = None,
            max_barrier: float = None,
            x_size_in: float = 8,
            y_size_in: float = 6) -> None:
        ''' Plot all results and automatically detect how many reactions and
        facet types exist.

        Parameters
        ----------
        plot_filename: str, optional
            a name of the file of generated plot,
            by default ``None``
        max_barrier: float, optional
            all barriers lower then max_barrier will be ploted,
            by default ``None``
        x_size_in: float, optional
            x size of the plot in inches, by default ``8``
        y_size_in: float, optional
            y size of the plot in inches, by default ``6``

        '''
        reaction_energies = self.get_reaction_energies_all()
        activation_barriers = self.get_barrier_all()
        all_rxn_names = IO().get_list_all_rxns_names(self.yamlfile)
        n_facets = len(self.facetpaths)
        n_rxns = len(all_rxn_names)
        _, axes = plt.subplots(n_facets, n_rxns, squeeze=False)

        # print(axes)
        for num, rxn in enumerate(self.reactions):
            for ax, facetpath, in zip(axes, self.facetpaths):
                rxn_name = IO().get_rxn_name(rxn)
                key = facetpath + '_' + rxn_name
                Results.plot_rxn(key,
                                 reaction_energies,
                                 activation_barriers,
                                 rxn_name,
                                 ax,
                                 num,
                                 plot_filename,
                                 max_barrier,
                                 x_size_in,
                                 y_size_in)

    @staticmethod
    def plot_rxn(
            key: str,
            reaction_energies: Dict[str, float],
            activation_barriers: Dict[str, Dict[str, str]],
            rxn_name: str,
            axes: plt.subplot,
            num: int,
            plot_filename: str = None,
            max_barrier: float = None,
            x_size_in: float = 8,
            y_size_in: float = 6) -> None:
        ''' Plot reaction energy diagram for a given reaction and facetpath

        Parameters
        ----------
        key: str
            a key to look up for a correct reactions, e.g.

            >>> key = 'Cu_111_OH_O+H'

        reaction_energies: Dict[str: float]
            a dictionary with reaction energies for a given facetpath and
            rxn_name, e.g.

            >>> reaction_energies = {'Cu_111_OH_O+H': 70.81}

        activation_barriers: Dict[str: Dict[str, str]]
            a dictionary with all barrier heights(kj/mol)
            in a format like below:

            >>> activation_barriers = {'Cu_111_OH_O+H':
                    {'TS_00': '155.27', 'TS_01': '157.97'}}

        axes: matplotlib.subplot object
            an axis for matplotlib.subplot
        num: int
            a num indicating column of the axes(an index of axes)
        plot_file_name: str
            provide a title for the plot, optional
        apply_max_barrier: bool
            specify whether to apply a filter for a max barrier,
            `False` by default

        '''
        if not plot_filename:
            plot_filename = 'plot.png'

        reaction_energy = float(reaction_energies[key])
        activation_barriers_rxn = activation_barriers[key]

        if max_barrier is not None:
            activation_barriers_rxn = {ts_name: float(barrier) for (
                ts_name, barrier) in activation_barriers_rxn.items()
                if float(barrier) < max_barrier}

        rxn_name_title = Results.rxn_title(rxn_name)
        energy_0 = 0
        rxn_ener_position = reaction_energy + 5
        rxn_ener_position_label = reaction_energy - 8

        reactants, products = rxn_name_title.split(' --> ')

        for ts_name, barrier in activation_barriers_rxn.items():
            barrier = float(barrier)
            x = np.arange(6)
            y = np.array([0, 0, barrier, barrier,
                          reaction_energy, reaction_energy])
            axes[num].plot(
                x, y, label='{}; {} kJ/mol'.format(ts_name, barrier))
            # axes[num].plot(
            #     x, y, label='{}'.format(ts_name))
            axes[num].hlines(barrier, 0, 2.0, linestyles='dotted')
            # add label with ener of the TS
            # barrier_position = barrier + 5
            # plt.annotate('{:.2f}'.format(barrier),
            #              (2.5, barrier_position), ha='center')

        # add lablel with the 0 ener for reactants
        axes[num].annotate('{:.2f}'.format(energy_0), (0.5, 5), ha='center')
        axes[num].annotate(reactants, (0.5, -8), ha='center')

        # add lablel with the reaction energy for products
        axes[num].annotate('{:.2f}'.format(reaction_energy),
                           (4.5, rxn_ener_position), ha='center')
        axes[num].annotate(
            products, (4.5, rxn_ener_position_label), ha='center')

        # plt.gca().axes.get_xaxis().set_visible(False)
        minor_locator = AutoMinorLocator(5)
        axes[num].margins(x=0)
        axes[num].xaxis.set_major_locator(plt.NullLocator())
        axes[num].yaxis.set_minor_locator(minor_locator)
        # labels
        axes[num].set_ylabel('E (kJ/mol)')
        axes[num].set_xlabel('reaction coordinate')

        axes[num].legend()
        axes[num].title.set_text(key)
        # plt.show()
        plt.tight_layout()
        figure = plt.gcf()
        figure.set_size_inches(x_size_in, y_size_in)
        plt.savefig(plot_filename)
