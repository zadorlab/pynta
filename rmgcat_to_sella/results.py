import re
from ase.io import read
from pathlib import Path
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
from numpy.core.fromnumeric import prod
from rmgcat_to_sella.io import IO


class Results():
    def __init__(
            self,
            ts_path,
            slab_path,
            yamlfile,
            facetpaths):
        '''
        Parameters:
        ___________
        ts_path : str
            a path to main TS directory
            e.g. 'Cu_111/TS_estimate_unique'
        slab_path : str
            a path to the slab
            e.g. 'Cu_111_slab_opt.xyz'
        reactant_list : list(str)
            a list with all reactants
            e.g. ['OH']
        product_list : list(str)
            a list with all products
            e.g. ['O', 'H']

        '''
        # self.minima_path = minima_path
        self.ts_path = ts_path
        self.slab_path = slab_path
        self.yamlfile = yamlfile
        self.facetpaths = facetpaths
        self.reactions = IO().open_yaml_file(self.yamlfile)
        self.ev_to_kjmol = 23.06035 * 4.184

    def get_reaction_energies_all(self):
        reaction_energies = {}
        for facetpath, slab_path in zip(self.facetpaths, self.slab_path):
            minima_path = os.path.join(facetpath, 'minima')
            for rxn in self.reactions:
                r_name_list, p_name_list, _ = IO().prepare_react_list(rxn)
                rxn_name = IO().get_rxn_name(rxn)
                key = facetpath+'_'+rxn_name
                reaction_energies[key] = self.get_reaction_energy(
                    minima_path, facetpath, r_name_list, p_name_list,
                    slab_path)
        return reaction_energies

    def get_reaction_energy(
            self,
            minima_path,
            facetpath,
            r_name_list,
            p_name_list,
            slab_path):
        ''' Calclate reaction energy as a difference between
            the most stable product and the most stable reactant

        Parameters:
        ___________
        minima_path : str
            a path to main minima directory
            e.g. Cu_111/minima
        reactant_list : list(str)
            a list with all reactants
            e.g. ['OH']
        product_list : list(str)
            a list with all products
            e.g. ['O', 'H']
        slab_path : str
            a path to the slab
            e.g. 'Cu_111_slab_opt.xyz'

        Returns:
        ________
        reaction_energy : str
            an energy of reaction calculated a difference between
            the most stable product and the most stable reactant

        '''
        r_ener_list, p_ener_list, slab_ener, nslabs = self.get_data(
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
        # print(reaction_energy)
        reaction_energy = '{:.2f}'.format(round(reaction_energy[0], 3))
        return reaction_energy

    def get_barrier_all(self):
        ts_ener = {}
        for facetpath, slab_path in zip(self.facetpaths, self.slab_path):
            minima_path = os.path.join(facetpath, 'minima')
            for rxn in self.reactions:
                r_name_list, p_name_list, _ = IO().prepare_react_list(rxn)
                rxn_name = IO().get_rxn_name(rxn)
                ts_path = os.path.join(
                    facetpath, rxn_name, 'TS_estimate_unique')
                ts_ener[facetpath+'_'+rxn_name] = self.get_barrier(
                    minima_path, ts_path, facetpath, r_name_list, p_name_list,
                    slab_path)
        return ts_ener

    def get_barrier(
            self,
            minima_path,
            ts_path,
            facetpath,
            r_name_list,
            p_name_list,
            slab_path):
        ''' Calculate reaction energy relatively to the most stable reactant

        Parameters:
        ___________
        minima_path : str
            a path to main minima directory
            e.g. Cu_111/minima
        ts_path : str
            a path to main TS directory
            e.g. 'Cu_111/TS_estimate_unique'
        reactant_list : list(str)
            a list with all reactants
            e.g. ['OH']
        product_list : list(str)
            a list with all products
            e.g. ['O', 'H']
        slab_path : str
            a path to the slab
            e.g. 'Cu_111_slab_opt.xyz'


        Returns:
        ________

        activation_barriers : dict('str'=float)
            a dictonary with all barrier heights relatively
            to the most stable reactant (in kJ/mol)

        '''

        r_ener_list, _, slab_ener, _ = self.get_data(
            minima_path, facetpath, r_name_list, p_name_list, slab_path)
        tss_ener = self.get_ts_ener(ts_path)
        tss_name = self.format_TS_name(ts_path)

        activation_barriers = {}
        for ts_ener, ts_name in zip(tss_ener, tss_name):
            # Should be valud for any type of TS, no matter how many reactants
            # take part in the reaction, i.e. A + B -> AB or AB -> A + B
            barrier = (ts_ener + slab_ener * (len(r_ener_list) - 1) -
                       sum(r_ener_list)) * self.ev_to_kjmol
            activation_barriers['TS_' +
                                ts_name] = '{:.2f}'.format(barrier)
        return activation_barriers

    def get_data(
            self,
            minima_path,
            facetpath,
            r_name_list,
            p_name_list,
            slab_path):
        ''' Returns the lowest energies lists for reactants and products.

        Parameters:
        ___________
        minima_path : str
            a path to main minima directory
            e.g. Cu_111/minima
        reactant_list : list(str)
            a list with all reactants
            e.g. ['OH']
        product_list : list(str)
            a list with all products
            e.g. ['O', 'H']
        slab_path : str
            a path to the slab
            e.g. 'Cu_111_slab_opt.xyz'


        Returns:
        ________
        r_ener_list : list(float)
            a list with the lowest energy conformer for each reactant
        p_ener_list : list(float)
            a list with the lowest energy conformer for each products
        slab_ener : float
            an energy of the representative slab
            of the size the same as reactants, TS, or product
        nslabs : int
            a number specifying how many additional slabs have to be considered
            in a stoichiometric reaction. Its defined as the absoulute value
            of the difference between amount of products and reactants
            e.g.
            O + H --> OH (1 slab)
            C + O + H --> COH (2 slabs)

        '''
        r_ener_list = []
        p_ener_list = []
        # get the lowest energy for all reactants
        for reactant in r_name_list:
            # bug to be fixed
            if reactant == 'OH':
                reactant = 'HO'
            lowest_reactant_ener = self.get_lowest_species_ener(
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
            if product == 'OH':
                product = 'HO'
            lowest_product_ener = self.get_lowest_species_ener(
                minima_path, product, facetpath)
            p_ener_list.append(lowest_product_ener)
        # check if .out files for products are copied
        if None in p_ener_list:
            print('----')
            print('Found None in products energy list. Missing .out files'
                  'for products.')
            print('----')
            raise TypeError

        slab_ener = self.get_slab_ener(slab_path)
        nslabs = abs(len(p_ener_list) - len(r_ener_list))

        return r_ener_list, p_ener_list, slab_ener, nslabs

    def get_slab_ener(
            self,
            slab_path):
        ''' Get energy of the slab

        Parameters:
        ___________
        slab_path : str
            a path to the slab
            e.g. 'Cu_111_slab_opt.xyz'

        Returns:
        ________
        slab_ener : float
            an energy of the slab in eV

        '''
        slab = read(slab_path)
        slab_ener = slab.get_potential_energy()
        return slab_ener

    def get_ts_ener(
            self,
            ts_path):
        ''' Get energy of all TSs

        Parameters:
        ___________
        ts_path : str
            a path to main TS directory
            e.g. 'Cu_111/TS_estimate_unique'

        Returns:
        ________
        ts_ener_list : list(float)
            a sorted list with energy value for each TS

        '''
        ts_ener_dict = {}
        tss = self.get_ts_out_files(ts_path)
        for ts in tss:
            with open(ts, 'r') as f:
                data = f.readlines()
                enerLine = data[-1]
                enerVal = enerLine.split()
                ts_ener_dict[ts] = float(enerVal[3])
        ts_ener_list = list(ts_ener_dict.values())
        return ts_ener_list

    def get_lowest_species_ener(
            self,
            minima_path,
            species,
            facetpath):
        ''' Get the lowest energy of the most stable species

        Parameters:
        ___________
        minima_path : str
            a path to main minima directory
            e.g. Cu_111/minima
        species : str
            a symbol of the species
            e.g. 'OH', 'H', 'O'

        Returns:
        ________
        lowest_species_ener : float
            conformer of the lowest energy among all
            calculated for the given species

        '''
        species_ener_dict = {}
        try:
            species_out_file_path_list = self.get_species_out_files(
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

    def get_ts_out_files(
            self,
            ts_path):
        ''' Get TS .out files

        Parameters:
        ___________
        ts_path : str
            a path to main TS directory
            e.g. 'Cu_111/TS_estimate_unique'

        Returns:
        ________
        ts_out_file_list : list(str)
            a sorted list with paths to Sella's .out files for each TSs

        '''
        ts_out_file_list = []
        ts_file_list = Path(ts_path).glob('*out')
        for ts_out_file in ts_file_list:
            ts_out_file_list.append(str(ts_out_file))
        return sorted(ts_out_file_list)

    def get_species_out_files(
            self,
            minima_path,
            species,
            facetpath):
        ''' Get .out files for each reactants

        Parameters:
        ___________
        minima_path : str
            a path to main minima directory
            e.g. Cu_111/minima
        species : str
            a symbol of the species
            e.g. 'OH', 'H', 'O'

        Returns:
        ________
        species_out_file_path_list : list(str)
            a list with paths to all minima Sella's *out files for given
            species
            e.g. ['Cu_111/minima/OH_01_relax.out',
            'Cu_111/minima/OH_00_relax.out',
            'Cu_111/minima/OH_03_relax.out', 'Cu_111/minima/OH_02_relax.out']

        '''
        species_out_file_path_list = []
        species = species + '_'
        outfile = '{}_{}*out'.format(facetpath, species)
        reactant_out_list = Path(minima_path).glob(outfile)
        for reactant_out_file in reactant_out_list:
            species_out_file_path_list.append(str(reactant_out_file))
        return species_out_file_path_list

    def format_TS_name(
            self,
            ts_path):
        ''' Function to get prefixes of TSs

        Parameters:
        ___________
        ts_path : str
            a path to main TS directory
            e.g. 'Cu_111/TS_estimate_unique'

        Returns:
        ________
        prefix_list : list(str)
            a list with all prefixes for TSs
        '''
        prefix_list = []
        ts_out_file_list = self.get_ts_out_files(ts_path)
        for ts_out_file in ts_out_file_list:
            prefix = os.path.split(ts_out_file)[1].split('_')[0]
            prefix_list.append(prefix)
        return prefix_list

    def rxn_title(
            self,
            rxn_name):
        ''' Return rxn name with arrow between reactants and products'''
        reactants, products = rxn_name.split('_')
        rxn_name_title = reactants + ' --> ' + products
        return rxn_name_title

    def plot(self):
        reaction_energies = self.get_reaction_energies_all()
        activation_barriers = self.get_barrier_all()
        _, (ax1, ax2) = plt.subplots(2, 2)
        for num, facetpath in enumerate(self.facetpaths):
            for ax, rxn, in zip((ax1, ax2), self.reactions):
                rxn_name = IO().get_rxn_name(rxn)
                key = facetpath+'_'+rxn_name
                self.plot_rxn(key, reaction_energies,
                              activation_barriers, rxn_name, ax, num)

    def plot_rxn(
            self,
            key,
            reaction_energies,
            activation_barriers,
            rxn_name,
            axes,
            num,
            plot_filename=None,
            apply_max_barrier=False):
        ''' Plot reaction energy diagram

        Parameters:
        ___________
        plot_title : str
            provide a title for the plot, optional
        apply_max_barrier : bool
            specify whether to apply a filter for a max barrier,
            default = False

        '''
        if not plot_filename:
            plot_filename = 'plot.png'

        reaction_energy = float(reaction_energies[key])
        activation_barriers_rxn = activation_barriers[key]

        if apply_max_barrier:
            activation_barriers_rxn = {ts_name: float(barrier) for (
                ts_name, barrier) in activation_barriers_rxn.items()
                if float(barrier) < 300}

        rxn_name_title = self.rxn_title(rxn_name)
        energy_0 = 0
        rxn_ener_position = reaction_energy + 5
        rxn_ener_position_label = reaction_energy - 8

        reactants, products = rxn_name_title.split(' --> ')

        for ts_name, barrier in activation_barriers_rxn.items():
            barrier = float(barrier)
            x = np.arange(6)
            y = np.array([0, 0, barrier, barrier,
                          reaction_energy, reaction_energy])
            axes[num].plot(x, y, label=ts_name)
            axes[num].hlines(barrier, 0, 2.0, linestyles='dotted')
            # add label with ener of the TS
            # barrier_position = barrier + 5
            # plt.annotate('{:.2f}'.format(barrier),
            #              (2.5, barrier_position), ha='center')

        # # add lablel with the 0 ener for reactants
        axes[num].annotate('{:.2f}'.format(energy_0), (0.5, 5), ha='center')
        axes[num].annotate(reactants, (0.5, -8), ha='center')

        # # add lablel with the reaction energy for products
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
        figure.set_size_inches(10, 10)
        plt.savefig(plot_filename)
