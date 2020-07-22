from ase.io import read
from pathlib import Path
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
# from scipy.interpolate import make_interp_spline


class Results():
    def __init__(self, minima_path, ts_path, slab_path,
                 reactants_list, products_list):
        '''
        Parameters:
        ___________
        minima_path : str
            a path to main minima directory
            e.g. Cu_111/minima
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
        self.minima_path = minima_path
        self.ts_path = ts_path
        self.slab_path = slab_path
        self.reactants_list = reactants_list
        self.products_list = products_list
        self.ev_to_kjmol = 23.06035 * 4.184

    def get_reaction_energy(self):
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
        r_ener_list, p_ener_list, slab_ener, nslabs = Results.get_data(self)
        # Depending how the reactants and products are defined,
        # there are three options here:
        # e.g. OH --> O + H
        if len(p_ener_list) > len(r_ener_list):
            reaction_energy = (sum(p_ener_list) - slab_ener * nslabs -
                               r_ener_list) * self.ev_to_kjmol
        # e.g. O + H --> OH
        elif len(p_ener_list) < len(r_ener_list):
            reaction_energy = (p_ener_list + slab_ener * nslabs -
                               sum(r_ener_list)) * self.ev_to_kjmol
        # e.g. A + B --> C + D
        else:
            reaction_energy = (sum(p_ener_list) -
                               sum(r_ener_list)) * self.ev_to_kjmol
        reaction_energy = '{:.2f}'.format(round(reaction_energy[0], 3))
        return reaction_energy

    def get_barrier(self):
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

        r_ener_list, p_ener_list, slab_ener, nslabs = Results.get_data(self)
        tss_ener = Results.get_ts_ener(self)
        tss_name = Results.format_TS_name(self)

        activation_barriers = {}
        for ts_ener, ts_name in zip(tss_ener, tss_name):
            # Depending how the reactants and products are defined,
            # there are three options here:
            # e.g. OH --> O + H
            if len(p_ener_list) >= len(r_ener_list):
                barrier = (ts_ener - sum(r_ener_list)) * self.ev_to_kjmol
                activation_barriers['TS_' + ts_name] = '{:.2f}'.format(barrier)
                # activation_barriers['TS_' + ts_name] = barrier
            # e.g. O + H --> OH
            elif len(p_ener_list) < len(r_ener_list):
                barrier = (ts_ener + slab_ener * nslabs -
                           sum(r_ener_list)) * self.ev_to_kjmol
                activation_barriers['TS_' + ts_name] = '{:.2f}'.format(barrier)
                # activation_barriers['TS_' + ts_name] = barrier
            else:
                raise NotImplementedError(
                    'Not tested if r_ener_list=p_ener_list')
        return activation_barriers

    def get_data(self):
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
        for reactant in self.reactants_list:
            lowest_reactant_ener = Results.get_lowest_species_ener(
                self, reactant)
            r_ener_list.append(lowest_reactant_ener)
        for product in self.products_list:
            lowest_product_ener = Results.get_lowest_species_ener(
                self, product)
            p_ener_list.append(lowest_product_ener)

        slab_ener = Results.get_slab_ener(self)
        nslabs = abs(len(p_ener_list) - len(r_ener_list))

        return r_ener_list, p_ener_list, slab_ener, nslabs

    def get_slab_ener(self):
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
        slab = read(self.slab_path)
        slab_ener = slab.get_potential_energy()
        return slab_ener

    def get_ts_ener(self):
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
        tss = Results.get_ts_out_files(self)
        for ts in tss:
            with open(ts, 'r') as f:
                data = f.readlines()
                enerLine = data[-1]
                enerVal = enerLine.split()
                ts_ener_dict[ts] = float(enerVal[3])
        ts_ener_list = list(ts_ener_dict.values())
        return ts_ener_list

    def get_lowest_species_ener(self, species):
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
            species_out_file_path_list = Results.get_species_out_files(self,
                                                                       species)
            for spiecies_out_file_path in species_out_file_path_list:
                with open(spiecies_out_file_path, 'r') as f:
                    data = f.readlines()
                    enerLine = data[-1]
                    enerVal = enerLine.split()
                    species_ener_dict[spiecies_out_file_path] = float(
                        enerVal[3])
                    f.close()
            lowest_species_ener = min(species_ener_dict.values())
            return lowest_species_ener
        except ValueError:
            print('Check if minima .out files copied successfully')

    def get_ts_out_files(self):
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
        ts_file_list = Path(self.ts_path).glob('*out')
        for ts_out_file in ts_file_list:
            ts_out_file_list.append(str(ts_out_file))
        return sorted(ts_out_file_list)

    def get_species_out_files(self, species):
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
        outfile = '{}*out'.format(species)
        reactant_out_list = Path(self.minima_path).glob(outfile)
        for reactant_out_file in reactant_out_list:
            species_out_file_path_list.append(str(reactant_out_file))
        return species_out_file_path_list

    def format_TS_name(self):
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
        ts_out_file_list = Results.get_ts_out_files(self)
        for ts_out_file in ts_out_file_list:
            prefix = os.path.split(ts_out_file)[1].split('_')[0]
            prefix_list.append(prefix)
        return prefix_list

    def rxn_title(self):
        ''' Return rxn name with arrow between reactants and products'''
        reactants = '+'.join([str(species) +
                              '*' for species in self.reactants_list])
        products = '+'.join([str(species) +
                             '*' for species in self.products_list])
        rxn_name = reactants + ' --> ' + products
        return rxn_name

    def plot(self, plot_title=None, plot_filename=None,
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

        reaction_energy = float(Results.get_reaction_energy(self))
        activation_barriers = Results.get_barrier(self)

        if apply_max_barrier:
            activation_barriers = {ts_name: float(barrier) for (
                ts_name, barrier) in activation_barriers.items()
                if float(barrier) < 300}

        rxn_name = Results.rxn_title(self)
        energy_0 = 0
        rxn_ener_position = reaction_energy + 5
        rxn_ener_position_label = reaction_energy - 8

        reactants, products = rxn_name.split(' --> ')

        for ts_name, barrier in activation_barriers.items():
            barrier = float(barrier)
            x = np.arange(6)
            y = np.array([0, 0, barrier, barrier,
                          reaction_energy, reaction_energy])
            plt.plot(x, y, label=ts_name)
            plt.hlines(barrier, 0, 2.0, linestyles='dotted')
            # add label with ener of the TS
            # barrier_position = barrier + 5
            # plt.annotate('{:.2f}'.format(barrier),
            #              (2.5, barrier_position), ha='center')

        # add lablel with the 0 ener for reactants
        plt.annotate('{:.2f}'.format(energy_0), (0.5, 5), ha='center')
        plt.annotate(reactants, (0.5, -8), ha='center')

        # add lablel with the reaction energy for products
        plt.annotate('{:.2f}'.format(reaction_energy),
                     (4.5, rxn_ener_position), ha='center')
        plt.annotate(products, (4.5, rxn_ener_position_label), ha='center')

        ax = plt.axes()
        # plt.gca().axes.get_xaxis().set_visible(False)
        minor_locator = AutoMinorLocator(5)
        plt.margins(x=0)
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_minor_locator(minor_locator)
        # labels
        ax.set_ylabel('E (kJ/mol)')
        ax.set_xlabel('reaction coordinate')

        plt.legend()
        # plt.title(rxn_name)
        plt.title(plot_title)
        # plt.show()
        plt.tight_layout()
        plt.savefig(plot_filename)

    # def get_latex_table(self):
    #     reaction_energy = Results.get_reaction_energy(self)
    #     activation_barriers = Results.get_barrier(self)
    #     rxn_name = Results.rxn_title(self)
    #     var = 1
    #     bra = '{'
    #     ket = '}'

    #     for i in range(10):
    #         print('{0}This is {2} in brackets{1}'.format(bra, ket, var))
