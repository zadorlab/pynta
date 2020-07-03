from ase.io import read
from pathlib import Path
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline


class Results():
    def __init__(self):
        self.ev_to_kjmol = 23.06035 * 4.184

    def get_reaction_energy(self, minima_path, reactants_list, products_list,
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
        reaction_energy : float
            an energy of reaction calculated a difference between
            the most stable product and the most stable reactant

        '''
        r_ener_list = []
        p_ener_list = []
        for reactant in reactants_list:
            lowest_reactant_ener = Results.get_lowest_species_ener(self,
                                                                   minima_path,
                                                                   reactant)
            r_ener_list.append(lowest_reactant_ener)
        for product in products_list:
            lowest_product_ener = Results.get_lowest_species_ener(self,
                                                                  minima_path,
                                                                  product)
            p_ener_list.append(lowest_product_ener)

        slab_ener = Results.get_slab_ener(self, slab_path)

        # Depending how the reactants and products are defined,
        # there are three options here:
        # e.g. OH --> O + H
        if len(p_ener_list) > len(r_ener_list):
            reaction_energy = (sum(p_ener_list) - slab_ener -
                               r_ener_list) * self.ev_to_kjmol
        # e.g. O + H --> OH
        elif len(p_ener_list) < len(r_ener_list):
            reaction_energy = (p_ener_list - slab_ener -
                               sum(r_ener_list)) * self.ev_to_kjmol
        # e.g. A + B --> C + D
        else:
            reaction_energy = (sum(p_ener_list) -
                               sum(r_ener_list)) * self.ev_to_kjmol
        reaction_energy = '{:.2f}'.format(round(reaction_energy[0], 3))
        return reaction_energy

    def get_barrier(self, minima_path, ts_path, species):
        ''' Calculate reaction energy relatively to the most stable reactant

        Parameters:
        ___________
        minima_path : str
            a path to main minima directory
            e.g. Cu_111/minima
        species : str
            a symbol of the species
            e.g. 'OH', 'H', 'O'
        ts_path : str
            a path to main TS directory
            e.g. 'Cu_111/TS_estimate_unique'

        Returns:
        ________

        activation_barriers : dict('str'='str')
            a dictonary with all barrier heights relatively
            to the most stable reactant (in kJ/mol)

        '''
        lowest_reactant_ener = Results.get_lowest_species_ener(self,
                                                               minima_path,
                                                               species)
        tss_ener = Results.get_ts_ener(self, ts_path)
        tss_name = Results.format_TS_name(self, ts_path)

        activation_barriers = {}
        for ts_ener, ts_name in zip(tss_ener, tss_name):
            barrier = (ts_ener - lowest_reactant_ener) * self.ev_to_kjmol
            activation_barriers['TS_' + ts_name] = '{:.2f}'.format(barrier)
        return activation_barriers

    def get_slab_ener(self, slab_path):
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

    def get_ts_ener(self, ts_path):
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
        tss = Results.get_ts_out_files(self, ts_path)
        for ts in tss:
            with open(ts, 'r') as f:
                data = f.readlines()
                enerLine = data[-1]
                enerVal = enerLine.split()
                ts_ener_dict[ts] = float(enerVal[3])
        ts_ener_list = list(ts_ener_dict.values())
        return ts_ener_list

    def get_lowest_species_ener(self, minima_path, species):
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
        species_out_file_path_list = Results.get_species_out_files(self,
                                                                   minima_path,
                                                                   species)
        for spiecies_out_file_path in species_out_file_path_list:
            with open(spiecies_out_file_path, 'r') as f:
                data = f.readlines()
                enerLine = data[-1]
                enerVal = enerLine.split()
                species_ener_dict[spiecies_out_file_path] = float(enerVal[3])
                f.close()
        lowest_species_ener = min(species_ener_dict.values())
        return lowest_species_ener

    def get_ts_out_files(self, ts_path):
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

    def get_species_out_files(self, minima_path, species):
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
            a list with paths to all minima Sella's *out files for given species
            e.g. ['Cu_111/minima/OH_01_relax.out', 'Cu_111/minima/OH_00_relax.out',
            'Cu_111/minima/OH_03_relax.out', 'Cu_111/minima/OH_02_relax.out']

        '''
        species_out_file_path_list = []
        species = species + '_'
        outfile = '{}*out'.format(species)
        reactant_out_list = Path(minima_path).glob(outfile)
        for reactant_out_file in reactant_out_list:
            species_out_file_path_list.append(str(reactant_out_file))
        return species_out_file_path_list

    def format_TS_name(self, ts_path):
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
        ts_out_file_list = Results.get_ts_out_files(self, ts_path)
        for ts_out_file in ts_out_file_list:
            prefix = os.path.split(ts_out_file)[1].split('_')[0]
            prefix_list.append(prefix)
        return prefix_list

    def plot(self, activation_barriers, reaction_energy):
        ''' Plot reaction energy diagram
        
        Parameters:
        ___________
        activation_barriers : dict('str'='str')
            a dictonary with all barrier heights relatively
            to the most stable reactant (in kJ/mol)
        reaction_energy : float

        '''
        for ts_name, barrier in activation_barriers.items():
            # x = [1, 2, 3]
            x = np.arange(6)
            y = np.array([0, 0, float(barrier), float(barrier),
                          float(reaction_energy, float(reaction_energy)])
            # x_new = np.linspace(1, 5, 50)
            # a_BSpline = make_interp_spline(x, y)
            # y_new = a_BSpline(x_new)
            plt.plot(x, y, label=ts_name)
        plt.legend()
        plt.title('OH --> O + H')
        plt.show()
