from pynta.excatkit.molecule import Molecule
from pynta.excatkit.gratoms import Gratoms
from typing import Dict, List, Tuple


class GeneralTSGuessesGenerator():
    def __init__(self,
                 ts_est,
                 rxn,
                 rxn_name,
                 easier_to_build,
                 scfactor):
        self.ts_est = ts_est
        self.rxn = rxn
        self.rxn_name = rxn_name
        self.easier_to_build = easier_to_build
        self.scfactor = scfactor
        self.reacting_species_connectivity = self.rxn[self.easier_to_build].split(
            '\n')

    def get_ts_guess_and_bonded_idx(self):
        ts_guess_list = self.build_ts_guess()
        # For diatimics, there is only one possible topology, so
        ts_guess_el = ts_guess_list[0]

        print(self.reacting_species_connectivity)

        s_bonded_idx, reacting_idxs = self.get_bondend_and_reacting_idxs()

        react_ind_1, react_ind_2 = reacting_idxs

        bondlen = ts_guess_el.get_distance(react_ind_1, react_ind_2)
        ts_guess_el.rotate(90, 'y')
        ts_guess_el.set_distance(
            react_ind_1, react_ind_2, bondlen * self.scfactor, fix=0)
        # print(ts_guess_el, s_bonded_idx)
        # for atom in ts_guess_el:
        #     print(atom)
        return ts_guess_el, s_bonded_idx

    def build_ts_guess(self) -> Gratoms:
        ''' Convert ts_est string into a list of Gratoms objects.

            Numner of elements in the list depends is equal to number of
            distinct topologies available for given ts_est.

            For diatomics there will always be one element in the list.
            For other types, more than one structure is possible, e.g.
            ts_est = 'COH' --> ts_guess_list = ['COH' (sp), 'CHO' (sp2)]

        Parameters
        ----------
        ts_est : str
            a string representing species that will be used to get ts_guess
            e.g. 'OH', 'COH'

        Returns
        -------
        ts_guess_list : List[Gratoms]
            a list of Gratoms objects with all distinct topologies for a given
            ts_est

        '''
        ts_guess_list = Molecule().molecule(self.ts_est)
        return ts_guess_list

    def get_bondend_and_reacting_idxs(self) -> int:
        ''' Get an index of surface bonded atom

        Parameters
        ----------
        ts_guess_el : Gratoms
            a Gratom object of ts_guess with the chosen topology, if more than
            one topologies are possible
        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        reacting_sp : str
            a key to rxn dictionary
            'reactant' or 'product' are the only avaiable options options

        Returns
        -------
        s_bonded_idx : int
            an int with index of atom bonded to the surface

        Raises
        ------
        NotImplementedError
            when there are more than one atoms connected to the surface

        '''
        s_bonded_idxs = self.get_s_bonded_idx()
        reacting_atoms_idx = self.get_reacting_atoms_idx()

        return s_bonded_idxs, reacting_atoms_idx

    def get_s_bonded_idx(self) -> int:
        surface_indicies = []
        s_bonded_idxs = []
        for line in self.reacting_species_connectivity:
            if 'X' in line:
                index = line.split()[0]
                surface_indicies.append(index)
        for index in surface_indicies:
            keyphrase = '{' + '{}'.format(index)
            for line in self.reacting_species_connectivity:
                if keyphrase in line:
                    s_bonded_idxs.append(line.split()[0])
        if len(s_bonded_idxs) > 1:
            raise NotImplementedError('Only monodendate type of adsorbtion is '
                                      'currently supported.')
        return int(s_bonded_idxs[0]) - 1

    def get_reacting_atoms_idx(self):
        reacting_idxs = []
        for num, line in enumerate(self.reacting_species_connectivity):
            if '*' in line and 'X' not in line:
                reacting_idxs.append(num - 1)
        return reacting_idxs
