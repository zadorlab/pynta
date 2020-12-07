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

    def build_ts_guess(ts_est: str) -> Gratoms:
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
        ts_guess_list = Molecule().molecule(ts_est)
        return ts_guess_list

    def get_surface_bonded_atom_idx(self) -> int:
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
        surface_bonded_atom_idx : int
            an int with index of atom bonded to the surface

        Raises
        ------
        NotImplementedError
            when there are more than one atoms connected to the surface

        '''
        reacting_sp_connectivity = self.rxn[self.easier_to_build].split('\n')
        surface_indicies = []
        surface_bonded_atom_idxs = []
        for line in reacting_sp_connectivity:
            if 'X' in line:
                index = line.split()[0]
                surface_indicies.append(index)
        for index in surface_indicies:
            keyphrase = '{' + '{}'.format(index)
            for line in reacting_sp_connectivity:
                if keyphrase in line:
                    surface_bonded_atom_idxs.append(line.split()[0])
        if len(surface_bonded_atom_idxs) > 1:
            raise NotImplementedError('Only monodendate type of adsorbtion is '
                                      'currently supported.')

        print(int(surface_bonded_atom_idxs[0]) - 1)
        return int(surface_bonded_atom_idxs[0]) - 1
