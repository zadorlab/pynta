from pynta.excatkit.molecule import Molecule
from pynta.excatkit.gratoms import Gratoms
from typing import Dict, List, Tuple


class GeneralTSGuessesGenerator():
    def __init__(self,
                 ts_est: str,
                 rxn: Dict[str, str],
                 rxn_name: str,
                 easier_to_build: str,
                 scfactor: float) -> None:
        '''[summary]

        Parameters
        ----------
        ts_est : str
            string representations of the species which is use as a TS guess
            skeleton - the one that is used to construct TS guess,
            e.g.
            'OH'
            e.g OH_O+H
        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        rxn_name : str
            a reaction name
            e.g OH_O+H
        easier_to_build : str
            a keys to ts_guess_skeleton; 'product' or 'reatant'
                scfator : float
            a scaling factor to scale a bond distance between
            atoms taking part in the reaction

        '''
        self.ts_est = ts_est
        self.rxn = rxn
        self.rxn_name = rxn_name
        self.easier_to_build = easier_to_build
        self.scfactor = scfactor
        self.ts_guess_skeleton = self.rxn[self.easier_to_build]
        self.reacting_species_connectivity = self.ts_guess_skeleton.split(
            '\n')

    def decide(
            self,
            reacting_idxs: List[int]) -> Tuple[Gratoms, int]:
        ''' Main method to decide to which family of reactions given
            reaction belongs

        Parameters
        ----------
        reacting_idxs : List[int]
            list with indicies of reacting atoms, in order as they appear in
            .yaml file but not necessary with the same indicies

        Returns
        -------
        ts_guess_el, s_bonded_idx : Tuple[Gratoms, int]
            ts_geuss_el is Gratoms representation of TS guess, whereas
            s_bonded_idx is the index of adsorbate atom that bonds to surface

        '''

        ts_guess_el, s_bonded_idx = None, None
        how_many_atoms_react = len(reacting_idxs)

        if how_many_atoms_react == 2:
            ts_guess_el, s_bonded_idx = Diatomic(
                self.ts_est,
                self.rxn,
                self.rxn_name,
                self.easier_to_build,
                self.scfactor).get_ts_guess_and_bonded_idx(reacting_idxs)

        elif how_many_atoms_react == 3:
            ts_guess_el, s_bonded_idx = Triatomic(
                self.ts_est,
                self.rxn,
                self.rxn_name,
                self.easier_to_build,
                self.scfactor).get_ts_guess_and_bonded_idx(reacting_idxs)

        return ts_guess_el, s_bonded_idx

    def build_ts_guess(self) -> Gratoms:
        ''' Convert ts_est string into a list of Gratoms objects.

            Numner of elements in the list depends is equal to number of
            distinct topologies available for given ts_est.

            For diatomics there will always be one element in the list.
            For other types, more than one structure is possible, e.g.
            ts_est = 'COH' --> ts_guess_list = ['COH' (sp), 'CHO' (sp2)]

        Returns
        -------
        ts_guess_list : List[Gratoms]
            a list of Gratoms objects with all distinct topologies for a given
            ts_est

        '''
        ts_guess_list = Molecule().molecule(self.ts_est)
        return ts_guess_list

    def get_s_bonded_idx(self) -> int:
        ''' Get index of the atoms that connects adsorbate to the surface

        Returns
        -------
        s_bonded_idx: int
            an index of an adsorbate atom that bonds it to surface

        Raises
        ------
        NotImplementedError
            if more then 1 atoms connects adsorbate to the surface - curently,
            only monodentate adsorbtion is supported.

        '''
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

        # surface reaction, at least one atom connected to the surface
        if s_bonded_idxs:
            if len(s_bonded_idxs) > 1:
                raise NotImplementedError('Only monodendate type of adsorbtion'
                                          ' is currently supported.')
            else:
                s_bonded_idx = int(s_bonded_idxs[0]) - 1

        # gas phase reaction
        else:
            s_bonded_idx = self.get_the_most_connected_atom()
        return s_bonded_idx

    def get_the_most_connected_atom(self) -> int:
        '''[summary]

        Returns
        -------
        max_conected_atom_idx : int
            an index of an atom that has the most connections to other atoms

        '''
        number_of_connections = {}
        n_surf_at_befor_ads = 0

        for line in (self.reacting_species_connectivity):
            if 'X' in line:
                n_surf_at_befor_ads += 1
            else:
                break

        for num, line in enumerate(self.reacting_species_connectivity):
            if 'X' not in line:
                number_of_connections[num] = line.count('{')

        for tmp_idx, n_connect in number_of_connections.items():
            if n_connect == max(number_of_connections.values()):
                max_conected_atom_idx = (tmp_idx - n_surf_at_befor_ads - 1)
                return max_conected_atom_idx


class Diatomic(GeneralTSGuessesGenerator):
    def get_ts_guess_and_bonded_idx(
            self,
            reacting_idxs: List[int]) -> Tuple[Gratoms, int]:
        '''[summary]

        Parameters
        ----------
        reacting_idxs : List[int]
            list with indicies of reacting atoms, in order as they appear in
            .yaml file but not necessary with the same indicies

        Returns
        -------
        ts_guess_el, s_bonded_idx : Tuple[Gratoms, int]
            ts_geuss_el is Gratoms representation of TS guess, whereas
            s_bonded_idx is the index of adsorbate atom that bonds to surface

        '''
        # Convert adsorbate (string) to a list of Gratoms object.
        ts_guess_list = self.build_ts_guess()

        # For two atoms in adsorbat, there is only one possible topology
        ts_guess_el = ts_guess_list[0]
        s_bonded_idx = self.get_s_bonded_idx()

        react_atom_idx_1, react_atom_idx_2 = reacting_idxs
        bondlen = ts_guess_el.get_distance(react_atom_idx_1, react_atom_idx_2)

        n_total_ads_atoms = len(ts_guess_el)

        # edge cases
        if n_total_ads_atoms == 2:
            ts_guess_el.rotate(90, 'y')

        elif n_total_ads_atoms == 3:
            remaining_atom_idx = n_total_ads_atoms - \
                (react_atom_idx_1 + react_atom_idx_2)

            if react_atom_idx_2 != 2:
                # Symetrically not important, but structure looks
                # better visually
                ts_guess_el.rotate(-90, 'z')
            else:
                ts_guess_el.rotate(90, 'z')

            # set angle that puts react_atom_idx_2 closer to the surface
            ts_guess_el.set_angle(
                remaining_atom_idx,
                react_atom_idx_1,
                react_atom_idx_2,
                30,
                indices=[remaining_atom_idx,
                         react_atom_idx_1,
                         react_atom_idx_2],
                add=True)

        else:
            # continue with defaults
            pass

        # scale the bond distance between reacting part
        ts_guess_el.set_distance(react_atom_idx_1, react_atom_idx_2,
                                 bondlen * self.scfactor, fix=0)

        return ts_guess_el, s_bonded_idx


class Triatomic(GeneralTSGuessesGenerator):
    def get_ts_guess_and_bonded_idx(self):
        raise NotImplementedError('Only diatomic reactions at this moment')
