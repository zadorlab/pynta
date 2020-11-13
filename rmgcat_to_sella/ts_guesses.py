from rmgcat_to_sella.excatkit.molecule import Molecule
from rmgcat_to_sella.excatkit.gratoms import Gratoms
from typing import List, Dict, Tuple
from collections import Counter
from ase.io import write


class TSGuessesGenerator():
    def decide(self,
               ts_est,
               rxn,
               rxn_name,
               reacting_sp,
               reacting_species, scfactor):
        ts_guess_list = TSEstimator.build_ts_guess(ts_est)
        ts_guess_check_atomicity = ts_guess_list[0]
        n_atoms = len(ts_guess_check_atomicity)

        if n_atoms == 2:
            print('Reaction {} is a diatomic reaction'.format(rxn_name))

            # get ts_guess (Gratom) and index of bonded atom (int)
            ts_guess, bonded_idx = Diatomic().get_ts_guess_and_bonded_idx(
                ts_est, rxn, reacting_sp, reacting_species, scfactor)

        elif n_atoms == 3:
            print('Reaction {} is a triatomic reaction'.format(rxn_name))

            ts_guess, bonded_idx = Triatomic().get_ts_guess_and_bonded_idx(
                ts_est, rxn, reacting_sp, reacting_species, scfactor)

        elif n_atoms == 4:
            print('Reaction {} is a tetraatomic reaction'.format(rxn_name))

            ts_guess, bonded_idx = Tetraatomic().get_ts_guess_and_bonded_idx(
                ts_est, rxn, reacting_sp, reacting_species, scfactor)

        else:
            raise NotImplementedError('Only di-, tri- and some tetraatomic '
                                      'reactions are supported at this moment')
        return ts_guess, bonded_idx


class TSEstimator():
    @staticmethod
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

    @staticmethod
    def get_surface_bonded_atom_idx(
            rxn: Dict[str, str],
            reacting_sp: str) -> int:
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
        reacting_sp_connectivity = rxn[reacting_sp].split('\n')
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

        return int(surface_bonded_atom_idxs[0]) - 1

    @staticmethod
    def convert_symbol_to_index(
            bonded: str,
            ts_guess_el: Gratoms) -> int:
        ''' Convert symbol of atom bonded to the surface to its index
            in a given ts_est

        Parameters
        ----------
        bonded : str
            a chemical symbol of surface bonded atom
        ts_guess_el : Gratoms
            a Gratom object of ts_guess with the chosen topology, if more than
            one topologies are possible

        Returns
        -------
        surface_bonded_atom_idx : int
            an int with index of atom bonded to the surface

        '''
        symbol = str(ts_guess_el.symbols)
        surface_bonded_atom_idx = symbol.find(bonded)
        return surface_bonded_atom_idx

    @staticmethod
    def get_atomic_connections(
            rxn: Dict[str, str],
            reacting_sp: str) -> Dict[str, int]:
        ''' Get a dictionary describing how mamny connections exists for
            a given atom of a structure that is used to estimate TS structure


        Parameters
        ----------
        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        reacting_sp : str
            a key to rxn dictionary
            'reactant' or 'product' are the only avaiable options options

        Returns
        -------
        atomic_connections: Dict[str, int]
            a dict with connectivity info for a given structure that is used
            to estimate TS structure

        '''
        atomic_connections = {}
        reacting_sp_connectivity = rxn[reacting_sp].split('\n')
        for line in reacting_sp_connectivity:
            connections = line.count('{')
            symbol = line.split()[2]
            atomic_connections[symbol] = connections
        return atomic_connections

    @staticmethod
    def get_reacting_atoms_indices(
            ts_guess_el: Gratoms,
            reacting_atoms: List[str]) -> Dict[str, int]:
        ''' Convert a list of reacting atoms (chemical symbols) into a dict
            where key is a reacting atom and value is its index

        Parameters
        ----------
        ts_guess_el : Gratoms
            a Gratom object of ts_guess with the chosen topology, if more than
            one topologies are possible
        reacting_atoms : List[str]
            a list with chemical symbols of all atoms taking part in
            the reaction

        Returns
        -------
        reacting_atom_indicies: Dict[str, int]
            a dict with reacting atoms symbol as a key and index as a value

        Raises
        ------
        NotImplementedError
            when more than 2 atoms is taking part in the reaction

        '''
        symbol = str(ts_guess_el.symbols)
        reacting_atom_indicies = {}
        visited_species = []
        atom_list = [atom.symbol for atom in ts_guess_el]
        for idx, species in enumerate(reacting_atoms):
            # deal with the edge case when there are two atoms the same typ
            # (eg. 'CO2', 'O2')
            add = 0
            if species not in visited_species:
                visited_species.append(species)
            else:
                add += visited_species.count(species) - 1

            key = '{}_{}'.format(species, idx)
            is_repeted_atom = TSEstimator.is_double_atom(atom_list, species)
            n_atom_list = len(atom_list)

            # deal with the edge case when reactions is not diatomic and there
            # is a repeted atom that does not take part in the reaction,
            # eg. 'CO2' -> in that case increase add by one again
            if is_repeted_atom and n_atom_list > 2:
                add += 1
            reacting_atom_indicies[key] = symbol.find(species) + add
            visited_species.append(species)

        if len(reacting_atom_indicies) > 2:
            raise NotImplementedError('Only two atoms can take part in '
                                      'reaction')
        return reacting_atom_indicies

    @staticmethod
    def is_double_atom(
            atom_list: List[str],
            species: str) -> bool:
        ''' Check if there are two atoms of the same type chemical symbol in
            ts_guess_el

        Parameters
        ----------
        atom_list : list(str)
            a list of atoms in ts_guess_el
        species : str
            a reacting species, e.g. 'OH', 'CO', 'CO2'

        Returns
        -------
        bool
            True if there are more then one atoms of the same chemical symbol.
            False otherwise.

        '''
        count = Counter(atom_list)
        if count[species] > 1:
            return True
        return False


class Diatomic(TSEstimator):
    def get_ts_guess_and_bonded_idx(
            self,
            ts_est: str,
            rxn: Dict[str, str],
            reacting_sp: str,
            reacting_atoms: List[str],
            scfactor: float) -> Tuple[Gratoms, int]:
        ''' Get a ts_guess Gratom object that is ready to be placed on the
            surface and index of atom connecting it to the surface

        Parameters
        ----------
        ts_est : str
            a string representing species that will be used to get ts_guess
            e.g. 'OH', 'COH'
        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        reacting_sp : str
            a key to rxn dictionary
            'reactant' or 'product' are the only avaiable options options
        reacting_atom_indicies: Dict[str, int]
            a dict with reacting atoms symbol as a key and index as a value
        scfactor : float
            a scaling factor used to scale bond distance to get ts_guess

        Returns
        -------
        ts_guess, surface_bonded_idx = Tuple[Gratoms, int]
            a tuple with final ts_guess ready to be placed on the surface
            and index of atom that connects ts_guess with the surface

        '''
        ts_guess_list = Diatomic.build_ts_guess(ts_est)
        # For diatimics, there is only one possible topology, so
        ts_guess_el = ts_guess_list[0]
        return self.rotate_and_scale(ts_guess_el, rxn, reacting_sp,
                                     reacting_atoms, scfactor)

    def rotate_and_scale(
            self,
            ts_guess_el: Gratoms,
            rxn: Dict[str, str],
            reacting_sp: str,
            reacting_atoms: List[str],
            scfactor: float) -> Tuple[Gratoms, int]:
        ''' Rotate ts_guess and scale the bond distance between reacting atoms
            to get a right structures to be placed on the surface

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
        reacting_atom_indicies: Dict[str, int]
            a dict with reacting atoms symbol as a key and index as a value
        scfactor : float
            a scaling factor used to scale bond distance to get ts_guess

        Returns
        -------
        ts_guess, surface_bonded_idx = Tuple[Gratoms, int]
            a tuple with final ts_guess ready to be placed on the surface
            and index of atom that connects ts_guess with the surface

        '''
        surface_bonded_atom_idx = self.get_surface_bonded_atom_idx(
            rxn, reacting_sp)
        reacting_atom_indicies = Diatomic.get_reacting_atoms_indices(
            ts_guess_el, reacting_atoms)

        react_ind_1, react_ind_2 = reacting_atom_indicies.values()

        bondlen = ts_guess_el.get_distance(react_ind_1, react_ind_2)
        ts_guess_el.rotate(90, 'y')
        ts_guess_el.set_distance(
            react_ind_1, react_ind_2, bondlen * scfactor, fix=0)
        return ts_guess_el, surface_bonded_atom_idx


class Triatomic(TSEstimator):
    def get_ts_guess_and_bonded_idx(
            self,
            ts_est: str,
            rxn: Dict[str, str],
            reacting_sp: str,
            reacting_atoms: List[str],
            scfactor: float,
            conf: str = None) -> Tuple[Gratoms, int]:
        '''[summary]

        Parameters
        ----------
        ts_est : str
            a string representing species that will be used to get ts_guess
            e.g. 'OH', 'COH'
        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        reacting_sp : str
            a key to rxn dictionary
            'reactant' or 'product' are the only avaiable options options
        reacting_atom_indicies: Dict[str, int]
            a dict with reacting atoms symbol as a key and index as a value
        scfactor : float
            a scaling factor used to scale bond distance to get ts_guess
        conf : str, optional
            specify which topology to choose to construct ts_guess, e.g 'sp'
            or None (None -> 'sp2')
            by default None

        Returns
        -------
        ts_guess, surface_bonded_idx = Tuple[Gratoms, int]
            a tuple with final ts_guess ready to be placed on the surface
            and index of atom that connects ts_guess with the surface

        '''
        ts_guess_list = Triatomic.build_ts_guess(ts_est)
        if conf == 'sp':
            ts_guess_el = ts_guess_list[1]
        else:
            ts_guess_el = ts_guess_list[0]

        return self.rotate_and_scale(ts_guess_el, rxn, reacting_sp,
                                     reacting_atoms, scfactor)

    def rotate_and_scale(
            self,
            ts_guess_el: Gratoms,
            rxn: Dict[str, str],
            reacting_sp: str,
            reacting_atoms: List[str],
            scfactor: float) -> Tuple[Gratoms, int]:
        ''' Rotate ts_guess and scale the bond distance between reacting atoms
            to get a right structures to be placed on the surface

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
        reacting_atom_indicies: Dict[str, int]
            a dict with reacting atoms symbol as a key and index as a value
        scfactor : float
            a scaling factor used to scale bond distance to get ts_guess

        Returns
        -------
        ts_guess, surface_bonded_idx = Tuple[Gratoms, int]
            a tuple with final ts_guess ready to be placed on the surface
            and index of atom that connects ts_guess with the surface

        '''
        surface_bonded_atom_idx = self.get_surface_bonded_atom_idx(
            rxn, reacting_sp)
        reacting_atom_indicies = Triatomic.get_reacting_atoms_indices(
            ts_guess_el, reacting_atoms)

        react_ind_1, react_ind_2 = reacting_atom_indicies.values()

        bondlen = ts_guess_el.get_distance(react_ind_1, react_ind_2)
        ts_guess_el.rotate(90, 'z')
        ts_guess_el.set_distance(
            react_ind_1, react_ind_2, bondlen * scfactor, fix=0)
        # hardcoded values based on empirical tests
        ts_guess_el.set_angle(
            0, 1, 2, -30, indices=[0, 1, 2], add=True)
        return ts_guess_el, surface_bonded_atom_idx


class Tetraatomic(TSEstimator):
    def get_ts_guess_and_bonded_idx(
            self,
            ts_est: str,
            rxn: Dict[str, str],
            reacting_sp: str,
            reacting_atoms: List[str],
            scfactor: float,
            conf: str = None) -> Tuple[Gratoms, int]:
        '''[summary]

        Parameters
        ----------
        ts_est : str
            a string representing species that will be used to get ts_guess
            e.g. 'OH', 'COH'
        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        reacting_sp : str
            a key to rxn dictionary
            'reactant' or 'product' are the only avaiable options options
        reacting_atom_indicies: Dict[str, int]
            a dict with reacting atoms symbol as a key and index as a value
        scfactor : float
            a scaling factor used to scale bond distance to get ts_guess
        conf : str, optional
            specify which topology to choose to construct ts_guess, e.g 'sp'
            or None (None -> 'sp2')
            by default None

        Returns
        -------
        ts_guess, surface_bonded_idx = Tuple[Gratoms, int]
            a tuple with final ts_guess ready to be placed on the surface
            and index of atom that connects ts_guess with the surface

        '''
        ts_guess_list = Tetraatomic.build_ts_guess(ts_est)
        if conf == 'sp':
            ts_guess_el = ts_guess_list[1]
        else:
            ts_guess_el = ts_guess_list[0]

        return self.rotate_and_scale(ts_guess_el, rxn, reacting_sp,
                                     reacting_atoms, scfactor)

    def rotate_and_scale(
            self,
            ts_guess_el: Gratoms,
            rxn: Dict[str, str],
            reacting_sp: str,
            reacting_atoms: List[str],
            scfactor: float) -> Tuple[Gratoms, int]:
        ''' Rotate ts_guess and scale the bond distance between reacting atoms
            to get a right structures to be placed on the surface

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
        reacting_atom_indicies: Dict[str, int]
            a dict with reacting atoms symbol as a key and index as a value
        scfactor : float
            a scaling factor used to scale bond distance to get ts_guess

        Returns
        -------
        ts_guess, surface_bonded_idx = Tuple[Gratoms, int]
            a tuple with final ts_guess ready to be placed on the surface
            and index of atom that connects ts_guess with the surface

        '''
        surface_bonded_atom_idx = self.get_surface_bonded_atom_idx(
            rxn, reacting_sp)
        reacting_atom_indicies = Triatomic.get_reacting_atoms_indices(
            ts_guess_el, reacting_atoms)

        react_ind_1, react_ind_2 = reacting_atom_indicies.values()

        bondlen = ts_guess_el.get_distance(react_ind_1, react_ind_2)
        ts_guess_el.rotate(150, 'z')
        ts_guess_el.set_distance(
            react_ind_1, react_ind_2, bondlen * scfactor, fix=0)
        # hardcoded values based on empirical tests
        ts_guess_el.set_angle(
            0, 1, 2, -30, indices=[0, 1, 2, 3, 4], add=True)
        return ts_guess_el, surface_bonded_atom_idx
