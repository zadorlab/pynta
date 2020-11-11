from rmgcat_to_sella.excatkit.molecule import Molecule
from rmgcat_to_sella.excatkit.gratoms import Gratoms
from typing import List, Dict, Tuple
import os
from ase.io import read, write


class Diatomic():
    def get_ts_guess_and_bonded_idx(
            self,
            ts_est: str,
            rxn: Dict[str, str],
            reacting_sp: str,
            scfactor: float) -> Tuple[Gratoms, int]:
        ts_guess_list = Diatomic.build_ts_guess(ts_est)
        ts_guess_el = ts_guess_list[0]
        ts_guess, surface_bonded_atom_idx = self.rotate_and_scale(
            ts_guess_el, rxn, reacting_sp, scfactor)
        return ts_guess, surface_bonded_atom_idx

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

    def deal_with_bonds(
            self,
            ts_guess_el,
            rxn,
            reacting_sp):
        surface_bonded_atoms = self.get_surface_bonded_atoms(
            rxn, reacting_sp)
        if len(surface_bonded_atoms) > 1:
            raise NotImplementedError(
                'Only one atom can be connectedto the surface. '
                'Support for many atoms will be added later.')
        else:
            bonded = surface_bonded_atoms[0]
            surface_bonded_atom_idx = self.get_bonded_index(
                bonded, ts_guess_el)
        return surface_bonded_atom_idx

    @staticmethod
    def get_surface_bonded_atoms(
            rxn,
            reacting_sp):
        atomic_connections = Diatomic.get_atomic_connections(rxn, reacting_sp)
        surface_bonded_atoms = []
        for k, v in atomic_connections.items():
            if v == max(atomic_connections.values()):
                surface_bonded_atoms.append(k)
        return surface_bonded_atoms

    @staticmethod
    def get_bonded_index(
            bonded,
            ts_guess_el):
        symbol = str(ts_guess_el.symbols)
        surface_bonded_atom_idx = symbol.find(bonded)
        return surface_bonded_atom_idx

    @staticmethod
    def get_atomic_connections(
            rxn,
            reacting_sp):
        atomic_connections = {}
        reacting_sp_connectivity = rxn[reacting_sp].split('\n')
        for line in reacting_sp_connectivity:
            if '*' in line:
                connections = line.count('{')
                symbol = line.split()[2]
                atomic_connections[symbol] = connections
        return atomic_connections

    def rotate_and_scale(
            self,
            ts_guess_el,
            rxn,
            reacting_sp,
            scfactor):
        surface_bonded_atom_idx = self.deal_with_bonds(
            ts_guess_el, rxn, reacting_sp)
        if surface_bonded_atom_idx == 0:
            other_atom_idx = 1
        else:
            other_atom_idx = 0
        bondlen = ts_guess_el.get_distance(
            surface_bonded_atom_idx, other_atom_idx)
        ts_guess_el.rotate(90, 'y')
        ts_guess_el.set_distance(
            surface_bonded_atom_idx, other_atom_idx, bondlen * scfactor, fix=0)
        return ts_guess_el, surface_bonded_atom_idx


class Triatomic(Diatomic):
    def get_ts_guess_and_bonded_idx(
            self,
            ts_est,
            rxn,
            reacting_sp,
            reacting_atoms,
            scfactor,
            conf=None):
        ts_guess_list = Triatomic.build_ts_guess(ts_est)
        if conf == 'sp':
            ts_guess_el = ts_guess_list[1]
        else:
            ts_guess_el = ts_guess_list[0]

        ts_guess, surface_bonded_atom_idx = self.rotate_and_scale_tri(
            ts_guess_el, rxn, reacting_sp, reacting_atoms, scfactor)
        return ts_guess, surface_bonded_atom_idx

    @staticmethod
    def get_reacting_atoms_indices(
            ts_guess_el,
            reacting_atoms):
        symbol = str(ts_guess_el.symbols)
        reacting_atom_idx = {}
        for species in reacting_atoms:
            reacting_atom_idx[species] = symbol.find(species)
        return reacting_atom_idx

    def rotate_and_scale_tri(
            self,
            ts_guess_el,
            rxn,
            reacting_sp,
            reacting_atoms,
            scfactor):
        surface_bonded_atom_idx = self.deal_with_bonds(
            ts_guess_el, rxn, reacting_sp)
        reacting_atom_idx = Triatomic.get_reacting_atoms_indices(
            ts_guess_el, reacting_atoms)

        if len(reacting_atom_idx) > 2:
            raise NotADirectoryError('Only two atoms can take part in '
                                     'reaction')

        react_ind_1, react_ind_2 = reacting_atom_idx.values()

        bondlen = ts_guess_el.get_distance(react_ind_1, react_ind_2)
        # ts_guess_el.rotate(90, 'y')
        ts_guess_el.set_distance(
            react_ind_1, react_ind_2, bondlen * scfactor, fix=0)
        return ts_guess_el, surface_bonded_atom_idx
