from os.path import normpath
from rmgcat_to_sella.excatkit.molecule import Molecule
from ase import Atoms


class Diatomic():
    def get_ts_candidate(self, ts_guess, rxn, reacting_sp, scfactor):
        ts_candidate_init = Diatomic.build_ts_candidate(ts_guess)
        ts_candidate, bonded_idx = self.rotate_and_scale(
            ts_candidate_init, rxn, reacting_sp, scfactor)
        return ts_candidate, bonded_idx

    @staticmethod
    def build_ts_candidate(ts_guess):
        ts_candidate = Molecule().molecule(ts_guess)
        return ts_candidate

    def deal_with_bonds(self, ts_candidate_init, rxn, reacting_sp):
        surface_bonded_atoms = self.get_surface_bonded_atoms(
            rxn, reacting_sp)
        if len(surface_bonded_atoms) > 1:
            raise NotImplementedError(
                'Only one atom can be connectedto the surface. '
                'Support for many atoms will be added later.')
        else:
            bonded = surface_bonded_atoms[0]
            bonded_idx = self.get_bonded_index(bonded, ts_candidate_init)
        return bonded_idx

    def get_surface_bonded_atoms(self, rxn, reacting_sp):
        atomic_connections = self.get_atomic_connections(rxn, reacting_sp)
        surface_bonded_atoms = []
        for k, v in atomic_connections.items():
            if v == max(atomic_connections.values()):
                surface_bonded_atoms.append(k)
        return surface_bonded_atoms

    @staticmethod
    def get_bonded_index(bonded, ts_candidate_init):
        symbol = str(ts_candidate_init[0].symbols)
        bonded_idx = symbol.find(bonded)
        return bonded_idx

    def get_atomic_connections(self, rxn, reacting_sp):
        atomic_connections = {}
        reacting_sp_connectivity = rxn[reacting_sp].split('\n')
        for line in reacting_sp_connectivity:
            if '*' in line:
                connections = line.count('{')
                symbol = line.split()[2]
                atomic_connections[symbol] = connections
        return atomic_connections

    def rotate_and_scale(self, ts_candidate_init, rxn, reacting_sp, scfactor):
        bonded_idx = self.deal_with_bonds(ts_candidate_init, rxn, reacting_sp)
        if bonded_idx == 0:
            other_atom_idx = 1
        else:
            other_atom_idx = 0
        bondlen = ts_candidate_init[0].get_distance(bonded_idx, other_atom_idx)
        ts_candidate_init[0].rotate(90, 'y')
        ts_candidate_init[0].set_distance(
            bonded_idx, other_atom_idx, bondlen * scfactor, fix=0)
        return ts_candidate_init[0], bonded_idx
