from rmgcat_to_sella.excatkit.molecule import Molecule


class Diatomic():
    def get_ts_guess_and_bonded_idx(
            self,
            ts_guess,
            rxn,
            reacting_sp,
            scfactor):
        ts_guess_init = Diatomic.build_ts_guess(ts_guess)
        ts_guess, bonded_idx = self.rotate_and_scale(
            ts_guess_init, rxn, reacting_sp, scfactor)
        return ts_guess, bonded_idx

    @staticmethod
    def build_ts_guess(ts_guess):
        ts_guess = Molecule().molecule(ts_guess)
        return ts_guess

    def deal_with_bonds(self, ts_guess_init, rxn, reacting_sp):
        surface_bonded_atoms = self.get_surface_bonded_atoms(
            rxn, reacting_sp)
        if len(surface_bonded_atoms) > 1:
            raise NotImplementedError(
                'Only one atom can be connectedto the surface. '
                'Support for many atoms will be added later.')
        else:
            bonded = surface_bonded_atoms[0]
            bonded_idx = self.get_bonded_index(bonded, ts_guess_init)
        return bonded_idx

    def get_surface_bonded_atoms(self, rxn, reacting_sp):
        atomic_connections = self.get_atomic_connections(rxn, reacting_sp)
        surface_bonded_atoms = []
        for k, v in atomic_connections.items():
            if v == max(atomic_connections.values()):
                surface_bonded_atoms.append(k)
        return surface_bonded_atoms

    @staticmethod
    def get_bonded_index(bonded, ts_guess_init):
        symbol = str(ts_guess_init[0].symbols)
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

    def rotate_and_scale(self, ts_guess_init, rxn, reacting_sp, scfactor):
        bonded_idx = self.deal_with_bonds(ts_guess_init, rxn, reacting_sp)
        if bonded_idx == 0:
            other_atom_idx = 1
        else:
            other_atom_idx = 0
        bondlen = ts_guess_init[0].get_distance(bonded_idx, other_atom_idx)
        ts_guess_init[0].rotate(90, 'y')
        ts_guess_init[0].set_distance(
            bonded_idx, other_atom_idx, bondlen * scfactor, fix=0)
        return ts_guess_init[0], bonded_idx
