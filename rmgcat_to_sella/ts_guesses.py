from rmgcat_to_sella.excatkit.molecule import Molecule


class Diatomic():
    def get_ts_candidate(self, ts_guess, rxn, reacting_sp):
        ts_candidate = Diatomic.build_ts_candidate(ts_guess)

        get_surface_bonded_atoms = self.get_surface_bonded_atoms(
            rxn, reacting_sp)
        if len(get_surface_bonded_atoms) > 1:
            raise NotImplementedError(
                'Only one atom can be connectedto the surface. '
                'Support for many atoms will be added later.')
        else:
            bonded = get_surface_bonded_atoms[0]
            bonded_idx = self.get_bonded_index(bonded, ts_candidate)
        return ts_candidate, bonded_idx

    @staticmethod
    def build_ts_candidate(ts_guess):
        ts_candidate = Molecule().molecule(ts_guess)
        return ts_candidate

    def get_surface_bonded_atoms(self, rxn, reacting_sp):
        atomic_connections = self.get_atomic_connections(rxn, reacting_sp)
        surface_bonded_atoms = []
        for k, v in atomic_connections.items():
            if v == max(atomic_connections.values()):
                surface_bonded_atoms.append(k)
        return surface_bonded_atoms

    def get_bonded_index(self, bonded, ts_candidate):
        symbol = str(ts_candidate[0].symbols)
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

    def rotate_ts_candidate(self):
        pass
