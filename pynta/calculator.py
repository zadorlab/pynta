from xtb.ase.calculator import XTB
import numpy as np
from ase.atoms import Atoms
from ase.units import Hartree, Bohr
from ase.geometry import get_distances
from ase.calculators import calculator

def get_energy_atom_bond(atoms,ind1,ind2,k,deq):
    bd,d = get_distances([atoms.positions[ind1]], [atoms.positions[ind2]], cell=atoms.cell, pbc=atoms.pbc)
    return k*(d-deq)**2

def get_forces_atom_bond(atoms,ind1,ind2,k,deq):
    forces = np.zeros(atoms.positions.shape)
    bd,d = get_distances([atoms.positions[ind1]], [atoms.positions[ind2]], cell=atoms.cell, pbc=atoms.pbc)
    if d != 0.0:
        forces[ind1] = 2.0*bd*(1.0-deq/d)
        forces[ind2] = -forces[ind1]
    else:
        forces[ind1] = bd
        forces[ind2] = bd
    return k*forces

def get_energy_forces_atom_bond(atoms,ind1,ind2,k,deq):
    forces = np.zeros(atoms.positions.shape)
    bd,d = get_distances([atoms.positions[ind1]], [atoms.positions[ind2]], cell=atoms.cell, pbc=atoms.pbc)
    if d != 0.0:
        forces[ind1] = 2.0*bd*(1.0-deq/d)
        forces[ind2] = -forces[ind1]
    else:
        forces[ind1] = bd
        forces[ind2] = bd
    energy = k*(d-deq)**2
    return energy,k*forces

def get_energy_site_bond(atoms,ind,site_pos,k,deq):
    bd,d = get_distances([atoms.positions[ind]], [site_pos], cell=atoms.cell, pbc=atoms.pbc)
    return k*(d-deq)**2

def get_forces_site_bond(atoms,ind,site_pos,k,deq):
    forces = np.zeros(atoms.positions.shape)
    bd,d = get_distances([atoms.positions[ind]], [site_pos], cell=atoms.cell, pbc=atoms.pbc)
    if d != 0.0:
        forces[ind] = 2.0*bd*(1.0-deq/d)
    else:
        forces[ind] = bd
    return k*forces

def get_energy_forces_site_bond(atoms,ind,site_pos,k,deq):
    forces = np.zeros(atoms.positions.shape)
    bd,d = get_distances([atoms.positions[ind]], [site_pos], cell=atoms.cell, pbc=atoms.pbc)
    if d != 0:
        forces[ind] = 2.0*bd*(1.0-deq/d)
    else:
        forces[ind] = bd
    energy = k*(d-deq)**2
    return energy,k*forces

class HarmonicallyForcedXTB(XTB):
    def get_energy_forces(self):
        energy = 0.0
        forces = np.zeros(self.atoms.positions.shape)
        if hasattr(self.parameters,"atom_bond_potentials"):
            for atom_bond_potential in self.parameters.atom_bond_potentials:
                E,F = get_energy_forces_atom_bond(self.atoms,**atom_bond_potential)
                energy += E
                forces += F

        if hasattr(self.parameters,"site_bond_potentials"):
            for site_bond_potential in self.parameters.site_bond_potentials:
                E,F = get_energy_forces_site_bond(self.atoms,**site_bond_potential)
                energy += E
                forces += F
        
        return energy[0][0],forces

    def calculate(self, atoms=None, properties=None, system_changes=calculator.all_changes):
        XTB.calculate(self,atoms=atoms,properties=properties,system_changes=system_changes)
        energy,forces = self.get_energy_forces()
        self.results["energy"] += energy
        self.results["free_energy"] += energy
        self.results["forces"] += forces
