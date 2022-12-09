from xtb.ase.calculator import XTB
import numpy as np
import ase
from ase.atoms import Atoms
from ase.units import Hartree, Bohr
from ase.geometry import get_distances
from ase.calculators import calculator
from ase.data import reference_states,chemical_symbols
from ase.build import bulk
from pynta.utils import name_to_ase_software
from sella import Sella, Constraints
import scipy.optimize as opt
import copy

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

def run_harmonically_forced_xtb(atoms,atom_bond_potentials,site_bond_potentials,nslab,method="GFN1-xTB",
                               constraints=[]):
    """
    Optimize TS guess using xTB + harmonic forcing terms determined by atom_bond_potentials and site_bond_potentials
    """
    cons = Constraints(atoms)

    for c in constraints:
        if isinstance(c,dict):
            add_sella_constraint(cons,c)
        elif c == "freeze slab":
            for i,atom in enumerate(atoms): #freeze the slab
                if i < nslab:
                    cons.fix_translation(atom.index)
        else:
            raise ValueError("Constraint {} not understood".format(c))


    hfxtb = HarmonicallyForcedXTB(method="GFN1-xTB",
                              atom_bond_potentials=atom_bond_potentials,
                             site_bond_potentials=site_bond_potentials)

    atoms.calc = hfxtb

    opt = Sella(atoms,constraints=cons,trajectory="xtbharm.traj",order=0)

    try:
        opt.run(fmax=0.02,steps=150)
    except Exception as e:
        return None,None,None

    Eharm,Fharm = atoms.calc.get_energy_forces()

    return atoms,Eharm,Fharm

def add_sella_constraint(cons,d):
    """
    construct a constraint from a dictionary that is the input to the constraint
    constructor plus an additional "type" key that indices the name of the constraint
    in this case for Sella the full Constraints object cons must be included in the inputs
    adds the constraint to Constraints and returns None
    """
    constraint_dict = copy.deepcopy(d)
    constructor = getattr(cons,constraint_dict["type"])
    del constraint_dict["type"]
    constructor(**constraint_dict)
    return

def get_lattice_parameter(metal,surface_type,software,software_kwargs,da=0.1,options={"xatol":1e-4}):
    soft = name_to_ase_software(software)(**software_kwargs)
    def f(a):
        slab = bulk(metal,surface_type[:3],a=a)
        slab.calc = soft
        slab.pbc = (True, True, True)
        return slab.get_potential_energy()

    a0 = reference_states[chemical_symbols.index(metal)]['a']
    avals = np.arange(a0-da,a0+da,0.01)
    outavals = []
    Evals = []
    print("a,E")
    for a in avals:
        try:
            E = f(a)
            outavals.append(a)
            Evals.append(E)
            print((a,E))
        except:
            pass
    print("a values:")
    print(outavals)
    print("E values:")
    print(Evals)
    inds = np.argsort(np.array(Evals))[:7]
    p = np.polyfit(np.array(outavals)[inds],np.array(Evals)[inds],2)
    a = -p[1]/(2.0*p[0])
    print("ASE reference a: {}".format(a0))
    print("Interpolated a: {}".format(a))
    out = opt.minimize_scalar(f,method='bounded',bounds=(a-0.01,a+0.01),options=options)
    print(out)
    print("Optimized a: {}".format(out.x))
    return out.x
