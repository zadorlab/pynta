from xtb.ase.calculator import XTB
import numpy as np
import ase
from ase.atoms import Atoms
from ase.units import Hartree, Bohr
from ase.geometry import get_distances
from ase.calculators import calculator
from ase.data import reference_states,chemical_symbols
from ase.build import bulk
from ase.constraints import *
from pynta.utils import name_to_ase_software
from sella import Sella, Constraints
import scipy.optimize as opt
import copy
from copy import deepcopy
import itertools
from pynta.utils import *

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

def run_harmonically_forced_xtb(atoms,atom_bond_potentials,site_bond_potentials,nslab,
        molecule_to_atom_maps=None,ase_to_mol_num=None,method="GFN1-xTB",constraints=[]):
    """
    Optimize TS guess using xTB + harmonic forcing terms determined by atom_bond_potentials and site_bond_potentials
    """
    out_constraints = []
    for c in constraints:
        if isinstance(c,dict):
            constraint = construct_constraint(c)
            out_constraints.append(constraint)
        elif c == "freeze half slab":
            out_constraints.append(FixAtoms([
                atom.index for atom in atoms if atom.position[2] < atoms.cell[2, 2] / 2.
            ]))
        elif c == "freeze slab":
            out_constraints.append(FixAtoms(
                indices=list(range(nslab))
                ))
        elif c.split()[0] == "freeze" and c.split()[1] == "all": #ex: "freeze all Cu"
            sym = c.split()[2]
            out_constraints.append(FixAtoms(
                indices=[atom.index for atom in atoms if atom.symbol == sym]
                ))
        elif c.split()[0] == "freeze" and c.split()[1] == "up" and c.split()[2] == "to":
            n = int(c.split()[3])
            out_constraints.append(FixAtoms(
                indices=list(range(n))
                ))
            
    atoms.set_constraint(out_constraints)
    
    hfxtb = HarmonicallyForcedXTB(method="GFN1-xTB",
                              atom_bond_potentials=atom_bond_potentials,
                             site_bond_potentials=site_bond_potentials)

    atoms.calc = hfxtb

    opt = Sella(atoms,trajectory="xtbharm.traj",order=0)

    try:
        opt.run(fmax=0.02,steps=150)
    except Exception as e: #no pbc fallback
        return run_harmonically_forced_xtb_no_pbc(atoms,atom_bond_potentials,site_bond_potentials,nslab,
                                       molecule_to_atom_maps=molecule_to_atom_maps,ase_to_mol_num=ase_to_mol_num,
                                               constraints=constraints,method=method,dthresh=4.0)

    Eharm,Fharm = atoms.calc.get_energy_forces()

    return atoms,Eharm,Fharm

def run_harmonically_forced_xtb_no_pbc(atoms,atom_bond_potentials,site_bond_potentials,nslab,
                               molecule_to_atom_maps,ase_to_mol_num=None,
                                       constraints=[],method="GFN1-xTB",dthresh=4.0):
    """
    This algorithm extends the slab and selects a section of the slab in which the adsorbates are closest
    together, truncates the slab around it and optimizes without pbc based on the truncated slab before
    translating the result back to the original periodic slab

    dthresh is the threshold x-y distance a slab atom must be away from an adsorbate atom to be truncated
    """
    if ase_to_mol_num is None: #assume only one adsorbate and translate that as a unit
        ase_to_mol_num = {i+nslab:0 for i in range(len(atoms) - nslab)}

    if not isinstance(molecule_to_atom_maps[0],list):
        molecule_to_atom_maps = [molecule_to_atom_maps]

    target = atoms.cell[0][:2] + atoms.cell[1][:2]

    og_pbc = atoms.pbc # Need to apply original pbc's on atoms object (depends on QC software)

    site_poss = []
    site_inds = []
    site_mols = []
    for site_bond_potential in site_bond_potentials:
        ind = site_bond_potential["ind"]
        site_pos = site_bond_potential["site_pos"]
        mol_ind = ase_to_mol_num[ind]

        if mol_ind not in site_mols:
            site_poss.append(site_pos)
            site_inds.append(ind)
            site_mols.append(mol_ind)
        else:
            i = site_mols.index(mol_ind)
            if np.linalg.norm(target-site_pos[:2]) > np.linalg.norm(target-site_poss[i][:2]):
                site_poss[i] = site_pos
                site_inds[i] = ind

    for mol_ind in set(ase_to_mol_num.values()): #handle gas phase species that don't have site_bond_potentials
        if mol_ind not in site_mols:
            for key,val in ase_to_mol_num.items():
                if val == mol_ind:
                    site_mols.append(mol_ind)
                    site_poss.append(atoms.positions[key])
                    break

    aposx = [a.position[0] for a in atoms[nslab:]] + [a[0] for a in site_poss]
    aposy = [a.position[1] for a in atoms[nslab:]] + [a[1] for a in site_poss]
    apos = [np.array([min(aposx),min(aposy)]),
            np.array([min(aposx),max(aposy)]),
            np.array([max(aposx),min(aposy)]),
            np.array([max(aposx),max(aposy)])]
    trans = get_best_translation(site_poss,apos,atoms.cell)

    mol_to_trans = {site_mols[i]: trans[i] for i in range(len(site_mols))}

    new_site_potentials = []
    for i,site_potential in enumerate(site_bond_potentials):
        sitep = deepcopy(site_potential)
        mol_ind = ase_to_mol_num[site_potential["ind"]]
        sitep["site_pos"] = (np.array(sitep["site_pos"]) + mol_to_trans[mol_ind]).tolist()
        new_site_potentials.append(sitep)


    slab = atoms[:nslab]
    ad = atoms[nslab:]
    bigslab = slab * (3,3,1)

    for ind,molind in ase_to_mol_num.items():
        ad.positions[ind-nslab] += mol_to_trans[molind]


    bigad = bigslab+ad

    delinds = []
    for i,atom in enumerate(bigslab):
        dist = np.inf
        for j,a in enumerate(ad):
            if np.linalg.norm(atom.position[:2]-a.position[:2]) < dist:
                dist = np.linalg.norm(atom.position[:2]-a.position[:2])
        for sitep in new_site_potentials:
            if np.linalg.norm(atom.position[:2]-sitep["site_pos"][:2]) < dist:
                dist = np.linalg.norm(atom.position[:2]-sitep["site_pos"][:2])
        if dist > dthresh:
            delinds.append(i)

    for ind in reversed(sorted(delinds)):
        del bigad[ind]

    bigad.pbc = (False,False,False)
    new_nslab = len(bigad) - len(ad)

    for site_potential in new_site_potentials:
        site_potential["ind"] += new_nslab-nslab

    new_atom_bond_potentials = []
    for atom_bond_potential in atom_bond_potentials:
        abpot = deepcopy(atom_bond_potential)
        abpot["ind1"] += new_nslab-nslab
        abpot["ind2"] += new_nslab-nslab
        new_atom_bond_potentials.append(abpot)

    new_constraints = []
    for constraint in constraints:
        if constraint == "freeze slab":
            new_constraints.append(constraint)
            continue
        else:
            c = deepcopy(constraint)
            if "indices" in c:
                for i in range(len(c["indices"])):
                    c["indices"][i] += new_nslab-nslab
            if "a1" in c:
                c["a1"] += new_nslab-nslab
            if "a2" in c:
                c["a2"] += new_nslab-nslab
            new_constraints.append(c)
    
    out_constraints = []
    for c in new_constraints:
        if isinstance(c,dict):
            constraint = construct_constraint(c)
            out_constraints.append(constraint)
        elif c == "freeze half slab":
            out_constraints.append(FixAtoms([
                atom.index for atom in bigad if atom.position[2] < bigad.cell[2, 2] / 2.
            ]))
        elif c == "freeze slab":
            out_constraints.append(FixAtoms(
                indices=list(range(new_nslab))
                ))
        elif c.split()[0] == "freeze" and c.split()[1] == "all": #ex: "freeze all Cu"
            sym = c.split()[2]
            out_constraints.append(FixAtoms(
                indices=[atom.index for atom in bigad if atom.symbol == sym]
                ))
        elif c.split()[0] == "freeze" and c.split()[1] == "up" and c.split()[2] == "to":
            n = int(c.split()[3])
            out_constraints.append(FixAtoms(
                indices=list(range(n))
                ))
    
    hfxtb = HarmonicallyForcedXTB(method="GFN1-xTB",
                                  atom_bond_potentials=new_atom_bond_potentials,
                                 site_bond_potentials=new_site_potentials)
    bigad.set_constraint(out_constraints)
    bigad.calc = hfxtb

    opt = Sella(bigad,trajectory="xtbharm.traj",order=0)

    try:
        opt.run(fmax=0.02,steps=150)
    except:
        return None,None,None

    Eharm,Fharm = bigad.calc.get_energy_forces()

    newad = bigad[new_nslab:]
    for ind,molind in ase_to_mol_num.items():
        newad.positions[ind-nslab] -= mol_to_trans[molind]

    outadslab = slab + newad

    outadslab.pbc = og_pbc

    return outadslab,Eharm,Fharm

def get_best_translation(poss,apos,cell):
    target = (cell[0][:2] + cell[1][:2])*1.5
    pos2ds = [np.array(pos[:2]) for pos in poss]
    translations = [np.zeros(2),cell[0][:2],cell[1][:2],cell[0][:2] + cell[1][:2],
                2*cell[0][:2], 2*cell[1][:2], 2*cell[0][:2] + cell[1][:2],cell[0][:2] + 2*cell[1][:2],2*cell[0][:2] + 2*cell[1][:2]]
    mindist = np.inf
    minpos2ddist = np.inf
    tranout = None
    for transs in itertools.product(translations,repeat=len(pos2ds)):
        dist = 0.0
        pos2ddist = 0.0
        for i in range(len(pos2ds)):
            for pos in apos:
                dist += np.linalg.norm(pos+transs-target)
            for j in range(i):
                pos2ddist += np.linalg.norm(pos2ds[i]+transs[i]-(pos2ds[j]+transs[j]))
        if abs(pos2ddist - minpos2ddist) > 0.1 and pos2ddist < minpos2ddist:
            tranout = transs
            mindist = dist
            minpos2ddist = pos2ddist
        elif abs(pos2ddist - minpos2ddist) < 0.1 and dist < mindist:
            tranout = transs
            mindist = dist
            minpos2ddist = pos2ddist

    return [np.array([x[0],x[1],0.0]) for x in tranout]

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

def get_lattice_parameter(metal,surface_type,software,software_kwargs,da=0.1,options={"xatol":1e-4},a0=None):
    soft = name_to_ase_software(software)(**software_kwargs)
    def f(a):
        slab = bulk(metal,surface_type[:3],a=a)
        slab.calc = soft
        slab.pbc = (True, True, True)
        return slab.get_potential_energy()
    if a0 is None:
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

def optimize_lattice_parameter(metal,surface_type,software,software_kwargs,da=0.1,options={"xatol":1e-4},a0=None):
    soft = name_to_ase_software(software)(**software_kwargs)
    def f(a):
        slab = bulk(metal, surface_type[:3], a=a, orthorhombic=True)
        slab.calc = soft
        slab.pbc = (True, True, True)
        return slab.get_potential_energy()

    # Get the initial lattice parameter a0
    a0 = reference_states[chemical_symbols.index(metal)]['a']
    
    # Initialize the difference to a large number
    diff = float('inf')
    
    # Loop until the difference between a0 and a is less than or equal to 0.001
    iteration = 0
    while diff > 0.001:
        # Generate a range of a values
        avals = np.arange(a0 - da, a0 + da, 0.01)
        outavals = []
        Evals = []
        
        print(f"Iteration #{iteration}: a0 = {a0}")
        
        # Calculate potential energies for the range of a values
        for a in avals:
            try:
                E = f(a)
                outavals.append(a)
                Evals.append(E)
                print((a, E))
            except Exception as e:
                print(f"Error calculating energy for a = {a}: {e}")
        
        # Print a values and E values
        print("a values:")
        print(outavals)
        print("E values:")
        print(Evals)
        
        # Find the indices of the 7 smallest Evals
        inds = np.argsort(np.array(Evals))[:7]
        
        # Fit a polynomial to the selected points and calculate the interpolated a
        p = np.polyfit(np.array(outavals)[inds], np.array(Evals)[inds], 2)
        a = -p[1] / (2.0 * p[0])
        
        # Print the initial and interpolated values
        print(f"Initial a0: {a0}, Interpolated a: {a}")
        
        # Save the current a as a0 for the next iteration
        previous_a0 = a0
        a0 = a
        
        # Minimize the function within the bounds
        out = opt.minimize_scalar(f, method='bounded', bounds=(a - 0.01, a + 0.01), options=options)
        
        # Print the optimization result
        print(out)
        print("Optimized a: {}".format(out.x))
        
        # Calculate the difference between the new a and the previous a
        diff = abs(out.x - previous_a0)
        
        # Print the results of this iteration
        print(f"#{iteration + 1}: a0 = {previous_a0}, a = {out.x}")
        
        # Increment the iteration counter
        iteration += 1

    print("Converged!")
    print(f"Final optimized a: {out.x}")
    return out.x
