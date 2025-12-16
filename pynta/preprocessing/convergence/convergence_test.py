
from ase.calculators.espresso import Espresso
from ase.constraints import FixAtoms

def freeze_bottom_half(atoms):
    """Fix atoms with tag > max_tag (e.g. bottom layers)."""
    z = atoms.positions[:, 2]
    z_mid = 0.5 * (z.max() + z.min())
    mask = z < z_mid
    atoms.set_constraint(FixAtoms(mask=mask))
    return atoms

def run_convergence(
    atoms,
    base_input,
    pseudopotentials,
    profile,
    ecut_values,
    kmesh_values
):
    """Run ecut and k-point convergence tests.

    Returns:
        results: list of (ecut, kpts, energy)
    """
    results = []

    for ecut in ecut_values:
        for kpts in kmesh_values:
            input_data = base_input.copy()
            input_data['ecutwfc'] = ecut
            input_data['ecutrho'] = 4 * ecut

            calc = Espresso(
                input_data=input_data,
                pseudopotentials=pseudopotentials,
                kpts=kpts,
                profile=profile
            )

            atoms.set_calculator(calc)
            E = atoms.get_potential_energy()
            print(f"ecut={ecut} Ry, kpts={kpts}, Energy={E:.6f} eV")
            results.append((ecut, kpts, E))

    return results
