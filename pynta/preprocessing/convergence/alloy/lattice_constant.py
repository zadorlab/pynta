import numpy as np
import scipy.optimize as opt
from math import gcd
from ase.build import bulk
from ase.calculators.espresso import Espresso


def infer_crystal_from_surface(surface_type, crystal=None):
    """
    surface_type examples: 'fcc111', 'fcc100', 'bcc110', 'hcp0001'
    Returns crystal: 'fcc', 'bcc', or 'hcp'
    """
    if crystal is not None:
        return crystal.lower()

    st = surface_type.lower()
    if st.startswith("fcc"):
        return "fcc"
    if st.startswith("bcc"):
        return "bcc"
    if st.startswith("hcp"):
        return "hcp"
    raise ValueError(f"Cannot infer crystal from surface_type='{surface_type}'. Provide crystal explicitly.")


def make_qe_calc(pseudos, kpts=(12, 12, 12), ecutwfc=60, ecutrho=480, degauss=0.02, qe_input_overrides=None):
    input_data = {
        "control": {"calculation": "scf", "tstress": True, "tprnfor": True},
        "system": {
            "ecutwfc": ecutwfc,
            "ecutrho": ecutrho,
            "occupations": "smearing",
            "smearing": "mp",
            "degauss": degauss,
        },
        "electrons": {"conv_thr": 1e-8},
    }
    if qe_input_overrides:
        for nl, d in qe_input_overrides.items():
            input_data.setdefault(nl, {}).update(d)
    return Espresso(pseudopotentials=pseudos, input_data=input_data, kpts=kpts)


def minimal_cubic_repeat_for_fraction(xA, basis_atoms, maxrep=6, tol=1e-12):
    """
    Find smallest (n,n,n) repeat such that xA * (basis_atoms*n^3) is an integer.
    basis_atoms = number of atoms in the *cubic conventional* cell for that lattice:
      - fcc cubic conventional: 4
      - bcc cubic conventional: 2
    """
    for n in range(1, maxrep + 1):
        N = basis_atoms * (n**3)
        nA = xA * N
        if abs(nA - round(nA)) < tol:
            return (n, n, n), int(round(nA)), N
    raise ValueError(
        f"Cannot represent xA={xA} exactly with (n,n,n) up to n={maxrep}. "
        "Increase maxrep or choose a different cell shape."
    )


def build_ordered_alloy_cubic(A, B, xA, crystal, a_guess, maxrep=6):
    """
    Deterministically build an ordered alloy on a cubic parent lattice (fcc or bcc),
    using the smallest cubic supercell that exactly matches the requested fraction.

    Returns:
      atoms0: ASE Atoms object (ordered alloy)
      sc: repeat tuple (n,n,n)
      counts: (nA, nB, N)
    """
    crystal = crystal.lower()
    if crystal == "fcc":
        basis_atoms = 4  # fcc conventional cubic cell
    elif crystal == "bcc":
        basis_atoms = 2  # bcc conventional cubic cell
    else:
        raise ValueError("This minimal-cell builder supports only 'fcc' and 'bcc'. Use a separate hcp workflow.")

    # choose minimal (n,n,n) so that nA is integer
    sc, nA, N = minimal_cubic_repeat_for_fraction(xA, basis_atoms=basis_atoms, maxrep=maxrep)

    # build cubic conventional cell and repeat
    atoms = bulk(A, crystalstructure=crystal, a=a_guess, cubic=True).repeat(sc)

    # deterministic decoration: first nA sites are A, rest are B
    symbols = [B] * N
    for i in range(nA):
        symbols[i] = A
    atoms.set_chemical_symbols(symbols)

    atoms.set_pbc(True)
    return atoms, sc, (nA, N - nA, N)


def optimize_a_by_energy_minimization(atoms0, calc, a_guess, da=0.10, scan_step=0.01):
    """
    Minimize energy vs lattice constant a for a cubic cell by isotropically scaling the cell.

    atoms0: ordered alloy structure at initial a_guess (or close).
    a_guess: the cubic conventional lattice constant you want to optimize around (Å).
    da: scan half-width in Å.
    """
    # a_len0: current cubic cell length along x for atoms0's *supercell*
    # For cubic cell, scaling factor s relates to a via: a = a_guess * s
    # We assume atoms0 was built using a_guess; then s=1 corresponds to a_guess.
    cell0 = atoms0.get_cell().array
    # take x-vector length as representative for cubic
    Lx0 = np.linalg.norm(cell0[0])

    # supercell repeat factor along x relative to conventional cell:
    # If atoms0 is repeat(sc) of a cubic conventional cell, then Lx0 = sc[0] * a_guess
    # So s = (a_target * sc[0]) / Lx0
    # But since Lx0 ≈ sc[0]*a_guess, this simplifies to a_target/a_guess if built consistently.
    # We'll compute s from a_target robustly:
    scx = int(round(Lx0 / a_guess))

    def energy_at_a(a_target):
        s = (a_target * scx) / Lx0
        atoms = atoms0.copy()
        atoms.set_cell(atoms.get_cell() * s, scale_atoms=True)
        atoms.calc = calc
        return atoms.get_potential_energy()

    # coarse scan
    avals = np.arange(a_guess - da, a_guess + da + 1e-12, scan_step)
    Evals = np.array([energy_at_a(a) for a in avals])

    # quadratic fit near minimum to seed bounded minimization
    inds = np.argsort(Evals)[:7] if len(Evals) >= 7 else np.argsort(Evals)
    p = np.polyfit(avals[inds], Evals[inds], 2)
    a_est = -p[1] / (2.0 * p[0])

    out = opt.minimize_scalar(
        energy_at_a,
        method="bounded",
        bounds=(a_est - 0.02, a_est + 0.02),
        options={"xatol": 1e-4}
    )

    return {
        "a0_A": float(out.x),
        "Emin_eV": float(out.fun),
        "success": bool(out.success),
        "message": out.message,
    }


def get_alloy_lattice_constant_min_energy(
    A, B, xA, surface_type,                 # e.g., A="Pt", B="Au", xA=0.75, surface_type="fcc111"
    calc,                                   # QE calculator
    a_guess,                                # initial guess for fcc/bcc lattice parameter (Å)
    crystal=None,                           # optional override: "fcc" or "bcc"
    maxrep=6,                               # max cubic repeat to match fraction exactly
    da=0.10, scan_step=0.01
):
    crystal = infer_crystal_from_surface(surface_type, crystal=crystal)
    if crystal not in ("fcc", "bcc"):
        raise ValueError("This function currently supports fcc/bcc only. For hcp use a separate (a,c) optimizer.")

    atoms0, sc, counts = build_ordered_alloy_cubic(A, B, xA, crystal, a_guess, maxrep=maxrep)
    optres = optimize_a_by_energy_minimization(atoms0, calc, a_guess, da=da, scan_step=scan_step)

    return {
        "crystal": crystal,
        "surface_type": surface_type,
        "supercell_repeat": sc,
        "counts": {"nA": counts[0], "nB": counts[1], "N": counts[2]},
        **optres
    }