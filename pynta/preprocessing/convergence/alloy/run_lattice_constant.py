# Example: Pt/Au alloy, 75% Pt and 25% Au, fcc111 -> optimize fcc bulk lattice constant a0

from ase.calculators.espresso import Espresso

# 1) Define QE pseudopotentials (filenames must be available in your run directory or pseudo_dir)
pseudos = {
    "Pt": "Pt.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Au": "Au.pbe-n-kjpaw_psl.1.0.0.UPF",
}

# 2) Make the QE calculator
calc = make_qe_calc(
    pseudos=pseudos,
    kpts=(12, 12, 12),
    ecutwfc=60,
    ecutrho=480,
    degauss=0.02,
    # optional overrides:
    # qe_input_overrides={"control": {"pseudo_dir": "/path/to/pseudos"}}
)

# 3) Call the function
res = get_alloy_lattice_constant_min_energy(
    A="Pt",                 # element A
    B="Au",                 # element B
    xA=0.75,                # 75% Pt, 25% Au
    surface_type="fcc111",  # used to infer the parent lattice (fcc)
    calc=calc,
    a_guess=3.95,           # initial guess for fcc lattice constant in Å
    maxrep=6,               # will choose smallest (n,n,n) supercell that matches 0.75 exactly
    da=0.12,                # scan +/- 0.12 Å around a_guess
    scan_step=0.01          # coarse scan step in Å
)

# 4) Read results
print("Crystal:", res["crystal"])
print("Supercell repeat:", res["supercell_repeat"])
print("Atom counts:", res["counts"])     # e.g. N=4, nA=3, nB=1 for 75/25 on fcc conventional cell
print("Optimized a0:", res["a0_A"], "Å")
print("Min energy:", res["Emin_eV"], "eV")
print("Success:", res["success"])
print("Message:", res["message"])