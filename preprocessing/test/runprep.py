# runprep.py
import numpy as np
from ase.io import read
from preprocessing.prep_workflow import Prep

def mypause(flag):
    print("\033[31;49;1m {0:20s}\033[0m".format(flag))

def run():
    command = '/home/shikim/miniconda3/envs/pynta_env/bin/mpirun pw.x -input < PREFIX.pwi > PREFIX.pwo'

    prep = Prep(
        # ---- base settings ----
        metal="Pd",                 # used for slab generation; for an alloy slab you'd likely change this logic
        surface_type="fcc111",
        a0=3.89,                    # initial guess for lattice constant (used as a_guess by alloy optimizer)
        software="Espresso",
        fmax=0.05,
        vacuum=10,

        # ---- turn on alloy lattice optimization step ----
        optimize_alloy_lattice=True, #default is false

        # ---- choose ONE of the two alloy specification modes ----
        # Mode 1: ordered alloy builder (A/B/xA)
        alloy_A="Pt",
        alloy_B="Au",
        alloy_xA=0.75,
        alloy_maxrep=6,
        alloy_scan_step=0.01,

        # Mode 2: read a bulk alloy structure from file (comment out Mode 1 if using this)
        # alloy_path="bulk_alloy.xyz",

        # ---- QE runtime for slab / other steps ----
        software_kwargs={
            "kpts": (3, 3, 1),
            "ecutwfc": 40,
            "ecutrho": 160,
            "occupations": "smearing",
            "pseudo_dir": "/home/shikim/espresso/pseudo",
            "smearing": "marzari-vanderbilt",
            "degauss": 0.01,
            "input_dft": "BEEF-VDW",
            "nosym": True,
            "conv_thr": 1e-6,
            "pseudopotentials": {
                "Pd": "pd_pbesol_v1.4.uspp.F.UPF",
                "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
                # include alloy pseudos too if you later generate alloy slabs with Espresso:
                "Pt": "Pt.pbe-spn-kjpaw_psl.1.0.0.UPF",
                "Au": "Au.pbe-n-kjpaw_psl.1.0.0.UPF",
            },
            "command": command,
        },

        # ---- QE runtime for bulk lattice optimization step ----
        lattice_opt_software_kwargs={
            "kpts": (25, 25, 25),
            "ecutwfc": 70,
            "ecutrho": 280,          # (optional) set explicitly; your default logic does 4*ecut
            "degauss": 0.02,
            "input_dft": "BEEF-VDW",
            "mixing_mode": "plain",
            "pseudo_dir": "/home/shikim/espresso/pseudo",
            "pseudopotentials": {
                "Pt": "Pt.pbe-spn-kjpaw_psl.1.0.0.UPF",
                "Au": "Au.pbe-n-kjpaw_psl.1.0.0.UPF",
                # if you also do pure Pd bulk optimizations elsewhere:
                "Pd": "pd_pbesol_v1.4.uspp.F.UPF",
            },
            "command": command,
        },
    )

    # ---- Run just the steps enabled by your args (includes calculate_alloy_lattice_parameters) ----
    # This will set prep.a to optimized alloy a0 and write provenance.json entries.
    prep.run_selected()

    print("Optimized lattice constant in Prep (prep.a):", prep.a)

    # ---- Continue with your existing flow (only works if those methods exist in your codebase) ----
    prep.conventional_slab_opt_convergence(
        ecut_values=range(30, 71, 10),
        kmesh_values=[(2,2,1), (3,3,1), (4,4,1)],
    )

    prep.slab = read("my_custom_slab.xyz")
    prep.custom_slab_opt_convergence(
        ecut_values=range(30, 71, 10),
        kmesh_values=[(2,2,1), (3,3,1), (4,4,1)],
    )

    print("Preprocessing finished")

if __name__ == "__main__":
    mypause("Before Run")
    run()
    mypause("After  Run")