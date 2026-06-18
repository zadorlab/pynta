
"""
runprep.py

Example driver script for preprocessing workflow.
Safe to re-run: automatically resumes from prep_convergence.json.
"""

from ase.io import read
from preprocessing import Prep


# ------------------------------------------------------------
# Optional: user-defined slab builder
# ------------------------------------------------------------

def my_custom_slab_builder(
    metal,
    a,
    repeats,
    vacuum,
    pbc,
    c=None,
):
    """
    Example user slab builder.
    Must return an ase.Atoms object.
    """
    from ase import build

    slab = build.fcc111(
        metal,
        size=repeats,
        a=a,
        vacuum=vacuum,
    )
    slab.center(axis=2)
    slab.pbc = pbc
    return slab


# ------------------------------------------------------------
# Main workflow
# ------------------------------------------------------------

def run():
    command = '/home/shikim/miniconda3/envs/pynta_env/bin/mpirun pw.x -input < PREFIX.pwi > PREFIX.pwo' #QE run command

    prep = Prep(
        # --------------------
        # System definition
        # --------------------
        metal="Pt",
        surface_type="fcc111",     # used if slab_builder not provided
        repeats=(3, 3, 3),
        vacuum=10.0,
        fmax=0.05,
        frozen_layers=3,

        # --------------------
        # Workflow flags
        # --------------------
        #optimize_lattice=True,     # run lattice optimization
        #generate_slab=True,        # build + optimize slab
        #run_convergence_adsorbate=False,

        # --------------------
        # Alloy parameters (optional)
        # --------------------
        #alloy_path=None,           # path to pre-built alloy structure file
        #alloy_A="Pt",              # majority element
        #alloy_B="Au",              # minority element
        #alloy_xA=0.75,             # fraction of A (e.g. 0.75 = 75% Pt)
        #alloy_maxrep=6,
        #alloy_scan_step=0.01,
        #optimize_alloy_lattice=True,

        # --------------------
        # Optional user inputs
        # --------------------
        # a=3.89,                # uncomment to skip lattice optimization
        # slab_builder=my_custom_slab_builder,  # comment out to use ASE builder

        # --------------------
        # Convergence parameters
        # --------------------
        ecut_values=range(30, 71, 10),
        kmesh_values=[(2, 2, 1), (3, 3, 1), (4, 4, 1)],

        # --------------------
        # QE / ASE calculator inputs
        # --------------------
        software="Espresso",
        software_kwargs={
            "kpts": (3, 3, 1),
            "ecutwfc": 40,
            "ecutrho": 160,
            "occupations": "smearing",
            "smearing": "marzari-vanderbilt",
            "degauss": 0.01,
            "input_dft": "BEEF-VDW",
            "nosym": True,
            "conv_thr": 1e-6,
            "pseudo_dir" : '/home/shikim/espresso/pseudo',
            "pseudopotentials": {
                "Pt": "Pt.pbe-spn-kjpaw_psl.1.0.0.UPF",
                "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
            },
            "command": command
        },

        lattice_opt_software_kwargs={
            "kpts": (25, 25, 25),
            "ecutwfc": 70,
            "ecutrho": 280,
            "degauss": 0.02,
            "input_dft": "BEEF-VDW",
            "mixing_mode": "plain",
            "pseudo_dir" : '/home/shikim/espresso/pseudo',
            "pseudopotentials": {
                "Pt": "Pt.pbe-spn-kjpaw_psl.1.0.0.UPF",
            },
            'command':command
        },
    )

    # --------------------------------------------------------
    # Run workflow (restart-safe)
    # --------------------------------------------------------
    prep.run_selected()

    print("Preprocessing finished successfully")


# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------

if __name__ == "__main__":
    run()