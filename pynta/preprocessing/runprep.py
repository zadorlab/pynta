# preprocessing_input.py
import os
import numpy as np
from ase.io import read
from preprocessing import Prep

###
def mypause(flag):
    print("\033[31;49;1m {0:20s}\033[0m".format(flag))


def run():

    command = '/home/shikim/miniconda3/envs/pynta_env/bin/mpirun pw.x -input < PREFIX.pwi > PREFIX.pwo'  # QE run command

    prep = Prep(
        metal="Pd",
        surface_type="fcc111",
        a0=3.89,
        software="Espresso",
        fmax=0.05,
        vacuum=10,
        software_kwargs={  # this keyword is to run convergence test
            "kpts": (3, 3, 1),
            "ecutwfc": 40,
            "ecutrho": 160,
            "occupations": "smearing",
            "pseudo_dir": '/home/shikim/espresso/pseudo',
            "smearing": "marzari-vanderbilt",
            "degauss": 0.01,
            "input_dft": "BEEF-VDW",
            "nosym": True,
            "conv_thr": 1e-6,
            "pseudopotentials": {
                "Pd": "pd_pbesol_v1.4.uspp.F.UPF",
                "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
            },
            "command": command,
        },
        lattice_opt_software_kwargs={  # this keyword is specifically for lattice constant optimization
            "kpts": (25, 25, 25),  # something large
            "ecutwfc": 70,  # you can change this
            "degauss": 0.02,
            "input_dft": "BEEF-VDW",
            "mixing_mode": "plain",
            "pseudo_dir": '/home/shikim/espresso/pseudo',
            "pseudopotentials": {
                "Pd": "pd_pbesol_v1.4.uspp.F.UPF",
                "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
            },
            "command": command,
        },
    )

    ecut_values = range(30, 71, 10)
    kmesh_values = [(2, 2, 1), (3, 3, 1), (4, 4, 1)]

    # Conventional slab: optimize lattice constant, generate slab, then convergence test
    print("Running conventional slab convergence tests...")
    prep.calculate_lattice_parameters()
    prep.run_lattice_constant_scan(da=0.1)
    prep.generate_slab()
    prep.run_convergence(ecut_values=ecut_values, kmesh_values=kmesh_values)

    # Custom slab convergence tests if the slab file is available
    custom_slab_path = "my_custom_slab.xyz"
    if os.path.exists(custom_slab_path):
        print("Running custom slab convergence tests...")
        prep.slab = read(custom_slab_path)
        prep.run_convergence(ecut_values=ecut_values, kmesh_values=kmesh_values)
    else:
        print(f"Custom slab file not found: {custom_slab_path}. Skipping custom slab convergence.")

    print("Preprocessing finished")


if __name__ == '__main__':
    mypause("Before Run")
    run()
    mypause("After  Run")

