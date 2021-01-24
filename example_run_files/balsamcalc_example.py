#!/usr/bin/env python3
from ase.build import molecule
from ase.optimize import QuasiNewton
from rmgcat_to_sella.balsamcalc import EspressoBalsam
from pathlib import Path

# Set up a small, simple system
atoms = molecule('CH4')
atoms.rattle()
atoms.center(vacuum=3)

cwd = Path.cwd().as_posix()

# The calculator
atoms.calc = EspressoBalsam(
    # Special BalsamCalculator keywords:
    workflow='qetest',
    job_kwargs={
        'num_nodes': 1,
        'ranks_per_node': 1,
        'threads_per_rank': 1,  # OMP threading
        'threads_per_core': 1,  # hyperthreads
        'user_workdir': cwd
    },
    # Regular ASE Calculator keywords:
    pseudopotentials={
        'H': 'H.pbe-kjpaw_psl.1.0.0.UPF',
        'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF'
    },
    ecutwfc=40,
    conv_thr=1e-11,
    nosym=True,
    tprnfor=True,
    label='test_espresso_balsam_calculator'
)


# Just do a normal geometry minimization
opt = QuasiNewton(atoms)
opt.run(fmax=0.01)
