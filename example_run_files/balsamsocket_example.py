#!/usr/bin/env python3
from ase.build import molecule
from ase.optimize import QuasiNewton
from rmgcat_to_sella.balsamcalc import EspressoBalsamSocketIO


# Set up a small, simple system
atoms = molecule('CH4')
atoms.rattle()
atoms.center(vacuum=3)

# The calculator
atoms.calc = EspressoBalsamSocketIO(
    # Special BalsamCalculator keywords:
    workflow='qetest',
    job_kwargs={
        'num_nodes': 1,
        'ranks_per_node': 6
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

atoms.calc.close()
