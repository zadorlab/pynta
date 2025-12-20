
import numpy as np
import yaml
import matplotlib.pyplot as plt
import os
import json
# ase general
from ase import Atoms
from ase.io import read, write, Trajectory
from ase.build import add_adsorbate, bulk
from ase.constraints import FixAtoms
from ase.optimize import BFGSLineSearch
from ase.calculators.espresso import Espresso
from ase.constraints import FixAtoms
from ase import build
from ase.data import chemical_symbols, reference_states
from scipy import optimize as opt

# call pynta module
from pynta.utils import name_to_ase_software
from pynta.tasks import *
from pynta.calculator import get_lattice_parameters
# fireworks
from fireworks import LaunchPad, Workflow
from fireworks.queue.queue_launcher import rapidfire as rapidfirequeue
from fireworks.utilities.fw_serializers import load_object_from_file
from fireworks.core.fworker import FWorker

class Prep:
    def __init__(self, metal='Pt', surface_type='fcc111', a0=3.96, software='Espresso', fmax=0.05, adsorbate=None, position=None, vacuum=10, slab=None,
                 launchpad_path=None, fworker_path=None, reset_launchpad=False, queue_adapter_path=None, queue=False, njobs_queue=0, num_jobs=25,
                 software_kwargs=None):
        if software_kwargs is None:
            software_kwargs = {
                'kpts': (3, 3, 1),
                'tprnfor': True,
                'occupations': 'smearing',
                'smearing': 'marzari-vanderbilt',
                'degauss': 0.01,
                'ecutwfc': 40,
                'nosym': True,
                'conv_thr': 1e-6,
                'mixing_mode': 'local-TF',
                "pseudopotentials": {
                    "Cu": 'Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                    "H": 'H.pbe-kjpaw_psl.1.0.0.UPF',
                    "O": 'O.pbe-n-kjpaw_psl.1.0.0.UPF',
                    "C": 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
                    "N": 'N.pbe-n-kjpaw_psl.1.0.0.UPF',
                },
            }
        
        self.metal = metal
        self.surface_type = surface_type
        self.a0 = a0
        self.software = software
        self.fmax = fmax
        self.adsorbate = adsorbate
        self.position = position
        self.vacuum = vacuum
        self.slab = slab
        self.launchpad_path = launchpad_path  # we are not using fireworks at the moment
        self.fworker_path = fworker_path  # we are not using fireworks at the moment
        self.reset_launchpad = reset_launchpad  # we are not using fireworks at the moment
        self.queue_adapter_path = queue_adapter_path  # we are not using fireworks at the moment
        self.queue = queue  # we are not using fireworks at the moment
        self.njobs_queue = njobs_queue  # we are not using fireworks at the moment
        self.num_jobs = num_jobs  # we are not using fireworks at the moment
        self.software_kwargs = software_kwargs

    def freeze_bottom_half(self):
        """Fix atoms in the bottom half of the slab."""
        z = self.slab.positions[:, 2]
        z_mid = 0.5 * (z.max() + z.min())
        mask = z < z_mid
        self.slab.set_constraint(FixAtoms(mask=mask))
        return self.slab

    def optimize_slab(self, a, c=None, slab_path=None, out_path="slab.xyz"):
        """
        Optimize a slab using the calculator specified by `self.software`.

        Parameters
        ----------
        a : float
            Lattice constant
        c : float or None
            c lattice parameter (if applicable)
        slab_path : str or None
            Path to existing slab
        out_path : str
            Output file for optimized slab
        """

        # ------------------------
        # 1. Set up calculator
        # ------------------------
        soft_cls = name_to_ase_software(self.software)
        calc = soft_cls(**self.software_kwargs)

        # ------------------------
        # 2. Load or build slab
        # ------------------------
        if slab_path is not None:
            self.slab = read(slab_path)
        else:
            slab_builder = getattr(build, self.surface_type)

            if c is not None:
                self.slab = slab_builder(
                    symbol=self.metal,
                    size=(3, 3, 4),
                    a=a,
                    c=c,
                    vacuum=self.vacuum
                )
            else:
                self.slab = slab_builder(
                    symbol=self.metal,
                    size=(3, 3, 4),
                    a=a,
                    vacuum=self.vacuum
                )

        self.slab.pbc = (True, True, True)
        self.slab.calc = calc

        # ------------------------
        # 3. Freeze bottom half
        # ------------------------
        if self.software != "XTB":
            self.freeze_bottom_half()

            dyn = BFGSLineSearch(
                self.slab,
                trajectory="slab.traj"
            )
            dyn.run(fmax=self.fmax, steps=200)

        # ------------------------
        # 4. Save and return
        # ------------------------
        write(out_path, self.slab)
        return self.slab

    def run_convergence(self, ecut_values, kmesh_values):
        """Run ecut and k-point convergence tests.

        Returns:
            results: list of (ecut, kpts, energy)
        """
        results = []

        for ecut in ecut_values:
            for kpts in kmesh_values:
                input_data = self.software_kwargs.copy()
                input_data['ecutwfc'] = ecut
                input_data['ecutrho'] = 4 * ecut

                calc = Espresso(
                    input_data=input_data,
                    pseudopotentials=self.software_kwargs['pseudopotentials'],
                    kpts=kpts,
                    profile=None  # You can modify this as needed
                )

                self.slab.set_calculator(calc)  # Use self.slab instead of atoms
                E = self.slab.get_potential_energy()
                print(f"ecut={ecut} Ry, kpts={kpts}, Energy={E:.6f} eV")
                results.append((ecut, kpts, E))

        return results