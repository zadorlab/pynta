import numpy as np
import os
import matplotlib.pyplot as plt
import json

from scipy import optimize as opt
from ase import Atoms, build
from ase.io import read, write
from ase.build import add_adsorbate, bulk
from ase.constraints import FixAtoms
from ase.optimize import BFGSLineSearch
from ase.data import chemical_symbols, reference_states
from ase.calculators.espresso import Espresso

from pynta.utils import name_to_ase_software
from pynta.calculator import get_lattice_parameters

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as opt


def fit_lattice_constant_from_scan(
    outavals,
    Evals,
    n_fit=7,
    bounds_width=0.01,
    title="Lattice constant optimization",
    show_plot=True,
):
    """
    Post-process lattice scan from get_lattice_parameters().

    Parameters
    ----------
    outavals : list or ndarray
        Lattice constants returned by get_lattice_parameters()
    Evals : list or ndarray
        Energies returned by get_lattice_parameters()
    n_fit : int
        Number of lowest-energy points used for quadratic fit
    bounds_width : float
        Half-width of minimization bounds (Å)
    title : str
        Plot title
    show_plot : bool
        Whether to display the plot

    Returns
    -------
    a_interp : float
        Quadratic minimum
    a_opt : float
        Refined numerical minimum
    """

    outavals = np.asarray(outavals)
    Evals = np.asarray(Evals)

    # -----------------------
    # 1. Select lowest energies
    # -----------------------
    inds = np.argsort(Evals)[:n_fit]

    # -----------------------
    # 2. Quadratic fit
    # -----------------------
    p = np.polyfit(outavals[inds], Evals[inds], 2)

    # Analytic minimum
    a_interp = -p[1] / (2.0 * p[0])

    # -----------------------
    # 3. Numerical refinement
    # -----------------------
    def f(x):
        return np.polyval(p, x)

    out = opt.minimize_scalar(
        f,
        method="bounded",
        bounds=(a_interp - bounds_width, a_interp + bounds_width),
        options={"xatol": 1e-4},
    )

    a_opt = out.x

    # -----------------------
    # 4. Plot
    # -----------------------
    if show_plot:
        x_vals = np.linspace(
            a_interp - 2 * bounds_width,
            a_interp + 2 * bounds_width,
            200,
        )
        y_vals = f(x_vals)

        plt.figure(figsize=(10, 6))
        plt.plot(x_vals, y_vals, label="Quadratic fit", linewidth=2)
        plt.scatter(outavals, Evals, color="orange", label="DFT data", zorder=3)

        plt.axvline(a_interp - bounds_width, linestyle="--", color="red", label="Lower bound")
        plt.axvline(a_interp + bounds_width, linestyle="--", color="green", label="Upper bound")

        plt.scatter(a_opt, f(a_opt), color="red", s=80, label="Optimized a", zorder=5)

        plt.xlabel("Lattice constant a (Å)")
        plt.ylabel("Total energy (eV)")
        plt.title(title)
        plt.grid(True)
        plt.legend()
        plt.show()

    return a_interp, a_opt


class Prep:
    def __init__(
        self,
        metal="Pt",
        surface_type="fcc111",
        a0=3.96,
        c=None,
        repeats=(3, 3, 4),
        pbc=(True, True, False),
        software="Espresso",
        fmax=0.05,
        adsorbate=None,
        position=None,
        vacuum=10,
        slab=None,
        software_kwargs=None,
        lattice_opt_software_kwargs=None,
        path=None,
    ):
        if software_kwargs is None:
            software_kwargs = {
                "kpts": (3, 3, 1),
                "tprnfor": True,
                "occupations": "smearing",
                "smearing": "marzari-vanderbilt",
                "degauss": 0.01,
                "ecutwfc": 40,
                "nosym": True,
                "conv_thr": 1e-6,
                "input_dft": "BEEF-VDW",
                "mixing_mode": "local-TF",
                "pseudopotentials": {
                    "Cu": "Cu.pbe-spn-kjpaw_psl.1.0.0.UPF",
                    "H": "H.pbe-n-kjpaw_psl.1.0.0.UPF",
                    "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
                    "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
                    "N": "N.pbe-n-kjpaw_psl.1.0.0.UPF",
                    "Pt": "Pt.pbe-spn-kjpaw_psl.1.0.0.UPF",
                },
            }

        if lattice_opt_software_kwargs is None:
            lattice_opt_software_kwargs = {
                "kpts": (25, 25, 25),
                "ecutwfc": 70,
                "ecutrho": 280,
                "degauss": 0.02,
                "input_dft": "BEEF-VDW",
                "mixing_mode": "plain",
                "pseudopotentials": software_kwargs["pseudopotentials"],
            }

        self.metal = metal
        self.surface_type = surface_type
        self.software = software
        self.fmax = fmax

        self.a0 = a0
        self.a = a0
        self.c = c

        self.repeats = repeats
        self.pbc = pbc
        self.vacuum = vacuum

        self.adsorbate = adsorbate
        self.position = position
        self.slab = slab

        self.software_kwargs = software_kwargs
        self.lattice_opt_software_kwargs = lattice_opt_software_kwargs

        self.path = path or os.getcwd()
        self.slab_path = None

    def calculate_lattice_parameters(self, da=0.1, a0=None, return_scan=True):
        """
        Bulk lattice parameter optimization.
        Returns a_opt and (optionally) scan arrays outavals, Evals.
        This is modification of get_lattice_parameters() in calculator.py. We can modify this further
        """
        soft = name_to_ase_software(self.software)(**self.lattice_opt_software_kwargs)

        if self.surface_type == "hcp0001":
            raise NotImplementedError("hcp0001 (a,c) scan not wired into this plotting workflow yet.")

        options = {"xatol": 1e-4}

        def f(a):
            b = bulk(self.metal, self.surface_type[:3], a=a, cubic=True) #run cube cell
            b.calc = soft
            b.pbc = (True, True, True)
            return b.get_potential_energy()

        if a0 is None:
            a0 = self.a0 if self.a0 is not None else reference_states[chemical_symbols.index(self.metal)]["a"]

        avals = np.arange(a0 - da, a0 + da, 0.01)
        outavals, Evals = [], []

        for a in avals:
            try:
                E = f(a)
                outavals.append(a)
                Evals.append(E)
            except Exception:
                pass

        inds = np.argsort(Evals)[:7]
        p = np.polyfit(np.array(outavals)[inds], np.array(Evals)[inds], 2)
        a_interp = -p[1] / (2.0 * p[0])

        out = opt.minimize_scalar(
            f,
            method="bounded",
            bounds=(a_interp - 0.01, a_interp + 0.01),
            options=options,
        )

        a_opt = out.x

        # ---- SAVE RESULTS (compute node) ----
        log_data = {
            "metal": self.metal,
            "surface_type": self.surface_type,
            "a0": a0,
            "a_opt": float(a_opt),
            "outavals": list(map(float, outavals)),
            "Evals": list(map(float, Evals)),
        }

        with open("lattice_constant_convergence.out", "w") as f:
            f.write("# Lattice constant convergence\n")
            f.write(f"# Metal: {self.metal}\n")
            f.write(f"# Surface type: {self.surface_type}\n")
            f.write(f"# Reference a0 (Å): {a0:.6f}\n")
            f.write(f"# Optimized a (Å): {a_opt:.6f}\n")
            f.write("#\n")
            f.write("# a (Å)        Energy (eV)\n")
            f.write("# ------------------------\n")

            for a, E in zip(outavals, Evals):
                f.write(f"{a:12.6f}  {E:16.8f}\n")
        
        if return_scan:
            return a_opt, outavals, Evals
        return a_opt

    def freeze_bottom_half(self):
        """Fix atoms in the bottom half of the slab."""
        z = self.slab.positions[:, 2]
        z_mid = 0.5 * (z.max() + z.min())
        mask = z < z_mid
        self.slab.set_constraint(FixAtoms(mask=mask))
        return self.slab

    def optimize_slab(self, a, c=None, slab_path=None, out_path="custom_slab.xyz"):
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

    def generate_slab(self, a=None):
        """
        Generate and optimize a slab locally using ASE (no FireWorks).
        This module is revised from pyn.main generate_slab()
        If `a` is provided, use it directly.
        Otherwise, determine lattice constant internally.        
        """

        # ------------------------
        # 1. Lattice constant
        # ------------------------
        if a is not None:
            self.a = a
        elif self.a is None:
            self.a, _, _ = self.calculate_lattice_parameters(return_scan=True)

        # ------------------------
        # 2. Build slab
        # ------------------------
        slab_type = getattr(build, self.surface_type)

        if self.c is not None:
            slab = slab_type(
                symbol=self.metal,
                size=self.repeats,
                a=self.a,
                c=self.c,
                vacuum=self.vacuum,
            )
        else:
            slab = slab_type(
                symbol=self.metal,
                size=self.repeats,
                a=self.a,
                vacuum=self.vacuum,
            )

        slab.pbc = self.pbc
        self.slab = slab

        # ------------------------
        # 3. Optimize slab
        # ------------------------
        if self.software != "XTB":
            self.slab.calc = name_to_ase_software(self.software)(**self.software_kwargs)
            self.freeze_bottom_half()

            dyn = BFGSLineSearch(self.slab, trajectory="slab.traj")
            dyn.run(fmax=self.fmax, steps=200)

        write("slab.xyz", self.slab)
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

    def run_convergence_adsorbate(
        self,
        ecut_values,
        kmesh_values,
        adsorbate=None,
        height=1.0,
        position="ontop",
    ):
        """
        Run ecut and k-point convergence tests using adsorption energy.

        Parameters
        ----------
        adsorbate : None | str | ase.Atoms --> given by user
            Adsorbate species. If None, defaults to atomic H.
        height : float
            Adsorption height in Å
        position : str or tuple
            Adsorption site

        Returns
        -------
        results : list of (ecut, kpts, E_ads)
        """

        results = []

        # -----------------------
        # 0. Default adsorbate
        # -----------------------
        if adsorbate is None:
            adsorbate = Atoms("H")
            is_hydrogen = True
        else:
            is_hydrogen = False
            if isinstance(adsorbate, str):
                adsorbate = Atoms(adsorbate)

        slab_clean = self.slab.copy()

        for ecut in ecut_values:
            for kpts in kmesh_values:

                # -----------------------
                # 1. Calculator
                # -----------------------
                input_data = self.software_kwargs.copy()
                input_data["ecutwfc"] = ecut
                input_data["ecutrho"] = 4 * ecut

                calc = Espresso(
                    input_data=input_data,
                    pseudopotentials=self.software_kwargs["pseudopotentials"],
                    kpts=kpts,
                    profile=None,
                )

                # -----------------------
                # 2. Bare slab
                # -----------------------
                slab = slab_clean.copy()
                slab.calc = calc
                E_slab = slab.get_potential_energy()

                # -----------------------
                # 3. Slab + adsorbate
                # -----------------------
                slab_ads = slab_clean.copy()
                add_adsorbate(
                    slab_ads,
                    adsorbate,
                    height=height,
                    position=position,
                )
                slab_ads.calc = calc
                E_slab_ads = slab_ads.get_potential_energy()

                # -----------------------
                # 4. Reference energy
                # -----------------------
                if is_hydrogen:
                    # H reference = 1/2 H2
                    ref = Atoms(
                        "H2",
                        positions=[[0, 0, 0], [0, 0, 0.74]],
                        cell=[10, 10, 10],
                        pbc=False,
                    )
                    ref.calc = calc
                    E_ref = 0.5 * ref.get_potential_energy()
                else:
                    # Gas-phase molecule reference
                    ref = adsorbate.copy()
                    ref.cell = [10, 10, 10]
                    ref.pbc = False
                    ref.calc = calc
                    E_ref = ref.get_potential_energy()

                # -----------------------
                # 5. Adsorption energy
                # -----------------------
                E_ads = E_slab_ads - E_slab - E_ref

                print(
                    f"ecut={ecut} Ry, kpts={kpts}, "
                    f"E_ads={E_ads:.6f} eV"
                )

                results.append((ecut, kpts, E_ads))

                # -----------------------
                # 6. Save results
                # -----------------------
                if save_results:
                    with open(outfile, "w") as f:
                        f.write("# ecut (Ry)   kpts        E_ads (eV)\n")
                        for ecut, kpts, E_ads in results:
                            f.write(f"{ecut:<12} {str(kpts):<12} {E_ads:.8f}\n")

        return results

