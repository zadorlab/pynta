import os
import json
import time
import inspect
import logging
import numpy as np

from scipy import optimize as opt
from ase import Atoms, build
from ase.io import read, write
from ase.build import add_adsorbate, bulk
from ase.constraints import FixAtoms
from ase.optimize import BFGSLineSearch
from ase.data import chemical_symbols, reference_states
from ase.calculators.espresso import Espresso
from ase.constraints import FixAtoms

from pynta.utils import name_to_ase_software


# ============================================================
# Logging
# ============================================================

logger = logging.getLogger("preprocessing")
logger.setLevel(logging.INFO)

fmt = logging.Formatter(
    "%(asctime)s | %(levelname)s | %(message)s"
)

fh = logging.FileHandler("preprocessing.log")
fh.setFormatter(fmt)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setFormatter(fmt)
logger.addHandler(ch)


# ============================================================
# Workflow decorators
# ============================================================

def requires(*argnames):
    def decorator(func):
        func._requires = set(argnames)
        return func
    return decorator


def step(order):
    def decorator(func):
        func._order = order
        return func
    return decorator

def fit_lattice_constant_from_scan(
    outavals,
    Evals,
    n_fit=7,
    bounds_width=0.01,
    title="Lattice constant optimization",
    show_plot=False,
    save_plot=True,
    plot_filename="lattice_constant_fit.png",
):
    """
    Post-process lattice scan from get_lattice_parameters().

    Safe for headless environments (HPC / batch jobs).

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
        Display plot interactively (default: False)
    save_plot : bool
        Save plot to disk (default: True)
    plot_filename : str
        Output plot filename

    Returns
    -------
    a_interp : float
        Quadratic minimum
    a_opt : float
        Refined numerical minimum
    """

    import numpy as np
    from scipy import optimize as opt

    # ---- Optional, safe matplotlib import ----
    plotting_available = False
    if show_plot or save_plot:
        try:
            import matplotlib
            if not show_plot:
                matplotlib.use("Agg")  # headless backend
            import matplotlib.pyplot as plt
            plotting_available = True
        except Exception:
            plotting_available = False

    outavals = np.asarray(outavals, dtype=float)
    Evals = np.asarray(Evals, dtype=float)

    # -----------------------
    # 1. Select lowest energies
    # -----------------------
    if len(outavals) < n_fit:
        raise ValueError(
            f"Not enough points for fit: {len(outavals)} < n_fit={n_fit}"
        )

    inds = np.argsort(Evals)[:n_fit]

    # -----------------------
    # 2. Quadratic fit
    # -----------------------
    p = np.polyfit(outavals[inds], Evals[inds], 2)

    if p[0] <= 0:
        raise RuntimeError(
            "Quadratic fit curvature is non-positive; "
            "lattice scan may be insufficient."
        )

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

    if not out.success:
        raise RuntimeError("Lattice constant minimization failed")

    a_opt = float(out.x)

    # -----------------------
    # 4. Plot (optional)
    # -----------------------
    if plotting_available:
        x_vals = np.linspace(
            a_interp - 2 * bounds_width,
            a_interp + 2 * bounds_width,
            200,
        )
        y_vals = f(x_vals)

        plt.figure(figsize=(10, 6))
        plt.plot(x_vals, y_vals, label="Quadratic fit", linewidth=2)
        plt.scatter(outavals, Evals, color="orange", label="DFT data", zorder=3)

        plt.axvline(
            a_interp - bounds_width,
            linestyle="--",
            color="red",
            label="Lower bound",
        )
        plt.axvline(
            a_interp + bounds_width,
            linestyle="--",
            color="green",
            label="Upper bound",
        )

        plt.scatter(
            a_opt,
            f(a_opt),
            color="red",
            s=80,
            label="Optimized a",
            zorder=5,
        )

        plt.xlabel("Lattice constant a (Å)")
        plt.ylabel("Total energy (eV)")
        plt.title(title)
        plt.grid(True)
        plt.legend()

        if save_plot:
            plt.savefig(plot_filename, dpi=300, bbox_inches="tight")

        if show_plot:
            plt.show()

        plt.close()

    return a_interp, a_opt


# ============================================================
# Prep class
# ============================================================

class Prep:

    # --------------------------------------------------------
    # Initialization + provenance loading
    # --------------------------------------------------------

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
        frozen_layers = None,
        **kwargs,

    ):

        # ---- Detect user-provided arguments ----
        sig = inspect.signature(self.__init__)
        bound = sig.bind_partial(
            metal=metal,
            surface_type=surface_type,
            a0=a0,
            c=c,
            repeats=repeats,
            pbc=pbc,
            software=software,
            fmax=fmax,
            adsorbate=adsorbate,
            position=position,
            vacuum=vacuum,
            slab=slab,
            software_kwargs=software_kwargs,
            lattice_opt_software_kwargs=lattice_opt_software_kwargs,
            path=path,
            frozen_layers = frozen_layers
        )

        self.user_args = {
            k: v
            for k, v in bound.arguments.items()
            if sig.parameters[k].default != v
        }
        self.user_args.update(kwargs)

        logger.info("User arguments:")
        for k, v in self.user_args.items():
            logger.info(f"  {k} = {v}")

        # ---- Assign attributes ----
        self.metal = metal
        self.surface_type = surface_type
        self.a0 = a0
        self.a = a0
        self.c = c
        self.repeats = repeats
        self.pbc = pbc
        self.software = software
        self.fmax = fmax
        self.adsorbate = adsorbate
        self.position = position
        self.vacuum = vacuum
        self.slab = slab
        self.layers = self.repeats[2]
        if self.slab is None:
            self.nslab = int(np.prod(np.array(self.repeats)))
        else:
            self.nslab = len(read(self.slab))
        self.frozen_layers = frozen_layers
        self.freeze_ind = int((self.nslab/self.layers)*self.frozen_layers)


        # ---- software_kwargs handling ----
        if software_kwargs is None:
            self.software_kwargs = {
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
            logger.info("Using default software_kwargs")
        else:
            self.software_kwargs = software_kwargs
            logger.info("Using user-provided software_kwargs")

        # ---- lattice optimization kwargs ----
        if lattice_opt_software_kwargs is None:
            self.lattice_opt_software_kwargs = {
                "kpts": (25, 25, 25),
                "ecutwfc": 70,
                "ecutrho": 280,
                "degauss": 0.02,
                "input_dft": "BEEF-VDW",
                "mixing_mode": "plain",
                "pseudopotentials": self.software_kwargs["pseudopotentials"],
            }
            logger.info("Using default lattice_opt_software_kwargs")
        else:
            self.lattice_opt_software_kwargs = lattice_opt_software_kwargs
            logger.info("Using user-provided lattice_opt_software_kwargs")

        self.path = path or os.getcwd()

        # ---- Load or initialize provenance ----
        self.provenance_file = "provenance.json"
        self.provenance = self._load_provenance()

        self.completed_steps = {
            step["step"]
            for step in self.provenance.get("steps", [])
            if step.get("status") == "completed"
        }

        self.restore_state_from_provenance()

    # --------------------------------------------------------
    # Provenance helpers
    # --------------------------------------------------------

    def _load_provenance(self):
        if os.path.exists(self.provenance_file):
            logger.info("Loading existing provenance.json")
            with open(self.provenance_file, "r") as f:
                return json.load(f)

        return {
            "created": time.strftime("%Y-%m-%d %H:%M:%S"),
            "user_arguments": self.user_args,
            "steps": [],
            "outputs": [],
        }

    def _save_provenance(self):
        with open(self.provenance_file, "w") as f:
            json.dump(self.provenance, f, indent=2)

    def _record_step(self, name, inputs=None, outputs=None, status="completed"):
        entry = {
            "step": name,
            "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "status": status,
            "inputs": inputs or {},
            "outputs": outputs or {},
        }
        self.provenance["steps"].append(entry)
        self._save_provenance()

    # --------------------------------------------------------
    # State restoration
    # --------------------------------------------------------

    def restore_state_from_provenance(self):
        logger.info("Restoring state from provenance")

        for step in self.provenance.get("steps", []):
            if step["step"] == "calculate_lattice_parameters":
                self.a = step["outputs"].get("lattice_constant", self.a)

            if step["step"] == "generate_slab":
                slab_file = step["outputs"].get("slab_file")
                if slab_file and os.path.exists(slab_file):
                    self.slab = read(slab_file)

    # --------------------------------------------------------
    # Workflow discovery & execution
    # --------------------------------------------------------

    def discover_runnable_methods(self):
        methods = []

        for _, method in inspect.getmembers(self, predicate=inspect.ismethod):
            if not hasattr(method, "_requires"):
                continue
            if method._requires.issubset(self.user_args.keys()):
                methods.append(method)
        #print methods
        print("Runnable methods", methods)

        return methods

    def run_selected(self, force=False):
        methods = self.discover_runnable_methods()
        methods.sort(key=lambda m: getattr(m, "_order", 999))

        logger.info("Workflow steps:")
        for m in methods:
            logger.info(f"  {m.__name__}")

        for method in methods:
            name = method.__name__

            if name in self.completed_steps and not force:
                logger.info(f"⏭ Skipping completed step: {name}")
                continue

            logger.info(f"▶ Running step: {name}")
            try:
                self._run_method(method)
                self.completed_steps.add(name)
            except Exception as e:
                logger.error(f"✖ Step failed: {name}")
                logger.exception(e)
                self._record_step(
                    name,
                    status="failed",
                    outputs={"error": str(e)}
                )
                raise

        logger.info("Preprocessing workflow finished successfully")

    def _run_method(self, method):
        sig = inspect.signature(method)
        kwargs = {
            k: self.user_args[k]
            for k in sig.parameters
            if k in self.user_args
        }
        return method(**kwargs)

    # ========================================================
    # Axiliary functions
    # ========================================================

    def apply_constraints(self, sp, constraints):
        """
        Apply user-defined constraints to a slab.

        Parameters
        ----------
        sp : ase.Atoms
            Slab structure
        constraints : list
            List of constraints, which may be:
            - dict (passed to construct_constraint)
            - str (human-readable constraint keywords)
        """

        out_constraints = []

        for c in constraints:

            # --------------------------
            # Dict-based constraint
            # --------------------------
            if isinstance(c, dict):
                constraint = construct_constraint(c)
                out_constraints.append(constraint)

            # --------------------------
            # Freeze half slab (by z) 
            # --------------------------
            elif c == "freeze half slab":
                z_mid = sp.cell[2, 2] / 2.0
                out_constraints.append(
                    FixAtoms(
                        indices=[
                            atom.index
                            for atom in sp
                            if atom.position[2] < z_mid
                        ]
                    )
                )

            # --------------------------
            # Freeze all atoms of an element
            # Example: "freeze all Cu"
            # --------------------------
            elif c.startswith("freeze all"):
                sym = c.split()[2]
                out_constraints.append(
                    FixAtoms(
                        indices=[
                            atom.index
                            for atom in sp
                            if atom.symbol == sym
                        ]
                    )
                )

            # --------------------------
            # Freeze bottom N layers (N from self.frozen_layers)
            # Keyword trigger only
            # --------------------------
            elif c == "freeze bottom layers":

                if self.frozen_layers is None:
                    raise ValueError(
                        "'freeze bottom layers' requested but frozen_layers is not set"
                    )

                nslab = len(sp)
                layers = getattr(self, "layers", None) or self.repeats[2]

                if self.frozen_layers > layers:
                    raise ValueError(
                        f"frozen_layers ({self.frozen_layers}) "
                        f"> total layers ({layers})"
                    )

                n_per_layer = nslab // layers
                n_freeze = n_per_layer * self.frozen_layers

                # Robust: sort atoms by z (bottom → top)
                z_sorted_indices = sorted(
                    range(nslab),
                    key=lambda i: sp.positions[i][2]
                )

                out_constraints.append(
                    FixAtoms(indices=z_sorted_indices[:n_freeze])
                )

            else:
                raise ValueError(f"Unknown constraint: {c}")

        if out_constraints:
            sp.set_constraint(out_constraints)

    # ========================================================
    # Actual preprocessing steps
    # ========================================================

    @step(1)
    @requires("optimize_lattice")
    def calculate_lattice_parameters(self, da=0.1):
        logger.info("Optimizing bulk lattice constant")

        calc = name_to_ase_software(self.software)(
            **self.lattice_opt_software_kwargs
        )

        def energy(a):
            atoms = bulk(self.metal, "fcc", a=a, cubic=True)
            atoms.calc = calc
            return atoms.get_potential_energy()

        avals = np.arange(self.a0 - da, self.a0 + da, 0.01)
        energies = [energy(a) for a in avals]

        p = np.polyfit(avals, energies, 2)
        self.a = -p[1] / (2 * p[0])

        logger.info(f"Optimized lattice constant: a = {self.a:.6f} Å")

        self._record_step(
            "calculate_lattice_parameters",
            inputs={"a0": self.a0, "da": da},
            outputs={"lattice_constant": self.a},
        )

        return self.a

    @step(2)
    @requires("generate_slab")
    def generate_slab(self, a=None):
        logger.info("Generating slab")

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
            self.freeze_bottom_n_layers()

            dyn = BFGSLineSearch(self.slab, trajectory="slab.traj")
            dyn.run(fmax=self.fmax, steps=200)

        write("slab.xyz", self.slab)
        return self.slab


    @step(3)
    @requires("ecut_values", "kmesh_values")
    def run_convergence(self, ecut_values, kmesh_values):
        logger.info("Running convergence tests")

        results = []

        for ecut in ecut_values:
            for kpts in kmesh_values:
                input_data = self.software_kwargs.copy()
                input_data["ecutwfc"] = ecut
                input_data["ecutrho"] = 4 * ecut

                calc = Espresso(
                    input_data=input_data,
                    pseudopotentials=self.software_kwargs["pseudopotentials"],
                    kpts=kpts,
                )

                self.slab.calc = calc
                E = self.slab.get_potential_energy()

                logger.info(
                    f"ecut={ecut} Ry | kpts={kpts} | E={E:.6f} eV"
                )

                results.append((ecut, kpts, E))

        self._record_step(
            "run_convergence",
            inputs={
                "ecut_values": list(ecut_values),
                "kmesh_values": list(kmesh_values),
            },
            outputs={"results": results},
        )

        return results

    @step(4)
    @requires("run_adsorbate_convergence", "ecut_values", "kmesh_values")
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
        """

        logger.info("Starting adsorbate convergence test")
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
                # Bare slab
                # -----------------------
                slab = slab_clean.copy()
                slab.calc = calc
                E_slab = slab.get_potential_energy()

                # -----------------------
                # Slab + adsorbate
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
                # Reference
                # -----------------------
                if is_hydrogen:
                    ref = Atoms(
                        "H2",
                        positions=[[0, 0, 0], [0, 0, 0.74]],
                        cell=[10, 10, 10],
                        pbc=False,
                    )
                    ref.calc = calc
                    E_ref = 0.5 * ref.get_potential_energy()
                else:
                    ref = adsorbate.copy()
                    ref.cell = [10, 10, 10]
                    ref.pbc = False
                    ref.calc = calc
                    E_ref = ref.get_potential_energy()

                E_ads = E_slab_ads - E_slab - E_ref

                logger.info(
                    f"ecut={ecut} Ry | kpts={kpts} | E_ads={E_ads:.6f} eV"
                )

                results.append((ecut, kpts, E_ads))

        self._record_step(
            "run_convergence_adsorbate",
            inputs={
                "ecut_values": list(ecut_values),
                "kmesh_values": list(kmesh_values),
                "adsorbate": str(adsorbate),
                "height": height,
                "position": position,
            },
            outputs={"results": results},
        )

        return results          