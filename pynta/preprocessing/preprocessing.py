# preprocessing.py
import os
import json
import time
import inspect
import logging
import numpy as np
import copy

from scipy import optimize as opt

from ase import Atoms, build
from ase.io import read, write, Trajectory
from ase.build import add_adsorbate, bulk
from ase.optimize import BFGSLineSearch
from ase.calculators.espresso import Espresso

from pynta.utils import name_to_ase_software
# NOTE: site utilities / vacancy logic live in site_analysis.py (no duplication)
import site_analysis as sa


# ============================================================
# Logging
# ============================================================

logger = logging.getLogger("preprocessing")
logger.setLevel(logging.INFO)

fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

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


# ============================================================
# Lattice constant scan post-processing
# ============================================================

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
    plotting_available = False
    if show_plot or save_plot:
        try:
            import matplotlib
            if not show_plot:
                matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            plotting_available = True
        except Exception:
            plotting_available = False

    outavals = np.asarray(outavals, dtype=float)
    Evals = np.asarray(Evals, dtype=float)

    if len(outavals) < n_fit:
        raise ValueError(f"Not enough points for fit: {len(outavals)} < n_fit={n_fit}")

    inds = np.argsort(Evals)[:n_fit]
    p = np.polyfit(outavals[inds], Evals[inds], 2)

    if p[0] <= 0:
        raise RuntimeError(
            "Quadratic fit curvature is non-positive; lattice scan may be insufficient."
        )

    a_interp = -p[1] / (2.0 * p[0])

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

    if plotting_available:
        x_vals = np.linspace(a_interp - 2 * bounds_width, a_interp + 2 * bounds_width, 200)
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

        if save_plot:
            plt.savefig(plot_filename, dpi=300, bbox_inches="tight")
        if show_plot:
            plt.show()
        plt.close()

    return a_interp, a_opt


# ============================================================
# Alloy lattice constant helpers
# ============================================================

def infer_crystal_from_surface(surface_type, crystal=None):
    if crystal is not None:
        return crystal.lower()

    st = surface_type.lower()
    if st.startswith("fcc"):
        return "fcc"
    if st.startswith("bcc"):
        return "bcc"
    if st.startswith("hcp"):
        return "hcp"
    raise ValueError(
        f"Cannot infer crystal from surface_type='{surface_type}'. Provide crystal explicitly."
    )


def minimal_cubic_repeat_for_fraction(xA, basis_atoms, maxrep=6, tol=1e-12):
    for n in range(1, maxrep + 1):
        N = basis_atoms * (n**3)
        nA = xA * N
        if abs(nA - round(nA)) < tol:
            return (n, n, n), int(round(nA)), N
    raise ValueError(
        f"Cannot represent xA={xA} exactly with (n,n,n) up to n={maxrep}. "
        "Increase maxrep or choose a different cell."
    )


def build_ordered_alloy_cubic(A, B, xA, crystal, a_guess, maxrep=6):
    crystal = crystal.lower()
    if crystal == "fcc":
        basis_atoms = 4
    elif crystal == "bcc":
        basis_atoms = 2
    else:
        raise ValueError("Only 'fcc' and 'bcc' are supported by this ordered cubic alloy builder.")

    sc, nA, N = minimal_cubic_repeat_for_fraction(xA, basis_atoms=basis_atoms, maxrep=maxrep)
    atoms = bulk(A, crystalstructure=crystal, a=a_guess, cubic=True).repeat(sc)

    symbols = [B] * N
    for i in range(nA):
        symbols[i] = A
    atoms.set_chemical_symbols(symbols)
    atoms.set_pbc(True)

    return atoms, sc, (nA, N - nA, N)


def optimize_a_by_energy_minimization(atoms0, calc, a_guess, da=0.10, scan_step=0.01):
    """
    Minimize energy vs lattice constant a for a cubic cell by isotropically scaling the cell.
    """
    cell0 = atoms0.get_cell().array
    Lx0 = np.linalg.norm(cell0[0])

    if a_guess <= 0:
        raise ValueError("a_guess must be > 0")

    scx = int(round(Lx0 / a_guess))
    scx = max(scx, 1)

    def energy_at_scale(s):
        atoms = atoms0.copy()
        atoms.set_cell(atoms.get_cell() * s, scale_atoms=True)
        atoms.calc = calc
        return atoms.get_potential_energy()

    s_vals = np.arange(
        (a_guess - da) / a_guess,
        (a_guess + da) / a_guess + 1e-12,
        scan_step / a_guess,
    )
    Evals = np.array([energy_at_scale(s) for s in s_vals])

    inds = np.argsort(Evals)[:7] if len(Evals) >= 7 else np.argsort(Evals)
    p = np.polyfit(s_vals[inds], Evals[inds], 2)
    if p[0] <= 0:
        raise RuntimeError("Quadratic fit curvature is non-positive; increase scan window/points.")
    s_est = -p[1] / (2.0 * p[0])

    out = opt.minimize_scalar(
        energy_at_scale,
        method="bounded",
        bounds=(s_est - 0.02, s_est + 0.02),
        options={"xatol": 1e-4},
    )

    s_opt = float(out.x)
    Lx_opt = Lx0 * s_opt
    a0 = Lx_opt / scx

    return {
        "a0_A": float(a0),
        "scale_opt": float(s_opt),
        "Emin_eV": float(out.fun),
        "success": bool(out.success),
        "message": str(out.message),
        "scan_scales": s_vals.tolist(),
        "scan_Evals_eV": Evals.tolist(),
        "scx_est": int(scx),
    }


def get_alloy_lattice_constant_min_energy(
    A, B, xA, surface_type, calc, a_guess,
    crystal=None, maxrep=6, da=0.10, scan_step=0.01
):
    crystal = infer_crystal_from_surface(surface_type, crystal=crystal)
    if crystal not in ("fcc", "bcc"):
        raise ValueError("Only fcc/bcc are supported by get_alloy_lattice_constant_min_energy().")

    atoms0, sc, counts = build_ordered_alloy_cubic(A, B, xA, crystal, a_guess, maxrep=maxrep)
    optres = optimize_a_by_energy_minimization(atoms0, calc, a_guess, da=da, scan_step=scan_step)

    return {
        "crystal": crystal,
        "surface_type": surface_type,
        "supercell_repeat": sc,
        "counts": {"nA": counts[0], "nB": counts[1], "N": counts[2]},
        **optres,
    }


# ============================================================
# High-level workflows (owned by preprocessing.py)
# ============================================================

def workflow_sites_no_defect(
    slab,
    nslab,
    adsorbate_height=1.0,
    site_bond_cutoff=1.5,
    surface_string="fcc332",
    traj_filename="unique_sites.traj",
    sites_json="sites.json",
    neighbor_json="neighbor_site_list.json",
    sites_graph_json="sites_graph.json",
    tag_symbol="He",
    verbose=True,
):
    """
    Wrapper around site_analysis that keeps high-level workflows in preprocessing.py.
    """
    return sa.workflow_no_defect_unique_sites(
        slab=slab,
        nslab=nslab,
        adsorbate_height=adsorbate_height,
        site_bond_cutoff=site_bond_cutoff,
        surface_string=surface_string,
        traj_filename=traj_filename,
        sites_json=sites_json,
        neighbor_json=neighbor_json,
        sites_graph_json=sites_graph_json,
        tag_symbol=tag_symbol,
        verbose=verbose,
    )


def workflow_sites_defect_vacancy_drop(
    slab,
    nslab,
    xyz_path,
    surface_obj,
    site_bond_cutoff=1.5,
    tag_symbol="Ne",
    dz=0.1,
    stable_steps=3,
    max_drop=10.0,
    margin=0.25,
    min_clearance=1.3,
    save_all_drop_steps=True,
    traj_geom_all="geom_all.traj",
    traj_drop="drop_steps.traj",
    traj_maxcn="maxcn_geoms.xyz",
    json_geom_all_sites="geom_all_sites_lists.json",
    neighbor_json="neighbor_site_list.json",
    verbose=True,
):
    return sa.workflow_defect_vacancy_drop(
        slab=slab,
        nslab=nslab,
        xyz_path=xyz_path,
        surface_obj=surface_obj,
        site_bond_cutoff=site_bond_cutoff,
        tag_symbol=tag_symbol,
        dz=dz,
        stable_steps=stable_steps,
        max_drop=max_drop,
        margin=margin,
        min_clearance=min_clearance,
        save_all_drop_steps=save_all_drop_steps,
        traj_geom_all=traj_geom_all,
        traj_drop=traj_drop,
        traj_maxcn=traj_maxcn,
        json_geom_all_sites=json_geom_all_sites,
        neighbor_json=neighbor_json,
        verbose=verbose,
    )


def workflow_sites_auto(
    xyz_path,
    n_layers=4,
    adsorbate_height=1.0,
    site_bond_cutoff=1.5,
    surface_string_no_defect="fcc332",
    tag_symbol="Ne",
    dz=0.1,
    stable_steps=3,
    max_drop=10.0,
    margin=0.25,
    min_clearance=1.3,
    save_all_drop_steps=True,
    verbose=True,
):
    return sa.workflow_auto(
        xyz_path=xyz_path,
        n_layers=n_layers,
        adsorbate_height=adsorbate_height,
        site_bond_cutoff=site_bond_cutoff,
        surface_string_no_defect=surface_string_no_defect,
        tag_symbol=tag_symbol,
        dz=dz,
        stable_steps=stable_steps,
        max_drop=max_drop,
        margin=margin,
        min_clearance=min_clearance,
        save_all_drop_steps=save_all_drop_steps,
        verbose=verbose,
    )


# ============================================================
# Prep class
# ============================================================

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
        frozen_layers=None,
        alloy_path=None,
        alloy_A=None,
        alloy_B=None,
        alloy_xA=None,
        alloy_maxrep=6,
        alloy_scan_step=0.01,
        **kwargs,
    ):
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
            frozen_layers=frozen_layers,
            alloy_path=alloy_path,
            alloy_A=alloy_A,
            alloy_B=alloy_B,
            alloy_xA=alloy_xA,
            alloy_maxrep=alloy_maxrep,
            alloy_scan_step=alloy_scan_step
        )

        self.user_args = {
            k: v for k, v in bound.arguments.items()
            if sig.parameters[k].default != v
        }
        self.user_args.update(kwargs)

        logger.info("User arguments:")
        for k, v in self.user_args.items():
            logger.info(f"  {k} = {v}")

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
        self.freeze_ind = int((self.nslab / self.layers) * self.frozen_layers) if self.frozen_layers is not None else None

        self.alloy_path = alloy_path
        self.alloy_A = alloy_A
        self.alloy_B = alloy_B
        self.alloy_xA = alloy_xA
        self.alloy_maxrep = alloy_maxrep
        self.alloy_scan_step = alloy_scan_step

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

        self.provenance_file = "provenance.json"
        self.provenance = self._load_provenance()

        self.completed_steps = {
            step_entry["step"]
            for step_entry in self.provenance.get("steps", [])
            if step_entry.get("status") == "completed"
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

        for step_entry in self.provenance.get("steps", []):
            if step_entry["step"] == "calculate_lattice_parameters":
                self.a = step_entry["outputs"].get("lattice_constant", self.a)

            if step_entry["step"] == "generate_slab":
                slab_file = step_entry["outputs"].get("slab_file")
                if slab_file and os.path.exists(slab_file):
                    self.slab = read(slab_file)

            if step_entry["step"] == "calculate_alloy_lattice_parameters":
                self.a = step_entry["outputs"].get("lattice_constant", self.a)

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
                logger.info(f"Skipping completed step: {name}")
                continue

            logger.info(f"Running step: {name}")
            try:
                self._run_method(method)
                self.completed_steps.add(name)
            except Exception as e:
                logger.error(f"Step failed: {name}")
                logger.exception(e)
                self._record_step(name, status="failed", outputs={"error": str(e)})
                raise

        logger.info("Preprocessing workflow finished successfully")

    def _run_method(self, method):
        sig = inspect.signature(method)
        kwargs = {k: self.user_args[k] for k in sig.parameters if k in self.user_args}
        return method(**kwargs)

    # ========================================================
    # Preprocessing steps
    # ========================================================

    @step(2)
    @requires("optimize_lattice")
    def calculate_lattice_parameters(self, da=0.1):
        logger.info("Optimizing bulk lattice constant")

        calc = name_to_ase_software(self.software)(**self.lattice_opt_software_kwargs)

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

    @step(3)
    @requires("generate_slab")
    def generate_slab(self, a=None, slab_xyz="slab.xyz"):
        logger.info("Generating slab")

        if a is not None:
            self.a = a
        elif self.a is None:
            self.a = self.calculate_lattice_parameters()

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

        if self.software != "XTB":
            self.slab.calc = name_to_ase_software(self.software)(**self.software_kwargs)

            if hasattr(self, "freeze_bottom_n_layers"):
                self.freeze_bottom_n_layers()

            dyn = BFGSLineSearch(self.slab, trajectory="slab.traj")
            dyn.run(fmax=self.fmax, steps=200)

        write(slab_xyz, self.slab)

        self._record_step(
            "generate_slab",
            inputs={"surface_type": self.surface_type, "repeats": list(self.repeats), "vacuum": self.vacuum},
            outputs={"slab_file": slab_xyz},
        )

        return self.slab

    @step(4)
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

                logger.info(f"ecut={ecut} Ry | kpts={kpts} | E={E:.6f} eV")
                results.append((ecut, kpts, E))

        self._record_step(
            "run_convergence",
            inputs={"ecut_values": list(ecut_values), "kmesh_values": list(kmesh_values)},
            outputs={"results": results},
        )

        return results

    @step(5)
    @requires("run_adsorbate_convergence", "ecut_values", "kmesh_values")
    def run_convergence_adsorbate(
        self,
        ecut_values,
        kmesh_values,
        adsorbate=None,
        height=1.0,
        position="ontop",
    ):
        logger.info("Starting adsorbate convergence test")
        results = []

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

                slab = slab_clean.copy()
                slab.calc = calc
                E_slab = slab.get_potential_energy()

                slab_ads = slab_clean.copy()
                add_adsorbate(slab_ads, adsorbate, height=height, position=position)
                slab_ads.calc = calc
                E_slab_ads = slab_ads.get_potential_energy()

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

                logger.info(f"ecut={ecut} Ry | kpts={kpts} | E_ads={E_ads:.6f} eV")
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

    @step(6)
    @requires("analyze_sites")
    def analyze_sites(
        self,
        xyz_path="slab.xyz",
        n_layers=4,
        adsorbate_height=1.0,
        site_bond_cutoff=1.5,
        surface_string_no_defect="fcc332",
        tag_symbol="Ne",
        dz=0.1,
        stable_steps=3,
        max_drop=10.0,
        margin=0.25,
        min_clearance=1.3,
        save_all_drop_steps=True,
        verbose=True,
    ):
        """
        High-level: run site_analysis workflow_auto on the current slab.
        Writes outputs (traj/json) according to site_analysis defaults.
        """
        if self.slab is None:
            if not os.path.exists(xyz_path):
                raise ValueError("No slab in memory and xyz_path does not exist.")
            slab = read(xyz_path)
        else:
            slab = self.slab
            write(xyz_path, slab)

        # IMPORTANT: site_analysis expects nslab to be "number of slab atoms"
        # For a clean slab file, len(slab) is fine.
        nslab = len(slab)

        out = workflow_sites_auto(
            xyz_path=xyz_path,
            n_layers=n_layers,
            adsorbate_height=adsorbate_height,
            site_bond_cutoff=site_bond_cutoff,
            surface_string_no_defect=surface_string_no_defect,
            tag_symbol=tag_symbol,
            dz=dz,
            stable_steps=stable_steps,
            max_drop=max_drop,
            margin=margin,
            min_clearance=min_clearance,
            save_all_drop_steps=save_all_drop_steps,
            verbose=verbose,
        )

        self._record_step(
            "analyze_sites",
            inputs={
                "xyz_path": xyz_path,
                "n_layers": int(n_layers),
                "adsorbate_height": float(adsorbate_height),
                "site_bond_cutoff": float(site_bond_cutoff),
                "surface_string_no_defect": str(surface_string_no_defect),
                "tag_symbol": str(tag_symbol),
                "dz": float(dz),
                "stable_steps": int(stable_steps),
                "max_drop": float(max_drop),
                "margin": float(margin),
                "min_clearance": float(min_clearance),
                "save_all_drop_steps": bool(save_all_drop_steps),
            },
            outputs={"workflow": "workflow_sites_auto"},
        )

        return out