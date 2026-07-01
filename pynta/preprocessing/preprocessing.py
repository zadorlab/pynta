import os
import sys
import json
import time
import inspect
import logging
import copy

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy import optimize as opt

from ase import Atoms, build
from ase.io import read, write, Trajectory
from ase.build import add_adsorbate, bulk
from ase.optimize import BFGSLineSearch
from ase.calculators.espresso import Espresso
from ase.constraints import FixAtoms
from ase.visualize import view

from pynta.utils import name_to_ase_software
from pynta.tasks import *
from pynta.calculator import get_lattice_parameters
from pynta.main import generate_slab

from fireworks import LaunchPad, Workflow
from fireworks.queue.queue_launcher import rapidfire as rapidfirequeue
from fireworks.core.rocket_launcher import rapidfire
from fireworks.features.multi_launcher import launch_multiprocess
from fireworks.utilities.fw_serializers import load_object_from_file
from fireworks.core.fworker import FWorker

# site_analysis.py must be in the same directory as this file or in custom_surfaces/
_this_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_this_dir, "custom_surfaces"))
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
# Plotting
# ============================================================

def plot_results(cutoff_energy, kpts_list, energies):
    """Plot optimized energy vs kinetic energy cutoff and vs k-points."""

    def plot_energy_vs_cutoff():
        ecut_values = list(cutoff_energy.keys())
        optimized_energies = list(cutoff_energy.values())
        plt.figure(figsize=(10, 6))
        plt.plot(ecut_values, optimized_energies, marker='o', linestyle='-', color='b')
        plt.xlabel('Kinetic Energy Cutoff (eV)', fontsize=14)
        plt.ylabel('Optimized Energy (eV)', fontsize=14)
        plt.title('Optimized Energy vs Kinetic Energy Cutoff', fontsize=16)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig('optimized_energy_vs_cutoff.png')
        plt.close()

    def plot_energy_vs_kpoints():
        kpts_values = [kpts[0] * kpts[1] * kpts[2] for kpts in kpts_list]
        plt.figure(figsize=(10, 6))
        plt.plot(kpts_values, energies, marker='o', linestyle='-', color='r')
        plt.xlabel('K-points', fontsize=14)
        plt.ylabel('Optimized Energy (eV)', fontsize=14)
        plt.title('Optimized Energy vs K-points', fontsize=16)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig('optimized_energy_vs_kpoints.png')
        plt.close()

    plot_energy_vs_cutoff()
    plot_energy_vs_kpoints()


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
    """Minimize energy vs lattice constant a for a cubic cell by isotropically scaling the cell."""
    cell0 = atoms0.get_cell().array
    Lx0 = np.linalg.norm(cell0[0])

    if a_guess <= 0:
        raise ValueError("a_guess must be > 0")

    scx = max(int(round(Lx0 / a_guess)), 1)

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
    a0 = (Lx0 * s_opt) / scx

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
# High-level site-analysis workflow wrappers
# ============================================================

def workflow_sites_no_defect(
    slab, nslab, adsorbate_height=1.0, site_bond_cutoff=1.5,
    surface_string="fcc332", traj_filename="unique_sites.traj",
    sites_json="sites.json", neighbor_json="neighbor_site_list.json",
    sites_graph_json="sites_graph.json", tag_symbol="He", verbose=True,
):
    return sa.workflow_no_defect_unique_sites(
        slab=slab, nslab=nslab, adsorbate_height=adsorbate_height,
        site_bond_cutoff=site_bond_cutoff, surface_string=surface_string,
        traj_filename=traj_filename, sites_json=sites_json,
        neighbor_json=neighbor_json, sites_graph_json=sites_graph_json,
        tag_symbol=tag_symbol, verbose=verbose,
    )


def workflow_sites_defect_vacancy_drop(
    slab, nslab, xyz_path, surface_obj, site_bond_cutoff=1.5, tag_symbol="Ne",
    dz=0.1, stable_steps=3, max_drop=10.0, margin=0.25, min_clearance=1.3,
    save_all_drop_steps=True, traj_geom_all="geom_all.traj",
    traj_drop="drop_steps.traj", traj_maxcn="maxcn_geoms.xyz",
    json_geom_all_sites="geom_all_sites_lists.json",
    neighbor_json="neighbor_site_list.json", verbose=True,
):
    return sa.workflow_defect_vacancy_drop(
        slab=slab, nslab=nslab, xyz_path=xyz_path, surface_obj=surface_obj,
        site_bond_cutoff=site_bond_cutoff, tag_symbol=tag_symbol, dz=dz,
        stable_steps=stable_steps, max_drop=max_drop, margin=margin,
        min_clearance=min_clearance, save_all_drop_steps=save_all_drop_steps,
        traj_geom_all=traj_geom_all, traj_drop=traj_drop, traj_maxcn=traj_maxcn,
        json_geom_all_sites=json_geom_all_sites, neighbor_json=neighbor_json,
        verbose=verbose,
    )


def workflow_sites_auto(
    xyz_path, n_layers=4, adsorbate_height=1.0, site_bond_cutoff=1.5,
    surface_string_no_defect="fcc332", tag_symbol="Ne", dz=0.1, stable_steps=3,
    max_drop=10.0, margin=0.25, min_clearance=1.3, save_all_drop_steps=True, verbose=True,
):
    return sa.workflow_auto(
        xyz_path=xyz_path, n_layers=n_layers, adsorbate_height=adsorbate_height,
        site_bond_cutoff=site_bond_cutoff, surface_string_no_defect=surface_string_no_defect,
        tag_symbol=tag_symbol, dz=dz, stable_steps=stable_steps, max_drop=max_drop,
        margin=margin, min_clearance=min_clearance, save_all_drop_steps=save_all_drop_steps,
        verbose=verbose,
    )


# ============================================================
# Prep class
# ============================================================

class Prep:
    def __init__(
        self,
        # Surface / slab params
        metal="Pt",
        surface_type="fcc111",
        a0=3.96,
        c=None,
        repeats=(3, 3, 4),
        pbc=(True, True, False),
        # Calculation params
        software="Espresso",
        fmax=0.05,
        adsorbate=None,
        position=None,
        vacuum=10,
        slab=None,
        software_kwargs=None,
        lattice_opt_software_kwargs=None,
        path=None,
        label="prep",
        frozen_layers=None,
        ecut_range=None,
        # FireWorks params (all optional; FW is used only when launchpad_path is set or queue=True)
        launchpad_path=None,
        fworker_path=None,
        reset_launchpad=False,
        queue_adapter_path=None,
        queue=False,
        njobs_queue=0,
        num_jobs=25,
        # Alloy params
        alloy_path=None,
        alloy_A=None,
        alloy_B=None,
        alloy_xA=None,
        alloy_maxrep=6,
        alloy_scan_step=0.01,
        **kwargs,
    ):
        # ---- FireWorks setup (optional) ----
        self.qadapter = None
        if launchpad_path is not None or queue:
            self.launchpad = LaunchPad.from_file(launchpad_path) if launchpad_path else LaunchPad()
            if reset_launchpad:
                self.launchpad.reset('', require_password=False)
            self.fworker = FWorker.from_file(fworker_path) if fworker_path else FWorker()
        else:
            self.launchpad = None
            self.fworker = None

        self.queue = queue
        if queue and queue_adapter_path:
            self.qadapter = load_object_from_file(queue_adapter_path)

        self.njobs_queue = njobs_queue
        self.num_jobs = num_jobs

        # ---- Surface / slab ----
        self.metal = metal
        self.surface_type = surface_type
        self.a0 = a0
        self.a = a0
        self.c = c
        self.repeats = repeats
        self.pbc = pbc
        self.vacuum = vacuum
        self.layers = repeats[2]

        if slab is not None:
            if isinstance(slab, str):
                self.slab = read(slab)
            else:
                self.slab = slab
        else:
            self.slab = None

        self.nslab = len(self.slab) if self.slab is not None else int(np.prod(np.array(repeats)))
        self.frozen_layers = frozen_layers
        self.freeze_ind = (
            int((self.nslab / self.layers) * frozen_layers) if frozen_layers is not None else None
        )

        # ---- Calculation ----
        self.software = software
        self.fmax = fmax
        self.adsorbate = adsorbate
        self.position = position
        self.ecut_range = ecut_range
        self.path = path or os.getcwd()
        self.label = label

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
        else:
            self.software_kwargs = software_kwargs

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
        else:
            self.lattice_opt_software_kwargs = lattice_opt_software_kwargs

        # ---- Alloy ----
        self.alloy_path = alloy_path
        self.alloy_A = alloy_A
        self.alloy_B = alloy_B
        self.alloy_xA = alloy_xA
        self.alloy_maxrep = alloy_maxrep
        self.alloy_scan_step = alloy_scan_step

        # ---- Provenance / workflow state ----
        sig = inspect.signature(self.__init__)
        bound = sig.bind_partial(
            metal=metal, surface_type=surface_type, a0=a0, c=c, repeats=repeats,
            pbc=pbc, software=software, fmax=fmax, adsorbate=adsorbate, position=position,
            vacuum=vacuum, slab=slab, software_kwargs=software_kwargs,
            lattice_opt_software_kwargs=lattice_opt_software_kwargs, path=path,
            label=label, frozen_layers=frozen_layers, ecut_range=ecut_range,
            launchpad_path=launchpad_path, fworker_path=fworker_path,
            reset_launchpad=reset_launchpad, queue_adapter_path=queue_adapter_path,
            queue=queue, njobs_queue=njobs_queue, num_jobs=num_jobs,
            alloy_path=alloy_path, alloy_A=alloy_A, alloy_B=alloy_B, alloy_xA=alloy_xA,
            alloy_maxrep=alloy_maxrep, alloy_scan_step=alloy_scan_step,
        )
        self.user_args = {
            k: v for k, v in bound.arguments.items()
            if sig.parameters[k].default != v
        }
        self.user_args.update(kwargs)

        logger.info("User arguments:")
        for k, v in self.user_args.items():
            logger.info("  %s = %s", k, v)

        self.provenance_file = os.path.join(self.path, "provenance.json")
        self.provenance = self._load_provenance()
        self.completed_steps = {
            e["step"] for e in self.provenance.get("steps", [])
            if e.get("status") == "completed"
        }
        self.restore_state_from_provenance()

    # --------------------------------------------------------
    # FireWorks execution
    # --------------------------------------------------------

    def _require_fw(self):
        if self.launchpad is None:
            raise RuntimeError(
                "FireWorks is not configured. Pass launchpad_path or set queue=True."
            )

    def launch(self, single_job=False):
        """Run queued FireWorks rockets (queue, single, or multi-process)."""
        self._require_fw()
        if self.queue:
            rapidfirequeue(self.launchpad, self.fworker, self.qadapter,
                           njobs_queue=self.njobs_queue, nlaunches="infinite")
        elif self.num_jobs == 1 or single_job:
            rapidfire(self.launchpad, self.fworker, nlaunches="infinite")
        else:
            launch_multiprocess(self.launchpad, self.fworker, "INFO",
                                "infinite", self.num_jobs, 5)

    def lattice_optimization(self, da=0.1, scan_step=0.02, centered=True,
                             inplane_only=None, n_fit=None, label=None,
                             skip_launch=False, wait=True, poll=2.0):
        """Run lattice-constant scan as a FireWorks workflow (parallel fan-out)."""
        self._require_fw()
        if self.slab is None:
            raise ValueError("Prep.slab is None; set a slab before lattice_optimization().")

        lbl = label or self.label or "lattice"
        software_kwargs = self.lattice_opt_software_kwargs or self.software_kwargs
        scan_dir = os.path.join(self.path, "lattice_scan")
        out_path = os.path.join(scan_dir, "lattice_constant.json")

        wf = lattice_optimization_workflow(
            slab=self.slab, a0=self.a0, software=self.software,
            software_kwargs=software_kwargs, da=da, scan_step=scan_step,
            centered=centered, inplane_only=inplane_only, pbc=self.pbc,
            label=lbl, scan_dir=scan_dir, out_path=out_path, n_fit=n_fit,
        )
        self.launchpad.add_wf(wf)
        logger.info("Added lattice-opt workflow '%s' (%d energy FWs + 1 collector)",
                    wf.name, len(wf.fws) - 1)

        if skip_launch:
            return wf

        self.launch(single_job=True)

        if not wait:
            return wf

        while not os.path.exists(out_path):
            time.sleep(poll)
        with open(out_path) as f:
            result = json.load(f)
        a_opt = result["lattice_constant"]
        self.a = a_opt
        if result.get("edge_of_window_warning"):
            logger.warning(
                "Fitted a=%.4f is at/outside the scan window %s; widen da%s.",
                a_opt, result.get("scan_window"),
                " or set centered=True" if not centered else "")
        logger.info("Optimized lattice constant: a = %.6f Angstrom", a_opt)
        return a_opt

    def opt_cutoff_energy(self, ecut_range=None, skip_launch=False):
        """Optimize cutoff energy via a FireWorks workflow (or locally if skip_launch)."""
        self._require_fw()

        if self.slab is None:
            self.slab = bulk(self.metal, self.surface_type[:3], a=self.a0)
            self.slab.pbc = (True, True, False)

        write(os.path.join(self.path, "test_bulk.xyz"), self.slab)
        ecut_range = ecut_range or self.ecut_range
        cutoff_energy = {}

        for ecut in range(*ecut_range):
            sw_kw = copy.deepcopy(self.software_kwargs)
            if self.software.lower() == 'gpaw':
                sw_kw['ecut'] = ecut
            elif self.software.lower() == 'espresso':
                sw_kw['ecutwfc'] = ecut
            elif self.software.lower() in ['nwchem', 'pwdft']:
                sw_kw['cutoff'] = ecut

            if skip_launch:
                calc = name_to_ase_software(self.software)(**sw_kw)
                self.slab.calc = calc
                E = self.slab.get_potential_energy()
                cutoff_energy[ecut] = E
                print(f'Calculator: {self.software}, Cutoff: {ecut} eV, Energy: {E:.6f} eV')
            else:
                out_json = os.path.join(self.path, f"energy_convergence_{ecut}.json")
                fwenergy = energy_firework(
                    os.path.join(self.path, "test_bulk.xyz"), self.software,
                    f"energy_convergence_{ecut}", software_kwargs=sw_kw, out_path=out_json,
                )
                wf = Workflow([fwenergy], name=f"{self.label}_ecut_{ecut}")
                self.launchpad.add_wf(wf)
                while not os.path.exists(out_json):
                    time.sleep(1)
                with open(out_json) as f:
                    E = json.load(f)
                cutoff_energy[ecut] = E
                print(f'Calculator: {self.software}, Cutoff: {ecut} eV, Energy: {E:.6f} eV')

        with open(os.path.join(self.path, 'cutoff_energy.json'), 'w') as f:
            json.dump(cutoff_energy, f, indent=4)

        return cutoff_energy

    def opt_kpoints(self, kpts_range=None, skip_launch=False):
        """Optimize k-points via a FireWorks workflow (or locally if skip_launch)."""
        self._require_fw()

        if self.slab is None:
            self.slab = self.make_test_slab()

        energies = []
        kpts_list = []
        results = {}

        for kpts in kpts_range:
            sw_kw = copy.deepcopy(self.software_kwargs)
            sw_kw['kpts'] = kpts

            if skip_launch:
                calc = name_to_ase_software(self.software)(**sw_kw)
                slab = self.slab.copy()
                slab.calc = calc
                slab.pbc = (True, True, False)
                E = slab.get_potential_energy()
                kpts_list.append(kpts)
                energies.append(E)
                results[str(kpts)] = E
            else:
                kpts_name = ''.join(map(str, kpts))
                out_json = os.path.join(self.path, f"kpoints_{kpts_name}_energy.json")
                fwkpt = energy_firework(
                    os.path.join(self.path, "test_slab.xyz"), self.software,
                    f"kpoints_{kpts_name}", software_kwargs=sw_kw,
                )
                wf = Workflow([fwkpt], name=f"{self.label}_kpt_{kpts_name}")
                self.launchpad.add_wf(wf)
                while not os.path.exists(out_json):
                    time.sleep(1)
                with open(out_json) as f:
                    E = json.load(f)
                kpts_list.append(kpts)
                energies.append(E)
                results[str(kpts)] = E

        with open(os.path.join(self.path, 'kpts.json'), 'w') as f:
            json.dump(results, f, indent=4)

        return kpts_list, energies

    def make_test_slab(self):
        """Create a small test slab for convergence calculations."""
        slab_builder = getattr(build, self.surface_type)
        test_slab = slab_builder(self.metal, size=self.repeats)
        test_slab.center(vacuum=self.vacuum, axis=2)
        mask = [atom.tag > 2 for atom in test_slab]
        test_slab.set_constraint(FixAtoms(mask=mask))
        if self.adsorbate is not None:
            add_adsorbate(test_slab, adsorbate=self.adsorbate, position=self.position)
        write(os.path.join(self.path, "test_slab.xyz"), test_slab)
        return test_slab

    # --------------------------------------------------------
    # Local (non-FireWorks) execution
    # --------------------------------------------------------

    @step(2)
    @requires("optimize_lattice")
    def calculate_lattice_parameters(self, da=0.1):
        """Optimize bulk lattice constant locally via energy scan + quadratic fit."""
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
        logger.info("Optimized lattice constant: a = %.6f Å", self.a)
        self._record_step(
            "calculate_lattice_parameters",
            inputs={"a0": self.a0, "da": da},
            outputs={"lattice_constant": self.a},
        )
        return self.a

    @step(3)
    @requires("generate_slab")
    def generate_slab(self, a=None, slab_xyz="slab.xyz"):
        """Build and locally optimize the slab structure."""
        logger.info("Generating slab")
        if a is not None:
            self.a = a
        elif self.a is None:
            self.a = self.calculate_lattice_parameters()

        slab_builder = getattr(build, self.surface_type)
        if self.c is not None:
            slab = slab_builder(symbol=self.metal, size=self.repeats, a=self.a,
                                c=self.c, vacuum=self.vacuum)
        else:
            slab = slab_builder(symbol=self.metal, size=self.repeats, a=self.a,
                                vacuum=self.vacuum)
        slab.pbc = self.pbc
        self.slab = slab

        if self.software != "XTB":
            self.slab.calc = name_to_ase_software(self.software)(**self.software_kwargs)
            dyn = BFGSLineSearch(self.slab, trajectory="slab.traj")
            dyn.run(fmax=self.fmax, steps=200)

        write(slab_xyz, self.slab)
        self._record_step(
            "generate_slab",
            inputs={"surface_type": self.surface_type, "repeats": list(self.repeats),
                    "vacuum": self.vacuum},
            outputs={"slab_file": slab_xyz},
        )
        return self.slab

    @step(4)
    @requires("ecut_values", "kmesh_values")
    def run_convergence(self, ecut_values, kmesh_values):
        """Run local ecut + kpoints convergence grid on self.slab."""
        logger.info("Running convergence tests")
        results = []
        for ecut in ecut_values:
            for kpts in kmesh_values:
                sw_kw = copy.deepcopy(self.software_kwargs)
                sw_kw["ecutwfc"] = ecut
                sw_kw["ecutrho"] = 4 * ecut
                calc = Espresso(
                    input_data=sw_kw,
                    pseudopotentials=sw_kw["pseudopotentials"],
                    kpts=kpts,
                )
                self.slab.calc = calc
                E = self.slab.get_potential_energy()
                logger.info("ecut=%s Ry | kpts=%s | E=%.6f eV", ecut, kpts, E)
                results.append((ecut, kpts, E))
        self._record_step(
            "run_convergence",
            inputs={"ecut_values": list(ecut_values), "kmesh_values": list(kmesh_values)},
            outputs={"results": results},
        )
        return results

    @step(5)
    @requires("run_adsorbate_convergence", "ecut_values", "kmesh_values")
    def run_convergence_adsorbate(self, ecut_values, kmesh_values,
                                  adsorbate=None, height=1.0, position="ontop"):
        """Run local adsorbate adsorption-energy convergence on self.slab."""
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
                sw_kw = copy.deepcopy(self.software_kwargs)
                sw_kw["ecutwfc"] = ecut
                sw_kw["ecutrho"] = 4 * ecut
                calc = Espresso(input_data=sw_kw,
                                pseudopotentials=sw_kw["pseudopotentials"], kpts=kpts)

                slab = slab_clean.copy()
                slab.calc = calc
                E_slab = slab.get_potential_energy()

                slab_ads = slab_clean.copy()
                add_adsorbate(slab_ads, adsorbate, height=height, position=position)
                slab_ads.calc = calc
                E_slab_ads = slab_ads.get_potential_energy()

                if is_hydrogen:
                    ref = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]],
                                cell=[10, 10, 10], pbc=False)
                    ref.calc = calc
                    E_ref = 0.5 * ref.get_potential_energy()
                else:
                    ref = adsorbate.copy()
                    ref.cell = [10, 10, 10]
                    ref.pbc = False
                    ref.calc = calc
                    E_ref = ref.get_potential_energy()

                E_ads = E_slab_ads - E_slab - E_ref
                logger.info("ecut=%s Ry | kpts=%s | E_ads=%.6f eV", ecut, kpts, E_ads)
                results.append((ecut, kpts, E_ads))

        self._record_step(
            "run_convergence_adsorbate",
            inputs={"ecut_values": list(ecut_values), "kmesh_values": list(kmesh_values),
                    "adsorbate": str(adsorbate), "height": height, "position": position},
            outputs={"results": results},
        )
        return results

    @step(6)
    @requires("analyze_sites")
    def analyze_sites(self, xyz_path="slab.xyz", n_layers=4, adsorbate_height=1.0,
                      site_bond_cutoff=1.5, surface_string_no_defect="fcc332",
                      tag_symbol="Ne", dz=0.1, stable_steps=3, max_drop=10.0,
                      margin=0.25, min_clearance=1.3, save_all_drop_steps=True, verbose=True):
        """Run site_analysis workflow_auto on the current slab."""
        if self.slab is None:
            if not os.path.exists(xyz_path):
                raise ValueError("No slab in memory and xyz_path does not exist.")
            slab = read(xyz_path)
        else:
            slab = self.slab
            write(xyz_path, slab)

        out = workflow_sites_auto(
            xyz_path=xyz_path, n_layers=n_layers, adsorbate_height=adsorbate_height,
            site_bond_cutoff=site_bond_cutoff, surface_string_no_defect=surface_string_no_defect,
            tag_symbol=tag_symbol, dz=dz, stable_steps=stable_steps, max_drop=max_drop,
            margin=margin, min_clearance=min_clearance,
            save_all_drop_steps=save_all_drop_steps, verbose=verbose,
        )
        self._record_step(
            "analyze_sites",
            inputs={"xyz_path": xyz_path, "n_layers": int(n_layers),
                    "adsorbate_height": float(adsorbate_height),
                    "site_bond_cutoff": float(site_bond_cutoff),
                    "surface_string_no_defect": str(surface_string_no_defect),
                    "tag_symbol": str(tag_symbol)},
            outputs={"workflow": "workflow_sites_auto"},
        )
        return out

    # --------------------------------------------------------
    # Workflow auto-dispatch
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
        """Auto-discover and run all workflow steps whose required args are present."""
        methods = self.discover_runnable_methods()
        methods.sort(key=lambda m: getattr(m, "_order", 999))
        logger.info("Workflow steps: %s", [m.__name__ for m in methods])
        for method in methods:
            name = method.__name__
            if name in self.completed_steps and not force:
                logger.info("Skipping completed step: %s", name)
                continue
            logger.info("Running step: %s", name)
            try:
                sig = inspect.signature(method)
                kwargs = {k: self.user_args[k] for k in sig.parameters if k in self.user_args}
                method(**kwargs)
                self.completed_steps.add(name)
            except Exception as e:
                logger.error("Step failed: %s", name)
                logger.exception(e)
                self._record_step(name, status="failed", outputs={"error": str(e)})
                raise
        logger.info("Preprocessing workflow finished successfully")

    # --------------------------------------------------------
    # Provenance
    # --------------------------------------------------------

    def _load_provenance(self):
        if os.path.exists(self.provenance_file):
            logger.info("Loading existing provenance.json")
            with open(self.provenance_file) as f:
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
        self.provenance["steps"].append({
            "step": name,
            "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "status": status,
            "inputs": inputs or {},
            "outputs": outputs or {},
        })
        self._save_provenance()

    def restore_state_from_provenance(self):
        logger.info("Restoring state from provenance")
        for entry in self.provenance.get("steps", []):
            if entry["step"] == "calculate_lattice_parameters":
                self.a = entry["outputs"].get("lattice_constant", self.a)
            if entry["step"] == "generate_slab":
                slab_file = entry["outputs"].get("slab_file")
                if slab_file and os.path.exists(slab_file):
                    self.slab = read(slab_file)
            if entry["step"] == "calculate_alloy_lattice_parameters":
                self.a = entry["outputs"].get("lattice_constant", self.a)
