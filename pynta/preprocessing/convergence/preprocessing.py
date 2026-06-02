import os
import json
import time
import inspect
import logging
import numpy as np
import yaml
import copy
from copy import deepcopy

from scipy import optimize as opt

from ase import Atoms, build
from ase.io import read, write, Trajectory
from ase.build import add_adsorbate, bulk
from ase.constraints import FixAtoms
from ase.optimize import BFGSLineSearch
from ase.data import chemical_symbols, reference_states
from ase.geometry import get_distances
from ase.geometry.analysis import Analysis
from ase.neighborlist import NeighborList

from pynta.utils import name_to_ase_software, get_unique_sym_struct_index_clusters
from pynta.mol import add_adsorbate_to_site

from molecule.molecule import Molecule, Atom, Bond


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
# Postprocessing / plotting utilities
# ============================================================

def plot_lattice_constant(lout, title="Lattice constant optimization"):
    import matplotlib.pyplot as plt

    a_opt = lout["lattice_constant"]
    avals = np.array(lout["scan_avals_A"])
    Evals = np.array(lout["scan_Evals_eV"])

    inds = np.argsort(Evals)[:7]
    p = np.polyfit(avals[inds], Evals[inds], 2)
    x_fit = np.linspace(avals.min(), avals.max(), 300)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(avals, Evals, color="orange", zorder=3, label="DFT scan")
    ax.plot(x_fit, np.polyval(p, x_fit), color="steelblue", label="Quadratic fit")
    ax.axvline(a_opt, color="red", linestyle="--", linewidth=1.5,
               label=f"a$_{{opt}}$ = {a_opt:.4f} Å")
    ax.set_xlabel("Lattice constant $a$ (Å)")
    ax.set_ylabel("Total energy (eV)")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()

    print(f"Optimized lattice constant: a = {a_opt:.6f} Å")
    return a_opt


def plot_ecut_convergence(results, ref_kpts=None):
    import matplotlib.pyplot as plt

    ecuts_all = sorted(set(r[0] for r in results))
    kpts_all = list(dict.fromkeys(tuple(r[1]) for r in results))
    if ref_kpts is None:
        ref_kpts = kpts_all[-1]

    ecut_data = sorted(
        [(r[0], r[2]) for r in results if tuple(r[1]) == tuple(ref_kpts)],
        key=lambda x: x[0],
    )
    ecut_x = [d[0] for d in ecut_data]
    ecut_E = np.array([d[1] for d in ecut_data])
    ecut_E_rel = np.abs(ecut_E - ecut_E[-1])
    kpts_label = "×".join(str(k) for k in ref_kpts)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].plot(ecut_x, ecut_E, marker="o", color="steelblue")
    axes[0].set_xlabel("Kinetic energy cutoff (Ry)")
    axes[0].set_ylabel("Energy (eV)")
    axes[0].set_title(f"ecut convergence  [kpts = {kpts_label}]")
    axes[0].grid(True, alpha=0.4)

    axes[1].plot(ecut_x, ecut_E_rel, marker="o", color="orange")
    axes[1].set_xlabel("Kinetic energy cutoff (Ry)")
    axes[1].set_ylabel("|ΔE| vs finest ecut (eV)")
    axes[1].set_title("ecut convergence (relative)")
    axes[1].set_yscale("log")
    axes[1].grid(True, which="both", alpha=0.4)

    plt.tight_layout()
    plt.show()

    for ec, E, dE in zip(ecut_x, ecut_E, ecut_E_rel):
        print(f"  ecut = {ec:3.0f} Ry  |  E = {E:.6f} eV  |  |ΔE| = {dE:.2e} eV")

    return ecuts_all, kpts_all


def plot_kpoint_convergence(results, ref_ecut=None):
    import matplotlib.pyplot as plt

    ecuts_all = sorted(set(r[0] for r in results))
    if ref_ecut is None:
        ref_ecut = ecuts_all[-1]

    kpt_data = sorted(
        [(tuple(r[1]), r[2]) for r in results if r[0] == ref_ecut],
        key=lambda x: x[0][0],
    )
    kpt_labels = ["×".join(str(k) for k in kp) for kp, _ in kpt_data]
    kpt_E = np.array([d[1] for d in kpt_data])
    kpt_E_rel = np.abs(kpt_E - kpt_E[-1])
    x_idx = range(len(kpt_labels))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].plot(x_idx, kpt_E, marker="o", color="steelblue")
    axes[0].set_xticks(x_idx)
    axes[0].set_xticklabels(kpt_labels, rotation=30)
    axes[0].set_xlabel("k-point mesh")
    axes[0].set_ylabel("Energy (eV)")
    axes[0].set_title(f"k-point convergence  [ecut = {ref_ecut} Ry]")
    axes[0].grid(True, alpha=0.4)

    axes[1].plot(x_idx, kpt_E_rel, marker="o", color="orange")
    axes[1].set_xticks(x_idx)
    axes[1].set_xticklabels(kpt_labels, rotation=30)
    axes[1].set_xlabel("k-point mesh")
    axes[1].set_ylabel("|ΔE| vs finest mesh (eV)")
    axes[1].set_title("k-point convergence (relative)")
    axes[1].set_yscale("log")
    axes[1].grid(True, which="both", alpha=0.4)

    plt.tight_layout()
    plt.show()

    for kl, E, dE in zip(kpt_labels, kpt_E, kpt_E_rel):
        print(f"  kpts = {kl:9s}  |  E = {E:.6f} eV  |  |ΔE| = {dE:.2e} eV")


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
    cell0 = atoms0.get_cell().array
    Lx0 = np.linalg.norm(cell0[0])
    scx = int(round(Lx0 / a_guess))
    scx = max(scx, 1)

    def energy_at_a(a_target):
        s = (a_target * scx) / Lx0
        atoms = atoms0.copy()
        atoms.set_cell(atoms.get_cell() * s, scale_atoms=True)
        atoms.calc = calc
        return atoms.get_potential_energy()

    avals = np.arange(a_guess - da, a_guess + da + 1e-12, scan_step)
    Evals = np.array([energy_at_a(a) for a in avals])

    inds = np.argsort(Evals)[:7] if len(Evals) >= 7 else np.argsort(Evals)
    p = np.polyfit(avals[inds], Evals[inds], 2)
    if p[0] <= 0:
        raise RuntimeError("Quadratic fit curvature is non-positive; increase scan window/points.")
    a_est = -p[1] / (2.0 * p[0])

    out = opt.minimize_scalar(
        energy_at_a,
        method="bounded",
        bounds=(a_est - 0.02, a_est + 0.02),
        options={"xatol": 1e-4},
    )

    return {
        "a0_A": float(out.x),
        "Emin_eV": float(out.fun),
        "success": bool(out.success),
        "message": out.message,
        "scan_avals_A": avals.tolist(),
        "scan_Evals_eV": Evals.tolist(),
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
# Site / geometry utilities (from site_analysis.py)
# ============================================================

def sites_match(site1, site2, slab, tol=0.5):
    _, dist = get_distances(
        [site1["position"]],
        [site2["position"]],
        cell=slab.cell,
        pbc=slab.pbc
    )
    return (
        dist < tol
        and site1["site"] == site2["site"]
        and site1["morphology"] == site2["morphology"]
    )


def get_occupied_sites(struct, sites, nslab, cutoff):
    occ = []
    for i in range(nslab, len(struct)):
        pos = struct.positions[i]
        best_site, best_dist = None, None
        for site in sites:
            _, dist = get_distances(
                [site["position"]],
                [pos],
                cell=struct.cell,
                pbc=struct.pbc
            )
            if best_dist is None or dist < best_dist:
                best_site, best_dist = site, dist
        if best_dist and best_dist < cutoff:
            s = dict(best_site)
            s["bonding_index"] = i
            s["bond_length"] = best_dist
            occ.append(s)
    return occ


def generate_unique_sites(slab, sites, nslab, site_bond_cutoff, adsorbate_height):
    occ = get_occupied_sites(slab, sites, nslab, site_bond_cutoff)
    unocc = [s for s in sites if not any(sites_match(s, o, slab) for o in occ)]

    geoms = []
    sites_per_geom = []

    for site in unocc:
        g = slab.copy()
        add_adsorbate_to_site(
            g,
            adsorbate=Atoms("He"),
            surf_ind=0,
            site=site,
            height=adsorbate_height,
            tilt_angle=25.24,
            offset=False
        )
        geoms.append(g)
        sites_per_geom.append([site])

    clusters = get_unique_sym_struct_index_clusters(geoms)

    return (
        [geoms[c[0]] for c in clusters],
        [sites_per_geom[c[0]] for c in clusters]
    )


def write_trajectory_pynta(
    slab, cas, nslab, site_bond_cutoff, adsorbate_height,
    trajectory_filename="unique_sites.traj"
):
    single_geoms, single_sites_lists = generate_unique_sites(
        slab,
        cas.get_sites(),
        nslab,
        site_bond_cutoff,
        adsorbate_height
    )

    print(f'There are {len(single_sites_lists)} unique sites out of {len(cas.get_sites())}.')

    traj = Trajectory(trajectory_filename, "w")
    all_sites = cas.get_sites()
    for g in single_geoms:
        traj.write(g)
    traj.close()

    save_sites_to_json(all_sites)


def save_sites_to_json(single_sites_lists, filename='sites.json'):
    def process_data(data):
        if isinstance(data, dict):
            return {key: process_data(value) for key, value in data.items()}
        elif isinstance(data, list):
            return [process_data(item) for item in data]
        elif isinstance(data, np.ndarray):
            return data.tolist()
        elif isinstance(data, (np.int64, np.float64)):
            return data.item()
        elif data is None:
            return "null"
        else:
            return data

    if isinstance(single_sites_lists, list) and all(isinstance(item, list) for item in single_sites_lists):
        processed_sites_data = process_data([item for sublist in single_sites_lists for item in sublist])
    else:
        processed_sites_data = process_data(single_sites_lists)

    with open(filename, "w") as f:
        json.dump(processed_sites_data, f, indent=4)

    print(f"Sites data saved to '{filename}' with None replaced by 'null'.")


def write_trajectory_for_acat(slab, cas, trajectory_filename):
    traj = Trajectory(trajectory_filename, 'w')
    for _, site in enumerate(cas.get_unique_sites()):
        my_slab = copy.deepcopy(slab)
        new_position = Atom('He', site["position"])
        my_slab.append(new_position)
        traj.write(my_slab)
    traj.close()


def save_neighbor_site_list_to_json(cas, filename='neighbor_site_list.json'):
    def convert_numpy_types(data):
        if isinstance(data, dict):
            return {key: convert_numpy_types(value) for key, value in data.items()}
        elif isinstance(data, list):
            return [convert_numpy_types(item) for item in data]
        elif isinstance(data, np.ndarray):
            return data.tolist()
        elif isinstance(data, (np.int64, np.float64)):
            return data.item()
        else:
            return data

    neighbor_site_list = cas.get_neighbor_site_list()
    neighbor_site_list = convert_numpy_types(neighbor_site_list)

    with open(filename, 'w') as f:
        json.dump(neighbor_site_list, f, indent=4)

    print(f"Neighbor site list saved to '{filename}'.")


def ase_to_rmg_symbol(sym):
    if sym in {"H", "C", "N", "O", "S", "F", "Cl", "Br", "I", "P", "Si", "Li"}:
        return sym
    if sym == "He":
        return "Li"
    return "X"


def lone_pairs_for(sym):
    return {"N": 1, "P": 1, "O": 2, "S": 2, "F": 3, "Cl": 3, "Br": 3, "I": 3}.get(sym, 0)


def generate_graph_slab_fixed(atoms):
    analysis = Analysis(atoms)
    adj = analysis.adjacency_matrix[0]

    adatoms = []
    for at in atoms:
        rsym = ase_to_rmg_symbol(at.symbol)
        adatoms.append(Atom(element=rsym, lone_pairs=lone_pairs_for(rsym)))

    mol = Molecule(atoms=adatoms)

    for i in range(len(adatoms)):
        for j in range(i + 1, len(adatoms)):
            if adj[i, j]:
                mol.add_bond(Bond(mol.atoms[i], mol.atoms[j], 1.0))

    mol.update_multiplicity()
    mol.update_connectivity_values()
    return mol


def classify_all_sites(single_geoms, single_sites_lists):
    admols = []
    geom_indices = []

    for i, (geom, sites) in enumerate(zip(single_geoms, single_sites_lists)):
        for s in sites:
            if s["site"]:
                admols.append(generate_graph_slab_fixed(geom))
                geom_indices.append(i)
                break

    return admols, geom_indices


def cluster_isomorphic_graphs(admols):
    n = len(admols)
    iso_mat = np.zeros((n, n), dtype=bool)

    for i in range(n):
        iso_mat[i, i] = True
        for j in range(i + 1, n):
            iso = admols[i].is_isomorphic(admols[j], strict=True)
            iso_mat[i, j] = iso_mat[j, i] = iso

    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(n):
        for j in range(i + 1, n):
            if iso_mat[i, j]:
                union(i, j)

    clusters = {}
    for i in range(n):
        clusters.setdefault(find(i), []).append(i)

    return iso_mat, clusters


def update_threefold_site_labels(single_sites_lists, clusters, geom_indices):
    type_map = {}

    for type_id, members in enumerate(clusters.values(), start=1):
        label = f"site{type_id}"
        for graph_idx in members:
            type_map[geom_indices[graph_idx]] = label

    for geom_idx, sites in enumerate(single_sites_lists):
        for s in sites:
            if s["site"] == "3fold":
                s["site"] = type_map.get(geom_idx, "3fold-unknown")

    return type_map


def update_site_labels(single_sites_lists, clusters, trajectory_filename="unique_site_graph.traj"):
    site_label_mapping = {}
    updated_sites_lists = []

    print("Clusters structure:", clusters)
    print("Single sites lists structure:", single_sites_lists)

    for cluster_index, members in clusters.items():
        if not isinstance(members, list):
            print(f"Warning: Expected a list for members in cluster {cluster_index}, got {type(members)}")
            continue

        first_member_index = members[0]
        representative_site = single_sites_lists[first_member_index][0]["site"]
        site_label_mapping[cluster_index] = representative_site

        updated_sites_lists.append(single_sites_lists[first_member_index])

        for member in members:
            single_sites_lists[member][0]["site"] = representative_site

    return updated_sites_lists, site_label_mapping


def write_sites_json(single_sites_lists, clusters, filename="sites_graph.json"):
    json_data = {
        "n_unique_geometries": len(single_sites_lists),
        "n_distinct_three_fold_types": len(clusters),
        "geometries": []
    }

    for sites in single_sites_lists:
        for s in sites:
            entry = {
                "site": s["site"],
                "surface": s.get("surface"),
                "morphology": s.get("morphology"),
                "position": list(map(float, s.get("position", []))),
                "normal": list(map(float, s.get("normal", []))),
                "indices": list(map(int, s.get("indices", []))),
                "composition": s.get("composition", "null"),
                "subsurf_index": s.get("subsurf_index", "null"),
                "subsurf_element": s.get("subsurf_element", "null"),
                "label": s.get("label", "null"),
                "topology": list(map(int, s.get("topology", [])))
            }
            json_data["geometries"].append(entry)

    with open(filename, "w") as f:
        json.dump(json_data, f, indent=4)


def write_trajectory_graph(updated_sites_lists, clusters, trajectory_filename="unique_sites_graph.traj"):
    print(f'There are {len(updated_sites_lists)} unique sites out of {len(clusters)}.')

    traj = Trajectory(trajectory_filename, "w")

    for site in updated_sites_lists:
        atoms_data = site[0]
        atoms = Atoms(
            symbols=atoms_data['composition'],
            positions=atoms_data['position'],
        )
        traj.write(atoms)

    traj.close()


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
        self.user_args["software"] = software

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

        self.output_file = "prep_convergence.json"
        self.output = self._load_output()

        self.completed_steps = {
            step["step"]
            for step in self.output.get("steps", [])
            if step.get("status") == "completed"
        }

        self.restore_state_from_output()

    # --------------------------------------------------------
    # Output helpers
    # --------------------------------------------------------

    def _load_output(self):
        if os.path.exists(self.output_file):
            logger.info("Loading existing output.json")
            with open(self.output_file, "r") as f:
                return json.load(f)

        return {
            "created": time.strftime("%Y-%m-%d %H:%M:%S"),
            "user_arguments": self.user_args,
            "steps": [],
            "outputs": [],
        }

    def _save_output(self):
        with open(self.output_file, "w") as f:
            json.dump(self.output, f, indent=2)

    def _record_step(self, name, inputs=None, outputs=None, status="completed"):
        entry = {
            "step": name,
            "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "status": status,
            "inputs": inputs or {},
            "outputs": outputs or {},
        }
        self.output["steps"].append(entry)
        self._save_output()

    # --------------------------------------------------------
    # State restoration
    # --------------------------------------------------------

    def restore_state_from_output(self):
        logger.info("Restoring state from output")

        for step_entry in self.output.get("steps", []):
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
                self._record_step(name, status="failed", outputs={"error": str(e)})
                raise

        logger.info("Preprocessing workflow finished successfully")

    def _run_method(self, method):
        sig = inspect.signature(method)
        kwargs = {k: self.user_args[k] for k in sig.parameters if k in self.user_args}
        return method(**kwargs)

    # ========================================================
    # Auxiliary functions
    # ========================================================

    def apply_constraints(self, sp, constraints):
        out_constraints = []

        for c in constraints:
            if c == "freeze half slab":
                z_mid = sp.cell[2, 2] / 2.0
                out_constraints.append(
                    FixAtoms(indices=[atom.index for atom in sp if atom.position[2] < z_mid])
                )

            elif c.startswith("freeze all"):
                sym = c.split()[2]
                out_constraints.append(
                    FixAtoms(indices=[atom.index for atom in sp if atom.symbol == sym])
                )

            elif c == "freeze bottom layers":
                if self.frozen_layers is None:
                    raise ValueError("'freeze bottom layers' requested but frozen_layers is not set")

                nslab = len(sp)
                layers = getattr(self, "layers", None) or self.repeats[2]

                if self.frozen_layers > layers:
                    raise ValueError(
                        f"frozen_layers ({self.frozen_layers}) > total layers ({layers})"
                    )

                n_per_layer = nslab // layers
                n_freeze = n_per_layer * self.frozen_layers
                z_sorted_indices = sorted(range(nslab), key=lambda i: sp.positions[i][2])
                out_constraints.append(FixAtoms(indices=z_sorted_indices[:n_freeze]))

            else:
                raise ValueError(f"Unknown constraint: {c}")

        if out_constraints:
            sp.set_constraint(out_constraints)

    # ========================================================
    # Actual preprocessing steps
    # ========================================================

    @step(1)
    @requires("optimize_alloy_lattice")
    def calculate_alloy_lattice_parameters(self, da=0.1):
        logger.info("Optimizing ALLOY bulk lattice constant")

        calc = name_to_ase_software(self.software)(**self.lattice_opt_software_kwargs)

        provenance_inputs = {
            "surface_type": self.surface_type,
            "a_guess": self.a0,
            "da": da,
            "scan_step": self.alloy_scan_step,
        }

        if self.alloy_path is not None:
            if not os.path.exists(self.alloy_path):
                raise FileNotFoundError(f"alloy_path not found: {self.alloy_path}")

            logger.info(f"Reading alloy structure from alloy_path: {self.alloy_path}")
            atoms0 = read(self.alloy_path)
            atoms0.set_pbc(True)

            provenance_inputs.update({"alloy_path": self.alloy_path})

            optres = optimize_a_by_energy_minimization(
                atoms0,
                calc=calc,
                a_guess=self.a0,
                da=da,
                scan_step=self.alloy_scan_step,
            )

            res = {
                "source": "alloy_path",
                "supercell_repeat": None,
                "counts": None,
                **optres
            }

        else:
            if self.alloy_A is None or self.alloy_B is None or self.alloy_xA is None:
                raise ValueError(
                    "alloy_path is None, so alloy_A, alloy_B, and alloy_xA must be provided."
                )

            logger.info("Building ordered alloy cell (minimal exact composition)")
            crystal = infer_crystal_from_surface(self.surface_type, crystal=None)
            if crystal not in ("fcc", "bcc"):
                raise ValueError("Only fcc/bcc are supported for ordered alloy building.")

            atoms0, sc, counts = build_ordered_alloy_cubic(
                self.alloy_A, self.alloy_B, self.alloy_xA,
                crystal=crystal, a_guess=self.a0, maxrep=self.alloy_maxrep
            )

            provenance_inputs.update({
                "A": self.alloy_A,
                "B": self.alloy_B,
                "xA": self.alloy_xA,
                "maxrep": self.alloy_maxrep,
                "crystal": crystal,
            })

            optres = optimize_a_by_energy_minimization(
                atoms0,
                calc=calc,
                a_guess=self.a0,
                da=da,
                scan_step=self.alloy_scan_step,
            )

            res = {
                "source": "ordered_build",
                "crystal": crystal,
                "supercell_repeat": sc,
                "counts": {"nA": counts[0], "nB": counts[1], "N": counts[2]},
                **optres
            }

        if not res["success"]:
            raise RuntimeError(f"Alloy lattice optimization failed: {res['message']}")

        self.a = res["a0_A"]
        logger.info(f"Optimized alloy lattice constant: a = {self.a:.6f} Å (source={res.get('source')})")

        self._record_step(
            "calculate_alloy_lattice_parameters",
            inputs=provenance_inputs,
            outputs={
                "lattice_constant": self.a,
                "source": res.get("source"),
                "supercell_repeat": res.get("supercell_repeat"),
                "counts": res.get("counts"),
                "Emin_eV": res.get("Emin_eV"),
                "scan_avals_A": res.get("scan_avals_A"),
                "scan_Evals_eV": res.get("scan_Evals_eV"),
            },
        )

        return self.a

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
    def generate_slab(self, a=None):
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

        write("slab.xyz", self.slab)
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

                calc = name_to_ase_software(self.software)(
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

                calc = name_to_ase_software(self.software)(
                    input_data=input_data,
                    pseudopotentials=self.software_kwargs["pseudopotentials"],
                    kpts=kpts,
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
