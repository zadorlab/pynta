# site_analysis.py
import numpy as np
import json
import copy
from copy import deepcopy

from ase import Atoms
from ase.io import Trajectory, read, write
from ase.geometry import get_distances
from ase.geometry.analysis import Analysis
from ase.neighborlist import NeighborList

from pynta.utils import get_unique_sym_struct_index_clusters
from pynta.mol import add_adsorbate_to_site

from molecule.molecule import Molecule, Atom, Bond
from acat.adsorption_sites import SlabAdsorptionSites


# ============================================================
# Site / geometry utilities
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
        and site1.get("site") == site2.get("site")
        and site1.get("morphology") == site2.get("morphology")
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
            s["bond_length"] = float(best_dist)
            occ.append(s)
    return occ


def generate_unique_sites(slab, sites, nslab, site_bond_cutoff, adsorbate_height, tag_symbol="He"):
    occ = get_occupied_sites(slab, sites, nslab, site_bond_cutoff)
    unocc = [s for s in sites if not any(sites_match(s, o, slab) for o in occ)]

    geoms = []
    sites_per_geom = []

    for site in unocc:
        g = slab.copy()
        add_adsorbate_to_site(
            g,
            adsorbate=Atoms(tag_symbol),
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


# ============================================================
# Graph construction (metals -> X)
# ============================================================

def ase_to_rmg_symbol(sym):
    if sym in {"H","C","N","O","S","F","Cl","Br","I","P","Si","Li"}:
        return sym
    if sym == "He":
        return "Li"
    return "X"


def lone_pairs_for(sym):
    return {"N":1,"P":1,"O":2,"S":2,"F":3,"Cl":3,"Br":3,"I":3}.get(sym, 0)


def generate_graph_slab_fixed(atoms):
    if not all(atoms.pbc[:2]) and _has_reasonable_cell(atoms):
        atoms = atoms.copy()
        atoms.pbc = (True, True, False)
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
            if s.get("site"):
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


def write_sites_json(single_sites_lists, clusters, filename="sites_graph.json"):
    def to_py(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, (np.integer, np.floating)):
            return x.item()
        if isinstance(x, tuple):
            return [to_py(v) for v in x]
        if isinstance(x, list):
            return [to_py(v) for v in x]
        if isinstance(x, dict):
            return {k: to_py(v) for k, v in x.items()}
        if x is None:
            return "null"
        return x

    json_data = {
        "n_unique_geometries": int(len(single_sites_lists)),
        "n_distinct_three_fold_types": int(len(clusters)),
        "geometries": []
    }

    for sites in single_sites_lists:
        for s in sites:
            entry = {
                "site": s.get("site"),
                "surface": s.get("surface", "null"),
                "morphology": s.get("morphology", "null"),
                "position": to_py(s.get("position", [])),
                "normal": to_py(s.get("normal", [])),
                "indices": to_py(s.get("indices", [])),
                "composition": s.get("composition", "null"),
                "subsurf_index": to_py(s.get("subsurf_index", "null")),
                "subsurf_element": s.get("subsurf_element", "null"),
                "label": s.get("label", "null"),
                "topology": to_py(s.get("topology", [])),
            }
            json_data["geometries"].append(to_py(entry))

    with open(filename, "w") as f:
        json.dump(json_data, f, indent=4)


# ============================================================
# Vacancy detect (coordination-based)
# ============================================================

def _pbc_wrap_frac(df):
    return df - np.round(df)

def _has_reasonable_cell(atoms, eps=1e-6):
    cell = atoms.cell.array
    return np.linalg.norm(cell[0]) > eps and np.linalg.norm(cell[1]) > eps

def _top_surface_indices(atoms, nslab, z_tol=None):
    slab = atoms[:nslab]
    z = slab.positions[:, 2]
    zmax = z.max()

    if z_tol is None:
        zr = zmax - z.min()
        z_tol = max(0.6, 0.05 * zr)

    top = np.where(z >= zmax - z_tol)[0]
    if len(top) < 3:
        k = max(3, int(0.10 * nslab))
        top = np.argsort(z)[-k:]
    return np.array(top, dtype=int)

def _estimate_nn_distance_xy(atoms, inds, use_pbc_xy=True):
    pos = atoms.positions[inds, :2]
    N = len(pos)
    if N < 2:
        return 2.5

    if use_pbc_xy and _has_reasonable_cell(atoms):
        cell = atoms.cell.array
        A = cell[:2, :2]
        Ainv = np.linalg.inv(A)

        frac = (pos @ Ainv.T)
        df = frac[:, None, :] - frac[None, :, :]
        df -= np.round(df)
        dxy = df @ A.T
    else:
        dxy = pos[:, None, :] - pos[None, :, :]

    d = np.sqrt((dxy ** 2).sum(axis=2))
    d[d < 1e-8] = np.inf
    nn = np.min(d, axis=1)
    nn_med = np.median(nn[np.isfinite(nn)])
    return float(nn_med if np.isfinite(nn_med) and nn_med > 0 else 2.5)

def _cluster_points(points, thresh):
    M = len(points)
    if M == 0:
        return []
    used = np.zeros(M, dtype=bool)
    clusters = []
    for i in range(M):
        if used[i]:
            continue
        q = [i]
        used[i] = True
        cl = [i]
        while q:
            a = q.pop()
            da = points - points[a]
            dist = np.sqrt((da ** 2).sum(axis=1))
            neigh = np.where((dist <= thresh) & (~used))[0]
            for j in neigh:
                used[j] = True
                q.append(j)
                cl.append(j)
        clusters.append(cl)
    return clusters

def _cluster_points_pbc(points_xy, cell, thresh):
    M = len(points_xy)
    if M == 0:
        return []
    A = np.array(cell)[:2, :2]
    Ainv = np.linalg.inv(A)
    frac = points_xy @ Ainv.T

    used = np.zeros(M, dtype=bool)
    clusters = []
    for i in range(M):
        if used[i]:
            continue
        q = [i]
        used[i] = True
        cl = [i]
        while q:
            a = q.pop()
            df = frac - frac[a]
            df = _pbc_wrap_frac(df)
            dxy = df @ A.T
            dist = np.sqrt((dxy ** 2).sum(axis=1))
            neigh = np.where((dist <= thresh) & (~used))[0]
            for j in neigh:
                used[j] = True
                q.append(j)
                cl.append(j)
        clusters.append(cl)
    return clusters

def detect_vacancy_sites_from_coordination(
    atoms, nslab,
    z_tol=None,
    nn_factor=1.25,
    degree_drop=1,
    cluster_factor=1.2,
    max_sites=20,
):
    top = _top_surface_indices(atoms, nslab, z_tol=z_tol)
    if len(top) < 3:
        return [], atoms

    use_pbc_xy = _has_reasonable_cell(atoms)

    d_nn = _estimate_nn_distance_xy(atoms, top, use_pbc_xy=use_pbc_xy)
    cutoff = nn_factor * d_nn

    cutoffs = [cutoff] * len(atoms)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True, skin=0.0)
    nl.update(atoms)

    top_set = set(top.tolist())

    deg = np.zeros(len(top), dtype=int)
    for ii, a in enumerate(top):
        neigh, _ = nl.get_neighbors(a)
        neigh_top = [j for j in neigh if j in top_set]
        deg[ii] = len(neigh_top)

    expected = int(np.median(deg))

    under_mask = deg <= (expected - degree_drop)
    under = top[under_mask]
    if len(under) == 0:
        return [], atoms

    cluster_thresh = cluster_factor * d_nn
    under_xy = atoms.positions[under, :2]

    if use_pbc_xy and atoms.pbc[0] and atoms.pbc[1]:
        clusters = _cluster_points_pbc(under_xy, atoms.cell.array, thresh=cluster_thresh)
    else:
        clusters = _cluster_points(under_xy, thresh=cluster_thresh)

    z_top = float(np.mean(atoms.positions[top, 2]))

    defect_sites = []
    for k, cl in enumerate(clusters):
        cl = np.array(cl, dtype=int)
        ring_inds = under[cl]
        ring_xy = atoms.positions[ring_inds, :2]
        center_xy = ring_xy.mean(axis=0)

        defect_sites.append({
            "position": np.array([center_xy[0], center_xy[1], z_top], dtype=float),
            "site": "defect",
            "name": f"vacancy_{k}",
        })

    sizes = np.array([len(cl) for cl in clusters], dtype=int)
    order = np.argsort(sizes)[::-1]
    defect_sites = [defect_sites[i] for i in order][:max_sites]

    return defect_sites, atoms


# ============================================================
# Vacancy drop logic (kept from your file)
# ============================================================

def generate_unique_site_additions_vacancy(
    geo, sites, slab, nslab, site_bond_cutoff, xyz_path,
    get_sites,
    site_bond_params_list=None,
    sites_list=None,
    tag_symbol=None,
    dz=0.1,
    max_drop=10.0,
    margin=0.25,
    min_clearance=1.3,
    stable_steps=2,
    noble_gases=("Ne", "Ar", "Kr", "Xe", "Rn"),
    save_all_drop_steps=True,
):
    if site_bond_params_list is None:
        site_bond_params_list = []
    if sites_list is None:
        sites_list = []

    def _estimate_first_neighbor_distance(atoms):
        idx = [i for i, a in enumerate(atoms) if a.symbol not in noble_gases]
        if len(idx) < 2:
            return 3.0
        pos = atoms.positions[idx]
        _, dmat = get_distances(pos, pos, cell=atoms.cell, pbc=atoms.pbc)
        dmat = np.array(dmat)
        np.fill_diagonal(dmat, np.inf)
        nn = np.min(dmat, axis=1)
        nn = nn[np.isfinite(nn)]
        return float(np.median(nn)) if len(nn) else 3.0

    def _build_nl(atoms, cutoff):
        nl = NeighborList([cutoff] * len(atoms), self_interaction=False, bothways=True, skin=0.0)
        nl.update(atoms)
        return nl

    def _tag_neighbors(atoms, i_tag, cutoff):
        nl = _build_nl(atoms, cutoff)
        neigh, _ = nl.get_neighbors(i_tag)
        neigh = [j for j in neigh if atoms[j].symbol not in noble_gases]
        return len(neigh), [int(x) for x in neigh]

    def _min_dist_to_non_noble(atoms, i_tag):
        idx = [j for j, a in enumerate(atoms) if j != i_tag and a.symbol not in noble_gases]
        if not idx:
            return np.inf
        _, d = get_distances(atoms.positions[i_tag:i_tag+1], atoms.positions[idx],
                             cell=atoms.cell, pbc=atoms.pbc)
        return float(np.min(d))

    def _drop_until_max_cn_with_traj(base_atoms, i_tag, cutoff, defect_id):
        atoms = base_atoms.copy()
        frames, meta = [], []

        max_cn = -1
        best_atoms = None
        best_step = None
        no_improve = 0

        nsteps = int(max_drop / dz) + 1
        for step_i in range(nsteps):
            cn, neigh = _tag_neighbors(atoms, i_tag, cutoff)
            mind = _min_dist_to_non_noble(atoms, i_tag)
            z = float(atoms.positions[i_tag, 2])

            if save_all_drop_steps:
                frames.append(atoms.copy())
                meta.append({
                    "defect_id": int(defect_id),
                    "step": int(step_i),
                    "z": float(z),
                    "cn": int(cn),
                    "neighbor_indices": neigh,
                    "min_dist": float(mind),
                    "cutoff": float(cutoff),
                    "dz": float(dz),
                })

            if mind < min_clearance:
                break

            if cn > max_cn:
                max_cn = int(cn)
                best_atoms = atoms.copy()
                best_step = int(step_i)
                no_improve = 0
            else:
                no_improve += 1

            if no_improve >= stable_steps:
                break

            atoms.positions[i_tag, 2] -= dz

        if best_atoms is None:
            best_atoms = atoms.copy()
            best_step = 0
            max_cn, _ = _tag_neighbors(best_atoms, i_tag, cutoff)

        found = {"max_cn": best_atoms, "max_cn_value": int(max_cn), "best_step": int(best_step)}
        return found, frames, meta

    slab_for_detect = read(xyz_path)
    defect_sites, _ = detect_vacancy_sites_from_coordination(slab_for_detect, nslab)

    allowed_tags = ("He", "Ne", "Ar", "Kr", "Xe", "Rn")
    if tag_symbol is None:
        nads = len(geo) - nslab
        sym = allowed_tags[nads % len(allowed_tags)]
    else:
        if tag_symbol not in allowed_tags:
            raise ValueError(f"tag_symbol must be one of {allowed_tags}, got {tag_symbol!r}")
        sym = tag_symbol
    tag = Atoms(sym, positions=[[0, 0, 0]])

    d1 = _estimate_first_neighbor_distance(slab_for_detect)
    cn_cutoff = d1 + margin

    maxcn_geoms = []
    maxcn_meta = []
    drop_geoms = []
    drop_meta = []

    for defect_id, ds in enumerate(defect_sites):
        g0 = geo.copy()
        tag_here = tag.copy()
        tag_here.positions[:] = np.array(ds["position"], dtype=float)
        g0 += tag_here
        i_tag = len(g0) - 1

        found, frames, meta = _drop_until_max_cn_with_traj(
            g0, i_tag=i_tag, cutoff=cn_cutoff, defect_id=defect_id
        )

        if save_all_drop_steps:
            drop_geoms.extend(frames)
            drop_meta.extend(meta)

        best = found["max_cn"]
        maxcn_geoms.append(best)
        maxcn_meta.append({
            "defect_id": int(defect_id),
            "max_cn": int(found["max_cn_value"]),
            "best_step": int(found["best_step"]),
            "site_pos": np.array(ds["position"]).tolist(),
            "cn_cutoff": float(cn_cutoff),
            "tag_symbol": sym,
        })

    occ = get_occupied_sites(geo, sites, nslab, site_bond_cutoff)
    unocc = [site for site in sites if not any(sites_match(site, osite, slab) for osite in occ)]

    geoms = []
    for site in unocc:
        g = geo.copy()
        add_adsorbate_to_site(g, adsorbate=tag, surf_ind=0, site=site, height=1.5)
        geoms.append(g)

    clusters = get_unique_sym_struct_index_clusters(geoms)
    geoms_unique = [geoms[c[0]] for c in clusters]

    geom_all = geoms_unique + maxcn_geoms

    sites_lists_all = []
    for i, atoms in enumerate(geom_all):
        sites_lists_all.append({
            "geom_index": i,
            "sites": get_sites(atoms),
        })

    site_bond_params_all = [deepcopy(site_bond_params_list) for _ in range(len(geom_all))]

    return (
        geom_all,
        site_bond_params_all,
        sites_lists_all,
        maxcn_geoms,
        maxcn_meta,
        drop_geoms,
        drop_meta,
    )


# ============================================================
# Workflows (called by preprocessing.py)
# ============================================================

def workflow_no_defect_unique_sites(
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
    cas = SlabAdsorptionSites(slab, surface_string, composition_effect=True)
    all_sites = cas.get_sites()

    single_geoms, single_sites_lists = generate_unique_sites(
        slab, all_sites, nslab, site_bond_cutoff, adsorbate_height, tag_symbol=tag_symbol
    )

    traj = Trajectory(traj_filename, "w")
    for g in single_geoms:
        traj.write(g)
    traj.close()

    save_sites_to_json(all_sites, filename=sites_json)
    save_neighbor_site_list_to_json(cas, filename=neighbor_json)

    admols, geom_indices = classify_all_sites(single_geoms, single_sites_lists)
    _, clusters = cluster_isomorphic_graphs(admols)
    write_sites_json(single_sites_lists, clusters, filename=sites_graph_json)

    if verbose:
        print(f"There are {len(single_sites_lists)} unique sites out of {len(all_sites)}.")
        print(f"Wrote: {traj_filename}")
        print(f"Wrote: {sites_json}")
        print(f"Wrote: {neighbor_json}")
        print(f"Wrote: {sites_graph_json}")

    return single_geoms, single_sites_lists, clusters


def workflow_defect_vacancy_drop(
    slab, nslab, xyz_path, surface_obj,
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
    geo = slab.copy()
    cas = SlabAdsorptionSites(slab, surface=surface_obj, composition_effect=True)
    sites = cas.get_sites()

    save_neighbor_site_list_to_json(cas, filename=neighbor_json)

    my_get_sites = lambda atoms: SlabAdsorptionSites(
        atoms[:nslab].copy(), surface=surface_obj, composition_effect=True
    ).get_sites()

    (geom_all,
     params_all,
     sites_lists_all,
     maxcn_geoms,
     maxcn_meta,
     drop_geoms,
     drop_meta) = generate_unique_site_additions_vacancy(
        geo=geo,
        sites=sites,
        slab=slab,
        nslab=nslab,
        site_bond_cutoff=site_bond_cutoff,
        xyz_path=xyz_path,
        get_sites=my_get_sites,
        tag_symbol=tag_symbol,
        dz=dz,
        stable_steps=stable_steps,
        max_drop=max_drop,
        margin=margin,
        min_clearance=min_clearance,
        save_all_drop_steps=save_all_drop_steps,
    )

    write(traj_geom_all, geom_all)
    write(traj_drop, drop_geoms)
    write(traj_maxcn, maxcn_geoms)
    save_sites_to_json(sites_lists_all, filename=json_geom_all_sites)

    if verbose:
        print(f"[defect] Number of identified sites: {len(sites_lists_all)}")
        print(f"Wrote: {traj_geom_all}")
        print(f"Wrote: {traj_drop}")
        print(f"Wrote: {traj_maxcn}")
        print(f"Wrote: {json_geom_all_sites}")
        print(f"Wrote: {neighbor_json}")

    return geom_all, sites_lists_all, maxcn_geoms, maxcn_meta, drop_geoms, drop_meta


def workflow_auto(
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
    from acat.settings import CustomSurface

    slab = read(xyz_path)
    nslab = len(slab)
    if not any(slab.pbc) and _has_reasonable_cell(slab):
        slab.pbc = (True, True, False)
    surface_obj = CustomSurface(slab, n_layers=n_layers)

    defect_sites, _ = detect_vacancy_sites_from_coordination(slab, nslab)

    if verbose:
        print(f"Found {len(defect_sites)} defect site(s)")

    if len(defect_sites) == 0:
        return workflow_no_defect_unique_sites(
            slab=slab,
            nslab=nslab,
            adsorbate_height=adsorbate_height,
            site_bond_cutoff=site_bond_cutoff,
            surface_string=surface_string_no_defect,
            tag_symbol="He",
            verbose=verbose,
        )
    else:
        return workflow_defect_vacancy_drop(
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
            verbose=verbose,
        )