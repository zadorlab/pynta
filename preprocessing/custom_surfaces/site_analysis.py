# site_analysis.py
import numpy as np
import yaml
from ase import Atoms
from ase.geometry import get_distances
from ase.geometry.analysis import Analysis

from pynta.utils import get_unique_sym_struct_index_clusters
from pynta.mol import add_adsorbate_to_site
from molecule.molecule import Molecule, Atom, Bond


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
            adsorbate=Atoms("Ne"),
            surf_ind=0,
            site=site,
            height=adsorbate_height
        )
        geoms.append(g)
        sites_per_geom.append([site])

    clusters = get_unique_sym_struct_index_clusters(geoms)

    return (
        [geoms[c[0]] for c in clusters],
        [sites_per_geom[c[0]] for c in clusters]
    )


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


# ============================================================
# 3-fold classification
# ============================================================

def classify_threefold_sites(single_geoms, single_sites_lists):
    admols = []
    geom_indices = []

    for i, (geom, sites) in enumerate(zip(single_geoms, single_sites_lists)):
        for s in sites:
            if s["site"] == "3fold":
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
        label = f"3fold{type_id}"
        for graph_idx in members:
            type_map[geom_indices[graph_idx]] = label

    for geom_idx, sites in enumerate(single_sites_lists):
        for s in sites:
            if s["site"] == "3fold":
                s["site"] = type_map.get(geom_idx, "3fold-unknown")

    return type_map


# ============================================================
# YAML writer (EXACT as requested)
# ============================================================

def write_sites_yaml(single_sites_lists, clusters, filename="sites.yaml"):
    yaml_data = {
        "n_unique_geometries": len(single_sites_lists),
        "n_distinct_three_fold_types": len(clusters),
        "geometries": []
    }

    for i, sites in enumerate(single_sites_lists):
        entry = {
            "geometry_index": i,
            "sites": []
        }
        for s in sites:
            entry["sites"].append({
                "site": s["site"],
                "morphology": s.get("morphology"),
                "indices": s.get("indices"),
                "position": list(map(float, s.get("position", []))),
                "surface": s.get("surface"),
                "label": s.get("label")
            })
        yaml_data["geometries"].append(entry)

    with open(filename, "w") as f:
        yaml.dump(yaml_data, f, sort_keys=False)