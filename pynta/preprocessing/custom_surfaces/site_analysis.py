# site_analysis.py
import numpy as np
import yaml
from ase import Atoms
from ase.geometry import get_distances
from ase.geometry.analysis import Analysis
from ase.neighborlist import NeighborList
import copy
from copy import deepcopy

import json
from ase.io import Trajectory, read, write
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
            height=adsorbate_height,
            tilt_angle=25.24, 
            offset=False)
        geoms.append(g)
        sites_per_geom.append([site])

    clusters = get_unique_sym_struct_index_clusters(geoms)

    return (
        [geoms[c[0]] for c in clusters],
        [sites_per_geom[c[0]] for c in clusters]
    )


def generate_all_sites(slab, sites, nslab, site_bond_cutoff, adsorbate_height):
    """Like generate_unique_sites but returns every unoccupied site without symmetry reduction."""
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
            height=adsorbate_height,
            tilt_angle=25.24,
            offset=False)
        geoms.append(g)
        sites_per_geom.append([site])

    return geoms, sites_per_geom


def suggest_n_layers(slab, dz=0.5, verbose=True):
    """
    Estimate the number of atomic layers in a slab by clustering z-coordinates.

    Atoms within `dz` Angstroms of each other in z are grouped into the same
    layer. The count of distinct groups is the suggested n_layers value for
    ACAT's CustomSurface.

    Parameters
    ----------
    slab    : ASE Atoms
    dz      : z-tolerance (Å) for grouping atoms into the same layer
    verbose : print layer summary

    Returns
    -------
    n_layers : int
    """
    import numpy as np

    z_sorted = np.sort(slab.get_positions()[:, 2])
    layers = []
    current = [z_sorted[0]]
    for zi in z_sorted[1:]:
        if zi - current[-1] < dz:
            current.append(zi)
        else:
            layers.append(current)
            current = [zi]
    layers.append(current)

    n_layers = len(layers)

    if verbose:
        print(f"Suggested n_layers = {n_layers}")
        print()
        print(f"  {'Layer':>6}  {'z mean (Å)':>10}  {'n atoms':>7}")
        print(f"  {'─'*30}")
        for i, layer in enumerate(layers):
            print(f"  {i+1:>6}  {np.mean(layer):>10.3f}  {len(layer):>7}")

    return n_layers


def _reduce_to_representatives(single_geoms, single_sites_lists):
    """Return one geometry and site-list per distinct site label (first occurrence)."""
    seen = {}
    for i, sites in enumerate(single_sites_lists):
        for s in sites:
            label = s.get("site")
            if label and label not in seen:
                seen[label] = i
    rep_indices = sorted(seen.values())
    return (
        [single_geoms[i] for i in rep_indices],
        [single_sites_lists[i] for i in rep_indices],
    )


def write_trajectory_pynta(slab, cas, nslab, site_bond_cutoff, adsorbate_height, trajectory_filename="unique_sites.traj"):
    # Generate unique sites and geometries
    single_geoms, single_sites_lists = generate_unique_sites(
        slab,
        cas.get_sites(),
        nslab,
        site_bond_cutoff,
        adsorbate_height
    )

    # Print the number of unique sites
    print(f'There are {len(single_sites_lists)} unique sites out of {len(cas.get_sites())}.')

    # Create a Trajectory object to write to the specified file
    traj = Trajectory(trajectory_filename, "w")
    all_sites = cas.get_sites()
    # Write each geometry to the trajectory
    for g in single_geoms:
        traj.write(g)
    
    # Close the trajectory file
    traj.close()

    # Save single_sites_lists to JSON
#    save_sites_to_json(single_sites_lists)
    save_sites_to_json(all_sites)

def save_sites_to_json(single_sits_lists, filename='sites.json'):
    """Save data to a JSON file with robust type conversion."""
    import json
    import numpy as np

    def process(x):
        if isinstance(x, dict):
            return {k: process(v) for k, v in x.items()}
        if isinstance(x, list):
            return [process(v) for v in x]
        if isinstance(x, tuple):
            return [process(v) for v in x]
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, (np.integer, np.floating)):
            v = x.item()
            return None if (isinstance(v, float) and (v != v)) else v  # NaN → null
        if x is None:
            return None
        if isinstance(x, float):
            return None if (x != x) else x  # NaN → null
        if isinstance(x, (str, int, bool)):
            return x
        # fallback for objects like CustomSurface, ASE objects, etc.
        return x.__class__.__name__

    with open(filename, "w") as f:
        json.dump(process(single_sits_lists), f, indent=4)

    print(f"Sites data saved to '{filename}'.")

#def save_sites_to_json(single_sites_lists, filename='sites.json'):
#    """Save unique sites to a JSON file."""
#    
#    def process_data(data):
#        """Recursively convert NumPy types to native Python types and replace None with 'null'."""
#        if isinstance(data, dict):
#            return {key: process_data(value) for key, value in data.items()}
#        elif isinstance(data, list):
#            return [process_data(item) for item in data]
#        elif isinstance(data, np.ndarray):
#            return data.tolist()  # Convert NumPy array to list
#        elif isinstance(data, (np.int64, np.float64)):  # Check for NumPy numeric types
#            return data.item()  # Convert to native Python int or float
#        elif data is None:  # Check for None (Python's equivalent of JSON null)
#            return "null"  # Replace with the string "null"
#        else:
#            return data  # Return the data as is if it's already a native type
#
#    # Flatten the input if it's a list of lists
#    if isinstance(single_sites_lists, list) and all(isinstance(item, list) for item in single_sites_lists):
#        # Flatten the list of lists into a single list of dictionaries
#        processed_sites_data = process_data([item for sublist in single_sites_lists for item in sublist])
#    else:
#        # If it's already a single list of dictionaries, process it directly
#        processed_sites_data = process_data(single_sites_lists)
#
#    # Save to JSON file
#    with open(filename, "w") as f:
#        json.dump(processed_sites_data, f, indent=4)
#
#    print(f"Sites data saved to '{filename}' with None replaced by 'null'.")


def write_trajectory_for_acat(slab, cas, trajectory_filename):
    # Create a Trajectory object to write to the specified file
    traj = Trajectory(trajectory_filename, 'w')
    
    # Iterate over all unique sites
    for site in cas.get_unique_sites():
        # Create a deep copy of the original slab
        my_slab = copy.deepcopy(slab)
        
        # Create a new Atom at the unique site's position
        new_position = Atom('He', site["position"])
        
        # Append the new Atom to the slab
        my_slab.append(new_position)
        
        # Write the current configuration to the trajectory
        traj.write(my_slab)
    
    # Close the trajectory file
    traj.close()



def save_neighbor_site_list_to_json(cas, filename='neighbor_site_list.json'):
    """Retrieve the neighbor site list from cas, convert NumPy types to native Python types, and save it to a JSON file."""
    
    def convert_numpy_types(data):
        """Recursively convert NumPy types to native Python types."""
        if isinstance(data, dict):
            return {key: convert_numpy_types(value) for key, value in data.items()}
        elif isinstance(data, list):
            return [convert_numpy_types(item) for item in data]
        elif isinstance(data, np.ndarray):
            return data.tolist()  # Convert NumPy array to list
        elif isinstance(data, (np.int64, np.float64)):  # Check for NumPy numeric types
            return data.item()  # Convert to native Python int or float
        else:
            return data  # Return the data as is if it's already a native type

    # Get the neighbor site list
    neighbor_site_list = cas.get_neighbor_site_list()

    # Convert NumPy types to native Python types
    neighbor_site_list = convert_numpy_types(neighbor_site_list)

    # Save the neighbor site list to a JSON file
    with open(filename, 'w') as f:
        json.dump(neighbor_site_list, f, indent=4)  # Use indent for pretty printing

    print(f"Neighbor site list saved to '{filename}'.")


# 
# ===================================================
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
# site classification
# ============================================================

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


def _classify_defect_sites_by_cn_and_graph(cn_geoms, cn_site_dicts):
    """
    Deduplicate CN>0 void-probe geometries for labeled_sites.json:

      1. Group by coordination number (cn key in each site dict).
      2. Within each CN group, cluster by RMG graph isomorphism.
      3. In each isomorphic cluster keep the geometry whose probe is
         closest to the median probe position of the cluster.

    Returns
    -------
    rep_geoms : list[Atoms]
    rep_sites : list[dict]   — site dicts with site="defect", morphology="defect"
    """
    if not cn_geoms:
        return [], []

    from collections import defaultdict

    _NOBLE = {"He", "Ne", "Ar", "Kr", "Xe", "Rn"}

    def _probe_pos(atoms):
        for i, a in enumerate(atoms):
            if a.symbol in _NOBLE:
                return atoms.positions[i].copy()
        return np.zeros(3)

    # group frame indices by CN
    cn_groups = defaultdict(list)
    for k, s in enumerate(cn_site_dicts):
        cn_groups[s.get("cn", 0)].append(k)

    rep_geoms, rep_sites = [], []

    for cn_val in sorted(cn_groups.keys()):
        members = cn_groups[cn_val]

        if len(members) == 1:
            rep_geoms.append(cn_geoms[members[0]])
            rep_sites.append(cn_site_dicts[members[0]])
            continue

        # build RMG graph for each frame in this CN group
        admols = [generate_graph_slab_fixed(cn_geoms[i]) for i in members]

        # union-find isomorphism clustering
        parent = list(range(len(admols)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i in range(len(admols)):
            for j in range(i + 1, len(admols)):
                if admols[i].is_isomorphic(admols[j], strict=True):
                    ri, rj = find(i), find(j)
                    if ri != rj:
                        parent[rj] = ri

        iso_clusters = defaultdict(list)
        for i in range(len(admols)):
            iso_clusters[find(i)].append(i)

        for cluster_local in iso_clusters.values():
            orig = [members[li] for li in cluster_local]
            probe_positions = np.array([_probe_pos(cn_geoms[oi]) for oi in orig])
            median_pos = np.median(probe_positions, axis=0)
            dists = np.linalg.norm(probe_positions - median_pos, axis=1)
            best = orig[int(np.argmin(dists))]

            g = cn_geoms[best].copy()
            s = dict(cn_site_dicts[best])
            s["position"] = probe_positions[int(np.argmin(dists))].tolist()
            rep_geoms.append(g)
            rep_sites.append(s)

    return rep_geoms, rep_sites


def update_site_labels_by_graph_and_type(single_sites_lists, clusters, geom_indices):
    """
    Assign labels of the form {site_type}{N} (e.g. "3fold0", "3fold1", "bridge0").

    Two sites receive the same label iff:
      1. Their local slab graphs are isomorphic (same cluster).
      2. They share the same "site" value.
      3. They share the same "morphology" value.

    N counts distinct (cluster, site, morphology) triplets per site_type,
    starting from 0.  Morphology is kept as a differentiator so that
    "3fold_terrace" and "3fold_sc-tc" always produce distinct labels even
    when they land in the same isomorphic cluster — but the morphology name
    does NOT appear in the label string itself.

    Returns
    -------
    geom_to_label : dict  {geom_idx: label_string}
    key_to_label  : dict  {(cluster_id, site, morphology): label_string}
    """
    def _first_site_morph(geom_idx):
        for s in single_sites_lists[geom_idx]:
            if s.get("site"):
                return s.get("site"), s.get("morphology")
        return None, None

    type_counters = {}   # site_val -> next available integer
    key_to_label  = {}   # (cluster_id, site_val, morph_val) -> label
    geom_to_label = {}

    for cluster_id, members in enumerate(clusters.values()):
        for graph_idx in members:
            geom_idx = geom_indices[graph_idx]
            site_val, morph_val = _first_site_morph(geom_idx)
            key = (cluster_id, site_val, morph_val)
            if key not in key_to_label:
                n = type_counters.get(site_val, 0)
                type_counters[site_val] = n + 1
                prefix = site_val if site_val else "site"
                key_to_label[key] = f"{prefix}{n}"
            geom_to_label[geom_idx] = key_to_label[key]

    for geom_idx, sites in enumerate(single_sites_lists):
        label = geom_to_label.get(geom_idx)
        if label is not None:
            for s in sites:
                if s.get("site"):
                    s["site"] = label

    return geom_to_label, key_to_label


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

#def update_site_labels(single_sites_lists, clusters):
#    # Create a mapping from cluster index to its representative site label
#    site_label_mapping = {}
#    
#    # Debugging: Print the structure of clusters
#    print("Clusters structure:", clusters)
#    print("Single sites lists structure:", single_sites_lists)
#
#    for cluster_index, members in clusters.items():
#        # Ensure members is a list of integers
#        if not isinstance(members, list):
#            print(f"Warning: Expected a list for members in cluster {cluster_index}, got {type(members)}")
#            continue
#        
#        # Debugging: Print the first member being accessed
#        print(f"Accessing member index: {members[0]} in single_sites_lists")
#        
#        # Get the site label from the first member of the cluster
#        representative_site = single_sites_lists[members[0]][0]["site"]  # Access the first element of the inner list
#        
#        # Store the representative site label in the mapping
#        site_label_mapping[cluster_index] = representative_site
#        
#        # Update the site label for all members in the cluster
#        for member in members:
#            single_sites_lists[member][0]["site"] = representative_site  # Access the first element of the inner list
#
#    return single_sites_lists, site_label_mapping

def update_site_labels(single_sites_lists, clusters, trajectory_filename="unique_site_graph.traj"):
    # Create a mapping from cluster index to its representative site label
    site_label_mapping = {}
    
    # New list to hold the updated single_sites_lists with only the first member of each cluster
    updated_sites_lists = []
    
    # Debugging: Print the structure of clusters
    print("Clusters structure:", clusters)
    print("Single sites lists structure:", single_sites_lists)

    for cluster_index, members in clusters.items():
        # Ensure members is a list of integers
        if not isinstance(members, list):
            print(f"Warning: Expected a list for members in cluster {cluster_index}, got {type(members)}")
            continue
        
        # Get the first member of the cluster
        first_member_index = members[0]
        
        # Get the site label from the first member of the cluster
        representative_site = single_sites_lists[first_member_index][0]["site"]
        
        # Store the representative site label in the mapping
        site_label_mapping[cluster_index] = representative_site
        
        # Add the first member to the updated list
        updated_sites_lists.append(single_sites_lists[first_member_index])
        
        # Update the site label for all members in the cluster
        for member in members:
            single_sites_lists[member][0]["site"] = representative_site  # Access the first element of the inner list

    return updated_sites_lists, site_label_mapping

# ============================================================
# json writer (EXACT as requested)
# ============================================================


def write_sites_json(single_sites_lists, clusters, filename="sites_graph.json"):
    def to_py(x):
        import numpy as np
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
                "normal": to_py(s.get("normal") or [0.0, 0.0, 1.0]),
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

def write_trajectory_graph(updated_sites_lists, clusters, trajectory_filename="unique_sites_graph.traj"):
    # Print the number of unique sites
    print(f'There are {len(updated_sites_lists)} unique sites out of {len(clusters)}.')

    # Create a Trajectory object to write to the specified file
    traj = Trajectory(trajectory_filename, "w")

    # Write each geometry to the trajectory
    for site in updated_sites_lists:
        # Assuming each site is a list containing a dictionary with atomic information
        # Extract the relevant data to create an Atoms object
        # Example: site[0] should contain the necessary data to create an Atoms object
        # You may need to adjust this based on the actual structure of your data
        atoms_data = site[0]  # Extract the first dictionary from the inner list
        # Create an Atoms object (you need to adapt this based on your data structure)
        atoms = Atoms(
            symbols=atoms_data['composition'],  # Assuming 'composition' contains the symbols
            positions=atoms_data['position'],    # Assuming 'position' contains the atomic positions
            # Add other necessary parameters like cell, pbc, etc. if needed
        )
        traj.write(atoms)  # Write the Atoms object to the trajectory

    # Close the trajectory file
    traj.close()

#============ vacancy detect ================

def _pbc_wrap_frac(df):
    """wrap fractional delta to [-0.5, 0.5)"""
    return df - np.round(df)

def _cluster_points_pbc(points_xy, cell, thresh):
    """
    PBC-aware clustering in xy using minimum-image distances.
    points_xy: (M,2) cartesian
    cell: ASE cell array (3,3) or (2,2) usable for xy
    """
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
    """
    Reference-free vacancy detection on the top surface using in-plane
    missing coordination (graph degree loss).

    Returns: (defect_sites, atoms)
      defect_sites: list of dicts with keys: position, site, name
    """
    # 1) top surface atoms
    top = _top_surface_indices(atoms, nslab, z_tol=z_tol)
    if len(top) < 3:
        return [], atoms

    use_pbc_xy = _has_reasonable_cell(atoms)

    # 2) estimate in-plane NN distance (sets cutoff scale)
    d_nn = _estimate_nn_distance_xy(atoms, top, use_pbc_xy=use_pbc_xy)
    cutoff = nn_factor * d_nn

    # 3) neighborlist (use same cutoff for all atoms)
    cutoffs = [cutoff] * len(atoms)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True, skin=0.0)
    nl.update(atoms)

    top_set = set(top.tolist())

    # 4) build in-plane adjacency and degrees among top atoms
    deg = np.zeros(len(top), dtype=int)
    for ii, a in enumerate(top):
        neigh, offsets = nl.get_neighbors(a)
        # keep only top-layer neighbors
        neigh_top = [j for j in neigh if j in top_set]
        deg[ii] = len(neigh_top)

    # robust "expected" coordination on this surface (no reference)
    expected = int(np.median(deg))

    # 5) undercoordinated atoms: degree at least `degree_drop` less than typical
    under_mask = deg <= (expected - degree_drop)
    under = top[under_mask]
    if len(under) == 0:
        return [], atoms

    # 6) cluster undercoordinated atoms (these tend to form a ring around vacancy)
    cluster_thresh = cluster_factor * d_nn
    under_xy = atoms.positions[under, :2]

    if use_pbc_xy and atoms.pbc[0] and atoms.pbc[1]:
        clusters = _cluster_points_pbc(under_xy, atoms.cell.array, thresh=cluster_thresh)
    else:
        clusters = _cluster_points(under_xy, thresh=cluster_thresh)

    # 7) place vacancy site for each cluster (centroid of ring in xy, z at top)
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

    # 8) score sites: bigger vacancy tends to create larger "undercoordination region".
    # Here a simple score = number of undercoordinated atoms in the cluster.
    # You can replace this with a geometric score if you want.
    sizes = np.array([len(cl) for cl in clusters], dtype=int)
    order = np.argsort(sizes)[::-1]
    defect_sites = [defect_sites[i] for i in order][:max_sites]

    return defect_sites, atoms

def _has_reasonable_cell(atoms, eps=1e-6):
    cell = atoms.cell.array
    return np.linalg.norm(cell[0]) > eps and np.linalg.norm(cell[1]) > eps

def _top_surface_indices(atoms, nslab, z_tol=None):
    slab = atoms[:nslab]
    z = slab.positions[:, 2]
    zmax = z.max()

    # auto z_tol: use small fraction of z-range if not provided
    if z_tol is None:
        zr = zmax - z.min()
        z_tol = max(0.6, 0.05 * zr)  # Å, decent default for most slabs

    top = np.where(z >= zmax - z_tol)[0]
    if len(top) < 3:
        # fallback: take top 10% atoms by z
        k = max(3, int(0.10 * nslab))
        top = np.argsort(z)[-k:]
    return np.array(top, dtype=int)

def _estimate_nn_distance_xy(atoms, inds, use_pbc_xy=True):
    """
    Estimate typical nearest-neighbor distance in XY among atoms[inds].
    Uses O(N^2) distances, OK for typical slabs.
    """
    pos = atoms.positions[inds, :2]
    N = len(pos)
    if N < 2:
        return 2.5

    if use_pbc_xy and _has_reasonable_cell(atoms):
        # fractional in xy only
        cell = atoms.cell.array
        A = cell[:2, :2]  # 2x2
        Ainv = np.linalg.inv(A)

        frac = (pos @ Ainv.T)  # N x 2
        # pairwise fractional diffs with minimum image
        df = frac[:, None, :] - frac[None, :, :]
        df -= np.round(df)  # wrap to [-0.5,0.5]
        dxy = df @ A.T
    else:
        dxy = pos[:, None, :] - pos[None, :, :]

    d = np.sqrt((dxy ** 2).sum(axis=2))
    d[d < 1e-8] = np.inf
    nn = np.min(d, axis=1)
    nn_med = np.median(nn[np.isfinite(nn)])
    return float(nn_med if np.isfinite(nn_med) and nn_med > 0 else 2.5)

def _distance_field_on_cell(atoms, top_inds, grid_n=80):
    """
    Compute distance-to-nearest-top-atom on an XY grid covering one unit cell.
    Uses PBC in xy if the cell is present; otherwise uses bounding box.
    Returns:
      X, Y (grid coords in cartesian), D (distance field), z_top (avg top z)
    """
    slab = atoms
    top_pos = slab.positions[top_inds]
    z_top = float(np.mean(top_pos[:, 2]))

    use_pbc_xy = _has_reasonable_cell(slab)
    if use_pbc_xy:
        cell = slab.cell.array
        A = cell[:2, :2]      # 2x2
        Ainv = np.linalg.inv(A)

        # grid in fractional coordinates [0,1)
        u = np.linspace(0, 1, grid_n, endpoint=False)
        v = np.linspace(0, 1, grid_n, endpoint=False)
        UU, VV = np.meshgrid(u, v, indexing="ij")
        frac_grid = np.stack([UU.ravel(), VV.ravel()], axis=1)  # Mx2
        xy_grid = frac_grid @ A  # Mx2 cartesian xy

        # top atoms in fractional xy
        top_xy = top_pos[:, :2]
        top_frac = top_xy @ Ainv.T  # Nx2

        # compute min distance with minimum image in fractional space
        df = frac_grid[:, None, :] - top_frac[None, :, :]
        df -= np.round(df)
        dxy = df @ A.T
        D = np.sqrt((dxy ** 2).sum(axis=2)).min(axis=1)

        X = xy_grid[:, 0].reshape(grid_n, grid_n)
        Y = xy_grid[:, 1].reshape(grid_n, grid_n)
        D = D.reshape(grid_n, grid_n)
        return X, Y, D, z_top, True
    else:
        # no cell: use bounding box of top atoms
        top_xy = top_pos[:, :2]
        xmin, ymin = top_xy.min(axis=0)
        xmax, ymax = top_xy.max(axis=0)

        u = np.linspace(xmin, xmax, grid_n)
        v = np.linspace(ymin, ymax, grid_n)
        UU, VV = np.meshgrid(u, v, indexing="ij")
        grid_xy = np.stack([UU.ravel(), VV.ravel()], axis=1)

        dxy = grid_xy[:, None, :] - top_xy[None, :, :]
        D = np.sqrt((dxy ** 2).sum(axis=2)).min(axis=1)

        X = UU
        Y = VV
        D = D.reshape(grid_n, grid_n)
        return X, Y, D, z_top, False

def _local_maxima_mask(D):
    """
    Simple 8-neighbor local maxima mask for 2D array D.
    """
    # pad with -inf so edges can still be maxima
    P = np.pad(D, 1, mode="constant", constant_values=-np.inf)
    center = P[1:-1, 1:-1]
    nbrs = [
        P[:-2, :-2], P[:-2, 1:-1], P[:-2, 2:],
        P[1:-1, :-2],              P[1:-1, 2:],
        P[2:, :-2],  P[2:, 1:-1],  P[2:, 2:]
    ]
    is_max = np.ones_like(center, dtype=bool)
    for n in nbrs:
        is_max &= center >= n
    return is_max

def _cluster_points(points, thresh):
    """
    BFS clustering of 2D points: cluster if within `thresh`.
    points: (M,2)
    returns list of clusters (list of point indices)
    """
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

def detect_defect_sites_from_xyz(xyz_path, nslab, grid_n=90,
                                 z_tol=None,
                                 hole_thresh_factor=0.55,
                                 cluster_factor=0.60,
                                 max_sites=20):
    """
    Detect 'hole centers' (vacancy/defect-like empty pockets) on top surface.
    Returns list of site dicts: {"position": np.array([x,y,z]), "site":"defect", "name":...}
    """
    slab = read(xyz_path)
    #slab = xyz_path
    top = _top_surface_indices(slab, nslab, z_tol=z_tol)

    use_pbc_xy = _has_reasonable_cell(slab)
    d_nn = _estimate_nn_distance_xy(slab, top, use_pbc_xy=use_pbc_xy)

    X, Y, D, z_top, _ = _distance_field_on_cell(slab, top, grid_n=grid_n)
    is_max = _local_maxima_mask(D)

    # "hole" if far enough from nearest top atom
    hole_thresh = hole_thresh_factor * d_nn
    cand = np.argwhere(is_max & (D >= hole_thresh))
    if len(cand) == 0:
        return [], slab

    # candidate xy points
    cand_xy = np.column_stack([X[cand[:, 0], cand[:, 1]], Y[cand[:, 0], cand[:, 1]]])

    # cluster them so one defect yields one site
    cluster_thresh = cluster_factor * d_nn
    clusters = _cluster_points(cand_xy, thresh=cluster_thresh)

    # pick cluster center at the maximum-D point within each cluster
    defect_sites = []
    for k, cl in enumerate(clusters):
        cl_inds = np.array(cl, dtype=int)
        # choose point with largest D
        # map back to D values by nearest cand index
        # (cand_xy indices align with `cand` order)
        best_local = cl_inds[np.argmax(D[cand[cl_inds, 0], cand[cl_inds, 1]])]
        x, y = cand_xy[best_local]
        defect_sites.append({
            "position": np.array([x, y, z_top], dtype=float),
            "site": "defect",
            "name": f"defect_{k}",
        })

    # sort by "depth" proxy: larger distance means bigger hole
    # (more likely to be a vacancy pocket)
    # recompute score using nearest top atom distance for each site
# --- NEW: keep only "real defects" (largest pockets) ---
    top_xy = slab.positions[top, :2]

    def pocket_score(site):
        # distance from site xy to nearest top atom (bigger = deeper pocket)
        d = np.sqrt(((top_xy - site["position"][:2]) ** 2).sum(axis=1)).min()
        return float(d)

    scores = np.array([pocket_score(s) for s in defect_sites], dtype=float)
    max_score = scores.max()

    # Keep sites that are close to the maximum pocket size.
    # For a single vacancy, this usually leaves exactly 1 site.
    keep_frac = 0.90  # try 0.90–0.95 for close-packed metals
    keep = scores >= keep_frac * max_score

    defect_sites = [s for s, k in zip(defect_sites, keep) if k]

    # Optional: if still many (e.g., multiple vacancies), cap by max_sites
    # and keep the largest ones.
    scores = np.array([pocket_score(s) for s in defect_sites], dtype=float)
    order = np.argsort(scores)[::-1]
    defect_sites = [defect_sites[i] for i in order][:max_sites]
    # --- end NEW ---

    return defect_sites, slab

def get_adsorbate_dist_from_center(atoms,nslab):
    cell_center = sum(atoms.cell)/3.5
    adatoms = atoms[nslab:]
    adcenter = sum([a.position for a in adatoms])/len(adatoms)
    return np.linalg.norm(adcenter - cell_center)

#===fix=== ``
#def generate_unique_site_additions_vacancy(geo, sites, slab, nslab, site_bond_cutoff, xyz_path,
#                                          site_bond_params_list=None,
#                                          sites_list=None):
#    if site_bond_params_list is None:
#        site_bond_params_list = []
#    if sites_list is None:
#        sites_list = []
#
#    # --- NEW: detect defect/hole sites from the xyz and add to sites ---
#    #defect_sites, _ = detect_defect_sites_from_xyz(xyz_path, nslab)
#    slab_for_detect = read(xyz_path)
#    defect_sites, _ = detect_vacancy_sites_from_coordination(slab_for_detect, nslab)
#
#    # add defect sites if not already present (uses your existing sites_match)
#    for ds in defect_sites:
#        ds = dict(ds)  # make sure it's a mutable copy
#
#        # REQUIRED for add_adsorbate_to_site()
#        ds.setdefault("normal", np.array([0.0, 0.0, 1.0]))
#
#        # REQUIRED for sites_match()
#        ds.setdefault("morphology", "defect")
#
#        # Optional: keep surface label consistent (if your other sites have it)
#        if "surface" not in ds and len(sites) > 0 and "surface" in sites[0]:
#            ds["surface"] = sites[0]["surface"]
#
#        if not any(sites_match(ds, s, slab) for s in sites):
#            sites = sites + [ds]
#    # --- end NEW ---
#
#    nads = len(geo) - nslab
#
#    # label sites with unique noble gas atoms
#    he = Atoms('He', positions=[[0, 0, 0]])
#    ne = Atoms('Ne', positions=[[0, 0, 0]])
#    ar = Atoms('Ar', positions=[[0, 0, 0]])
#    kr = Atoms('Kr', positions=[[0, 0, 0]])
#    xe = Atoms('Xe', positions=[[0, 0, 0]])
#    rn = Atoms('Rn', positions=[[0, 0, 0]])
#    site_tags = [he, ne, ar, kr, xe, rn]
#    tag = site_tags[nads]
#
#    occ = get_occupied_sites(geo, sites, nslab, site_bond_cutoff)
#    unocc = [site for site in sites if not any(sites_match(site, osite, slab) for osite in occ)]
#
#    site_bond_params_lists = [deepcopy(site_bond_params_list) for _ in range(len(unocc))]
#    sites_lists = [deepcopy(sites_list) for _ in range(len(unocc))]
#
#    geoms = []
#    for i, site in enumerate(unocc):
#        geom = geo.copy()
#        add_adsorbate_to_site(geom, adsorbate=tag, surf_ind=0, site=site, height=1.5)
#
#        pos = np.array(site["position"]).tolist()
#        params = {"site_pos": pos, "ind": None, "k": 100.0, "deq": 0.0}
#        site_bond_params_lists[i].append(params)
#        sites_lists[i].append(site)
#        geoms.append(geom)
#
#    indclusters = get_unique_sym_struct_index_clusters(geoms)
#
#    inds = []
#    for cluster in indclusters:        
#        min_dist = np.inf
#        indout = None
#        for ind in cluster:
#            d = get_adsorbate_dist_from_center(geoms[ind], nslab)
#            if d < min_dist:
#                indout = ind
#                min_dist = d
#        inds.append(indout)
#
#    return ([geoms[ind] for ind in inds],
#            [site_bond_params_lists[ind] for ind in inds],
#            [sites_lists[ind] for ind in inds])

def generate_unique_site_additions_vacancy(
    geo, sites, slab, nslab, site_bond_cutoff, xyz_path,
    get_sites,                          # <-- REQUIRED: function that returns site dicts for an Atoms
    site_bond_params_list=None,
    sites_list=None,

    # --- Noble gas selection ---
    tag_symbol=None,

    # --- Drop controls ---
    dz=0.1,
    max_drop=10.0,
    margin=0.25,
    min_clearance=1.3,
    stable_steps=2,                     # "patience": stop after this many non-improving cn steps
    noble_gases=("Ne", "Ar", "Kr", "Xe", "Rn"),

    # --- Output controls ---
    save_all_drop_steps=True,
):
    """
    Returns:
      geom_all,                 # geoms_unique + maxcn_geoms
      site_bond_params_all,     # placeholder list aligned with geom_all (empty params by default)
      sites_lists_all,          # sites_lists_all[i] == get_sites(geom_all[i])
      maxcn_geoms,
      maxcn_meta,
      drop_geoms,
      drop_meta
    """
    if site_bond_params_list is None:
        site_bond_params_list = []
    if sites_list is None:
        sites_list = []

    import numpy as np
    from copy import deepcopy
    from ase.io import read
    from ase import Atoms
    from ase.neighborlist import NeighborList
    from ase.geometry import get_distances

    # ---------- helpers ----------
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
        for step in range(nsteps):
            cn, neigh = _tag_neighbors(atoms, i_tag, cutoff)
            mind = _min_dist_to_non_noble(atoms, i_tag)
            z = float(atoms.positions[i_tag, 2])

            if save_all_drop_steps:
                frames.append(atoms.copy())
                meta.append({
                    "defect_id": int(defect_id),
                    "step": int(step),
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
                best_step = int(step)
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

    # ---------- detect vacancies (ONLY for locating where to drop) ----------
    slab_for_detect = read(xyz_path)
    defect_sites, _ = detect_vacancy_sites_from_coordination(slab_for_detect, nslab)

    # ---------- choose noble gas tag ----------
    allowed_tags = ("He", "Ne", "Ar", "Kr", "Xe", "Rn")
    if tag_symbol is None:
        nads = len(geo) - nslab
        sym = allowed_tags[nads % len(allowed_tags)]
    else:
        if tag_symbol not in allowed_tags:
            raise ValueError(f"tag_symbol must be one of {allowed_tags}, got {tag_symbol!r}")
        sym = tag_symbol
    tag = Atoms(sym, positions=[[0, 0, 0]])

    # ---------- cutoff for cn counting ----------
    d1 = _estimate_first_neighbor_distance(slab_for_detect)
    cn_cutoff = d1 + margin

    # ---------- build maxcn_geoms by dropping tag at each detected vacancy ----------
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

    # ---- ADD THIS HERE (after the loop) ----
    # Keep only ONE defect geometry: the one with the highest max CN
    if len(maxcn_meta) > 1:
        best_j = int(np.argmax([m["max_cn"] for m in maxcn_meta]))
        maxcn_geoms = [maxcn_geoms[best_j]]
        maxcn_meta  = [maxcn_meta[best_j]]
    # ---- end add ----

    # ---------- "main" geoms: keep your original behavior for normal sites ----------
    # Use original `sites` list (terrace sites) to generate site-tagged structures.
    # IMPORTANT: we do NOT add ds as defect sites to `sites` anymore.
    occ = get_occupied_sites(geo, sites, nslab, site_bond_cutoff)
    unocc = [site for site in sites if not any(sites_match(site, osite, slab) for osite in occ)]

    geoms = []
    for site in unocc:
        g = geo.copy()
        add_adsorbate_to_site(g, adsorbate=tag, surf_ind=0, site=site, height=1.5)
        geoms.append(g)

    # all terrace sites, no symmetry reduction
    geoms_unique = geoms
    unocc_rep_sites = unocc

    # ---------- NEW: geom_all = geoms_unique + maxcn_geoms ----------
    geom_all = geoms_unique + maxcn_geoms

    # Build per-geometry "which site is this geometry?" mapping
    sites_lists_all = []
    
    # For geoms_unique: one entry per symmetry-unique terrace site
    for i, s in enumerate(unocc_rep_sites):
        s2 = dict(s)
        if "topology" not in s2:
            s2["topology"] = s2.get("indices", [])
        if "surface" in s2 and s2["surface"] is not None and not isinstance(s2["surface"], str):
            s2["surface"] = s2["surface"].__class__.__name__
    
        sites_lists_all.append({
            "geom_index": i,
            "sites": [s2],
        })
    
    # For maxcn_geoms: append ONCE per defect geometry (usually 1 after your "keep best" filter)
    offset = len(geoms_unique)
    for j, meta in enumerate(maxcn_meta):
        s_def = {
            "site": "defect",
            "surface": "CustomSurface",
            "morphology": "defect",
            "position": meta["site_pos"],
            "normal": [0.0, 0.0, 1.0],
            "indices": [],
            "composition": "null",
            "subsurf_index": "null",
            "subsurf_element": "null",
            "label": meta.get("tag_symbol", "null"),
            "topology": [],
        }
        sites_lists_all.append({
            "geom_index": offset + j,
            "sites": [s_def],
        })

#----old---
    # reduce by symmetry like before    
    #clusters = get_unique_sym_struct_index_clusters(geoms)
    #geoms_unique = [geoms[c[0]] for c in clusters]

    # ---------- NEW: geom_all = geoms_unique + maxcn_geoms ----------
    #geom_all = geoms_unique + maxcn_geoms

    # ---------- NEW: sites_lists are computed from get_sites(geom_all[i]) ----------
    # This ensures morphology/site labels come from get_sites() (not "defect").
    #sites_lists_all = []
    #for i, atoms in enumerate(geom_all):
    #    sites_lists_all.append({
    #        "geom_index": i,
    #        "sites": get_sites(atoms),
    #    })
    # Sites are slab-defined and identical for every geometry in geom_all, so store once.
    #sites_lists_all = [{
    #    "geom_index": 0,
    #    "sites": get_sites(geo),   # clean slab
    #}]
    # If your downstream expects params aligned with images, create placeholders
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

NOBLE_GASES = {"Ne", "Ar", "Kr", "Xe", "Rn"}

def tag_index(atoms, noble_gases=NOBLE_GASES):
    """Return the index of the single noble-gas tag atom in this Atoms."""
    inds = [i for i, a in enumerate(atoms) if a.symbol in noble_gases]
    if len(inds) != 1:
        raise ValueError(f"Expected exactly 1 noble gas tag atom, found {len(inds)}")
    return inds[0]

def unique_drop_frames_by_tag_position(frames, tol=0.25, noble_gases=NOBLE_GASES):
    """
    Keep only frames that have unique tag positions (PBC-aware), within tol (Å).
    frames: list[ase.Atoms] for ONE defect_id (same cell/pbc).
    """
    if not frames:
        return []

    uniq = []
    for at in frames:
        it = tag_index(at, noble_gases=noble_gases)
        keep = True
        for u in uniq:
            iu = tag_index(u, noble_gases=noble_gases)
            _, d = get_distances([at.positions[it]], [u.positions[iu]],
                                 cell=at.cell, pbc=at.pbc)
            if float(d) < tol:
                keep = False
                break
        if keep:
            uniq.append(at)
    return uniq

def write_unique_sites_with_drop_traj(drop_geoms, drop_meta, out_traj="vacancy_unique_sites_with_drop.traj",
                                     tol=0.25, noble_gases=NOBLE_GASES):
    """
    Build a trajectory of UNIQUE sites encountered during noble-gas dropping into vacancy,
    including inside the defect. Uniqueness is based on tag position under PBC.

    drop_geoms, drop_meta come from generate_unique_site_additions_vacancy().
    """
    defect_ids = sorted({m["defect_id"] for m in drop_meta})
    all_unique = []

    for did in defect_ids:
        inds = [i for i, m in enumerate(drop_meta) if m["defect_id"] == did]
        frames = [drop_geoms[i] for i in inds]

        uniq_frames = unique_drop_frames_by_tag_position(frames, tol=tol, noble_gases=noble_gases)

        for k, at in enumerate(uniq_frames):
            at2 = at.copy()
            at2.info["defect_id"] = int(did)
            at2.info["drop_unique_id"] = int(k)
            all_unique.append(at2)

    write(out_traj, all_unique)
    return all_unique

# ============================================================
# High-level workflows (keeps notebooks clean)
# ============================================================

def workflow_detect_vacancies(slab, nslab, verbose=True):
    """Detect vacancy sites and optionally print a summary."""
    defect_sites, _ = detect_vacancy_sites_from_coordination(slab, nslab)
    if verbose:
        print(f"Found {len(defect_sites)} defect site(s)")
        for s in defect_sites:
            x, y, z = s["position"]
            print(f"{s.get('name','defect')}: x={x:.3f}  y={y:.3f}  z={z:.3f}  site={s.get('site')}")
    return defect_sites


def workflow_no_defect_unique_sites(
    slab,
    nslab,
    adsorbate_height=1.0,
    site_bond_cutoff=1.5,
    surface="fcc111",          # str (e.g. "fcc111") or CustomSurface instance
    traj_filename="unique_sites.traj",
    sites_json="sites.json",
    labeled_sites_json="labeled_sites.json",
    neighbor_json="neighbor_site_list.json",
    verbose=True,
):
    """
    No-defect workflow.

    `surface` can be:
      - a conventional ACAT surface string ("fcc111", "fcc332", ...)
      - an ACAT CustomSurface instance

    In both cases the pipeline is:
      1. Get all ACAT sites → sites.json
      2. Generate one tagged geometry per unoccupied site
      3. Build RMG graphs, cluster by isomorphism + (site, morphology)
      4. Assign labels: 3fold0, 3fold1, bridge0, ...
      5. Reduce to one representative per distinct label → labeled_sites.json + traj
    """
    from ase.io.trajectory import Trajectory
    from acat.settings import CustomSurface

    use_custom = isinstance(surface, CustomSurface)
    if use_custom:
        cas = SlabAdsorptionSites(slab, surface=surface, composition_effect=True)
    else:
        cas = SlabAdsorptionSites(slab, surface, composition_effect=True)

    all_sites = cas.get_sites()
    save_sites_to_json(all_sites, filename=sites_json)
    save_neighbor_site_list_to_json(cas, filename=neighbor_json)

    single_geoms, single_sites_lists = generate_all_sites(
        slab, all_sites, nslab, site_bond_cutoff, adsorbate_height
    )

    # graph clustering + label assignment
    admols, geom_indices = classify_all_sites(single_geoms, single_sites_lists)
    _, clusters = cluster_isomorphic_graphs(admols)
    _, key_to_label = update_site_labels_by_graph_and_type(
        single_sites_lists, clusters, geom_indices
    )

    # reduce: one representative geometry per distinct label
    rep_geoms, rep_sites = _reduce_to_representatives(single_geoms, single_sites_lists)

    traj = Trajectory(traj_filename, "w")
    for g in rep_geoms:
        traj.write(g)
    traj.close()

    labeled_sites_data = [
        {"geom_index": i, "sites": sites}
        for i, sites in enumerate(rep_sites)
    ]
    save_sites_to_json(labeled_sites_data, filename=labeled_sites_json)

    if verbose:
        surface_type = "CustomSurface" if use_custom else f'"{surface}"'
        print(f"Surface             : {surface_type}")
        print(f"Total ACAT sites    : {len(all_sites)}")
        print(f"Distinct labels     : {len(rep_geoms)}")
        print()
        label_to_key = {v: k for k, v in key_to_label.items()}
        print(f"  {'Label':<12} {'Original site':<14} Morphology")
        print(f"  {'-'*44}")
        for entry in rep_sites:
            for s in entry:
                label = s.get("site", "")
                _, orig_site, morph = label_to_key.get(label, (None, "", ""))
                print(f"  {label:<12} {orig_site:<14} {morph}")
        print()
        print(f"Wrote: {traj_filename}")
        print(f"Wrote: {sites_json}")
        print(f"Wrote: {neighbor_json}")
        print(f"Wrote: {labeled_sites_json}")

    return rep_geoms, rep_sites, clusters

# ============================================================
# Workflows (for clean notebooks)
# ============================================================

def workflow_detect_vacancies(slab, nslab, verbose=True):
    defect_sites, _ = detect_vacancy_sites_from_coordination(slab, nslab)
    if verbose:
        print(f"Found {len(defect_sites)} defect site(s)")
        for s in defect_sites:
            x, y, z = s["position"]
            print(f"{s.get('name','defect')}: x={x:.3f}  y={y:.3f}  z={z:.3f}  site={s.get('site')}")
    return defect_sites


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
    sites_json="sites.json",
    labeled_sites_json="labeled_sites.json",
    neighbor_json="neighbor_site_list.json",
    verbose=True,
):
    from ase.io import write

    geo = slab.copy()
    cas = SlabAdsorptionSites(slab, surface=surface_obj, composition_effect=True)
    sites = cas.get_sites()

    # all raw ACAT sites + neighbor list
    save_sites_to_json(sites, filename=sites_json)
    save_neighbor_site_list_to_json(cas, filename=neighbor_json)

    # rebuild CAS per-geometry to compute sites lists
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

    # ── collect all unique probe positions where CN > 0 ──────────────────────
    # Replace the single maxcn entry with every distinct void position that has
    # at least one slab neighbour — each becomes a candidate adsorption site.
    _NOBLE = {"He", "Ne", "Ar", "Kr", "Xe", "Rn"}

    def _probe_pos(atoms):
        idx = [i for i, a in enumerate(atoms) if a.symbol in _NOBLE]
        return atoms.positions[idx[0]].copy() if idx else None

    n_terrace       = len(geom_all) - len(maxcn_geoms)
    geom_all        = list(geom_all[:n_terrace])
    sites_lists_all = [e for e in sites_lists_all if e["geom_index"] < n_terrace]

    # ── step 1: collect unique CN>0 frames (coarse dedup by probe distance) ──
    cn_geoms, cn_site_dicts = [], []
    for f, m in zip(drop_geoms, drop_meta):
        if m.get("cn", 0) <= 0:
            continue
        pos_f = _probe_pos(f)
        if pos_f is None:
            continue
        keep = True
        for u in cn_geoms:
            pos_u = _probe_pos(u)
            if pos_u is not None:
                _, d = get_distances([pos_f], [pos_u], cell=f.cell, pbc=f.pbc)
                if float(d) < 0.25:
                    keep = False
                    break
        if not keep:
            continue
        cn_geoms.append(f.copy())
        cn_site_dicts.append({
            "site":            "defect",
            "surface":         "CustomSurface",
            "morphology":      "defect",
            "cn":              int(m["cn"]),
            "position":        pos_f.tolist(),
            "normal":          [0.0, 0.0, 1.0],
            "indices":         [],
            "composition":     "null",
            "subsurf_index":   "null",
            "subsurf_element": "null",
            "label":           f'cn{int(m["cn"])}',
            "topology":        [],
        })

    # ── step 2: group by CN, compare graphs within each group, keep median ───
    rep_geoms, rep_sites = _classify_defect_sites_by_cn_and_graph(
        cn_geoms, cn_site_dicts
    )

    offset = len(geom_all)
    for k, (g, s) in enumerate(zip(rep_geoms, rep_sites)):
        geom_all.append(g)
        sites_lists_all.append({"geom_index": offset + k, "sites": [s]})
    # ── end collect ──────────────────────────────────────────────────────────

    # build RMG graph for each geometry, cluster by isomorphism, assign labels
    single_sites_lists = [entry["sites"] for entry in sites_lists_all]
    admols, geom_indices = classify_all_sites(geom_all, single_sites_lists)
    _, clusters = cluster_isomorphic_graphs(admols)
    update_site_labels_by_graph_and_type(single_sites_lists, clusters, geom_indices)

    # reduce to one representative per distinct label → labeled_sites.json
    rep_geoms, rep_site_lists = _reduce_to_representatives(geom_all, single_sites_lists)
    rep_data = [{"geom_index": i, "sites": s} for i, s in enumerate(rep_site_lists)]
    save_sites_to_json(rep_data, filename=labeled_sites_json)

    if verbose:
        print(f"[defect] CN>0 unique void sites : {len(cn_geoms)}")
        print(f"[defect] Total sites            : {len(sites_lists_all)}")
        print(f"[defect] Distinct labels        : {len(rep_geoms)}")
        print(f"Wrote: {traj_geom_all}")
        print(f"Wrote: {traj_drop}")
        print(f"Wrote: {traj_maxcn}")
        print(f"Wrote: {sites_json}")
        print(f"Wrote: {labeled_sites_json}")
        print(f"Wrote: {neighbor_json}")

    return geom_all, sites_lists_all, maxcn_geoms, maxcn_meta, drop_geoms, drop_meta


def workflow_auto(
    xyz_path,
    n_layers=4,
    adsorbate_height=1.0,
    site_bond_cutoff=1.5,
    facet=None,            # e.g. "fcc111", "fcc332" — uses ACAT string; None → CustomSurface
    tag_symbol="Ne",
    dz=0.1,
    stable_steps=3,
    max_drop=10.0,
    margin=0.25,
    min_clearance=1.3,
    save_all_drop_steps=True,
    verbose=True,
):
    from ase.io import read
    from acat.settings import CustomSurface

    slab = read(xyz_path)
    nslab = len(slab)
    surface_obj = CustomSurface(slab, n_layers=n_layers)

    # surface passed to no-defect workflow: explicit facet string OR CustomSurface
    surface_no_defect = facet if facet is not None else surface_obj

    defect_sites = workflow_detect_vacancies(slab, nslab, verbose=verbose)

    if len(defect_sites) == 0:
        return workflow_no_defect_unique_sites(
            slab=slab,
            nslab=nslab,
            adsorbate_height=adsorbate_height,
            site_bond_cutoff=site_bond_cutoff,
            surface=surface_no_defect,
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


# ============================================================
# Visualization
# ============================================================

def visualize_labeled_sites(
    labeled_sites_json="labeled_sites.json",
    sites_json="sites.json",
    traj_file=None,
):
    """
    Load labeled_sites.json, re-run the graph pipeline, and display:
      - Figure 1 : one RMG graph per unique site label
      - Figure 2 : pairwise isomorphism matrix
      - Summary  : total ACAT sites vs. unique sites, label table

    Parameters
    ----------
    labeled_sites_json : path to labeled_sites.json written by workflow_auto
    sites_json         : path to sites.json (raw ACAT sites, for total count)
    traj_file          : traj to load geometries from; None → auto-detect
                         (tries geom_all.traj then unique_sites.traj)
    """
    import json, re
    import networkx as nx
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from ase.io import read

    # ── load geometries ───────────────────────────────────────────────────────
    if traj_file is not None:
        geoms = read(traj_file, index=":")
    else:
        try:
            geoms = read("geom_all.traj", index=":")
        except Exception:
            geoms = read("unique_sites.traj", index=":")

    # ── load labeled sites ────────────────────────────────────────────────────
    with open(labeled_sites_json) as f:
        raw = json.load(f)

    try:
        with open(sites_json) as f:
            n_total_acat = len(json.load(f))
    except Exception:
        n_total_acat = None

    single_sites_lists = [entry["sites"] for entry in raw]

    # ── graph pipeline ────────────────────────────────────────────────────────
    admols, geom_indices = classify_all_sites(geoms, single_sites_lists)
    iso_mat, clusters    = cluster_isomorphic_graphs(admols)
    _, key_to_label      = update_site_labels_by_graph_and_type(
        single_sites_lists, clusters, geom_indices
    )
    save_sites_to_json(raw, filename=labeled_sites_json)

    label_to_key = {v: k for k, v in key_to_label.items()}

    def _sort_key(lbl):
        m = re.match(r'^(.*?)(\d+)$', lbl)
        return (m.group(1), int(m.group(2))) if m else (lbl, 0)

    labels_sorted = sorted(label_to_key, key=_sort_key)

    # one representative admol per label
    label_to_admol = {}
    for graph_idx, geom_idx in enumerate(geom_indices):
        for s in single_sites_lists[geom_idx]:
            label = s.get("site")
            if label and label not in label_to_admol:
                label_to_admol[label] = admols[graph_idx]
                break

    # ── RMG mol → networkx graph ──────────────────────────────────────────────
    def rmg_to_nx(mol):
        G = nx.Graph()
        atom_to_idx = {a: i for i, a in enumerate(mol.atoms)}
        for i, atom in enumerate(mol.atoms):
            sym = atom.element.symbol if hasattr(atom.element, "symbol") else str(atom.element)
            G.add_node(i, element=sym)
        for i, atom in enumerate(mol.atoms):
            for other in atom.edges:
                j = atom_to_idx[other]
                if i < j:
                    G.add_edge(i, j)
        return G

    ELEM_COLOR = {
        "X": "#aaaaaa", "Ne": "#00bfff", "Li": "#ffd700",
        "C": "#404040", "O": "#e04040", "N": "#4040e0", "H": "#eeeeee",
    }

    def node_colors(G):
        return [ELEM_COLOR.get(G.nodes[n]["element"], "#cccccc") for n in G.nodes()]

    # ── Figure 1 : one graph per label ───────────────────────────────────────
    n      = len(labels_sorted)
    ncols  = min(5, n)
    nrows  = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    for ax, label in zip(axes_flat, labels_sorted):
        mol = label_to_admol.get(label)
        if mol is None:
            ax.set_visible(False)
            continue
        G   = rmg_to_nx(mol)
        pos = nx.kamada_kawai_layout(G)
        nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#666666", width=1.5, alpha=0.7)
        nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors(G),
                               node_size=300, edgecolors="#333333", linewidths=0.8)
        nx.draw_networkx_labels(G, pos, ax=ax,
                                labels={nd: G.nodes[nd]["element"] for nd in G.nodes()},
                                font_size=6)
        _, orig_site, morph = label_to_key[label]
        ax.set_title(f"{label}\n{orig_site} | {morph}", fontsize=8)
        ax.axis("off")

    for ax in axes_flat[n:]:
        ax.set_visible(False)

    legend_patches = [mpatches.Patch(color=c, label=e) for e, c in ELEM_COLOR.items()]
    fig.legend(handles=legend_patches, loc="lower center", ncol=len(ELEM_COLOR),
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.01))
    plt.suptitle(
        "One graph per unique site label\n"
        "(same label  ⟺  isomorphic graph + same site type + same morphology)",
        fontsize=11,
    )
    plt.tight_layout()
    plt.show()

    # ── Figure 2 : isomorphism matrix ────────────────────────────────────────
    n_g = len(admols)
    _, ax2 = plt.subplots(figsize=(max(5, n_g * 0.28), max(4, n_g * 0.28)))
    im = ax2.imshow(iso_mat.astype(float), cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax2.set_title("Graph isomorphism matrix  (dark = isomorphic pair)", fontsize=11)
    ax2.set_xlabel("Site index")
    ax2.set_ylabel("Site index")
    plt.colorbar(im, ax=ax2, shrink=0.7, label="isomorphic (1=yes)")
    plt.tight_layout()
    plt.show()

    # ── Summary ───────────────────────────────────────────────────────────────
    n_unique = len(labels_sorted)
    print("=" * 54)
    if n_total_acat is not None:
        print(f"  Total ACAT sites (sites.json)    : {n_total_acat}")
    print(f"  Unique sites     (labeled_sites) : {n_unique}")
    print()
    print("  Each distinct label = one unique site environment.")
    print("  'Unique sites' = distinct labels (same thing).")
    print("=" * 54)
    print()
    print(f"  {'Label':<12} {'Original site':<14} Morphology")
    print(f"  {'─' * 44}")
    for label in labels_sorted:
        _, site_val, morph_val = label_to_key[label]
        print(f"  {label:<12} {site_val:<14} {morph_val}")

# Add to site_analysis.py (reuses classify_all_sites + cluster_isomorphic_graphs).
# Shape = site type (ontop=circle, bridge=square, 3fold=triangle);
# colour = equivalence (grey=unique, coloured+number=shared via graph isomorphism).
# Legends sit at the bottom in two rows, like the ACAT scatter plots.

def plot_site_equivalence_xy(
    labeled_sites_json="labeled_sites.json",
    traj_file=None,
    site_types=("ontop", "bridge", "3fold", "4fold", "defect"),
    clusters=None,
    geom_indices=None,
    single_sites_lists=None,
    title="All sites — same vs different (graph isomorphism)",
    save_path="site_equivalence.png",
):
    """
    Single-panel same-vs-different map.
      grey marker            = unique environment (different from every other
                               site of its type)
      coloured marker        = shares its environment with same-shape markers
                               of the same colour (these are 'the same' site)
    'Same' = isomorphic slab+probe graphs (classify_all_sites ->
    cluster_isomorphic_graphs). Pass clusters / geom_indices /
    single_sites_lists to reuse values from visualize_labeled_sites.

    Returns report dict {site_type: {"n_sites","n_distinct","duplicate_groups"}}.
    """
    import json
    from collections import defaultdict, Counter
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    from ase.io import read

    def _stype(name):
        return (name or "").rstrip("0123456789")

    MARKER = {"ontop": "o", "bridge": "s", "3fold": "^", "4fold": "D", "defect": "X"}

    if single_sites_lists is None:
        with open(labeled_sites_json) as f:
            raw = json.load(f)
        single_sites_lists = [entry["sites"] for entry in raw]

    if clusters is None or geom_indices is None:
        if traj_file is not None:
            geoms = read(traj_file, index=":")
        else:
            try:
                geoms = read("geom_all.traj", index=":")
            except Exception:
                geoms = read("unique_sites.traj", index=":")
        admols, geom_indices = classify_all_sites(geoms, single_sites_lists)
        _, clusters = cluster_isomorphic_graphs(admols)

    graph_to_cluster = {}
    for cid, members in enumerate(clusters.values()):
        for graph_idx in members:
            graph_to_cluster[graph_idx] = cid
    geom_to_cluster = {geom_indices[gi]: cid for gi, cid in graph_to_cluster.items()}

    records = []  # (type, label, x, y, cluster_id)
    for geom_idx, sites in enumerate(single_sites_lists):
        cid = geom_to_cluster.get(geom_idx)
        if cid is None:
            continue
        for s in sites:
            if s.get("site"):
                x, y = s["position"][0], s["position"][1]
                records.append((_stype(s["site"]), s["site"], x, y, cid))
                break

    present = [t for t in site_types if any(r[0] == t for r in records)]
    tc_counts = Counter((t, cid) for (t, _, _, _, cid) in records)
    GROUP_COLORS = (list(plt.cm.tab10.colors) + list(plt.cm.Set2.colors)
                    + list(plt.cm.tab20.colors))

    fig, ax = plt.subplots(figsize=(7.5, 7))
    gi = 0
    group_color = {}
    for (t, lab, x, y, cid) in records:
        if t not in present:
            continue
        mk = MARKER.get(t, "P")
        if tc_counts[(t, cid)] == 1:
            ax.scatter(x, y, marker=mk, color="#d9d9d9", s=45,
                       alpha=0.9, edgecolors="none", zorder=2)
        else:
            if (t, cid) not in group_color:
                group_color[(t, cid)] = GROUP_COLORS[gi % len(GROUP_COLORS)]
                gi += 1
            ax.scatter(x, y, marker=mk, color=group_color[(t, cid)], s=45,
                       alpha=0.95, edgecolors="none", zorder=3)

    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)")
    ax.set_aspect("equal")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(title, fontsize=12, fontweight="bold")

    # ── two bottom rows: site type, then equivalence (like ACAT plots) ────────
    type_handles = [mlines.Line2D([], [], marker=MARKER.get(t, "P"), color="w",
                    markerfacecolor="#bbb", markeredgecolor="#555",
                    markersize=9, label=t) for t in present]
    state_handles = [
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor="#d9d9d9",
                      markeredgecolor="#999", markersize=9, label="unique (different)"),
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor="#e0709a",
                      markeredgecolor="none", markersize=9, label="shared (same colour = same site)"),
    ]
    fig.legend(handles=type_handles, loc="lower center", ncol=len(present),
               frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.02))
    fig.legend(handles=state_handles, loc="lower center", ncol=2,
               frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.08))
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()

    report = {}
    for t in present:
        by_cluster = defaultdict(list)
        for (st, lab, _, _, cid) in records:
            if st == t:
                by_cluster[cid].append(lab)
        dup = [labs for labs in by_cluster.values() if len(labs) > 1]
        report[t] = {"n_sites": sum(len(v) for v in by_cluster.values()),
                     "n_distinct": len(by_cluster), "duplicate_groups": dup}

    print("=" * 60)
    for t, info in report.items():
        if not info["duplicate_groups"]:
            print(f"  {t:7s}: all {info['n_sites']} sites are DISTINCT.")
        else:
            print(f"  {t:7s}: {info['n_sites']} sites -> {info['n_distinct']} distinct; "
                  f"{len(info['duplicate_groups'])} group(s) are 'the same':")
            for g in info["duplicate_groups"]:
                print(f"           {' = '.join(g)}")
    print("=" * 60)
    print(f"Saved: {save_path}")
    return report