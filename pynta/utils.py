import shutil
import os
import ase
from ase import Atoms 
from molecule.quantity import ScalarQuantity
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from ase.io import write, read
import ase.constraints
from ase.geometry import get_distances
from ase.calculators.mixing import SumCalculator
from copy import deepcopy
from importlib import import_module
import numpy as np
import copy
import json

def sites_match(site1,site2,slab,tol=0.5):
    """determine if two sites match

    Args:
        site1 (dict): site dictionary for first site
        site2 (dict): site dictionary for second site
        slab (ase.Atoms): the ase.Atoms object for the slab the sites are on
        tol (float, optional): Distance tolerance defaults to 0.5 Angstroms.

    Returns:
        _type_: _description_
    """
    v,dist = get_distances([site1["position"]], [site2["position"]], cell=slab.cell, pbc=slab.pbc)
    if dist > tol:
        return False
    elif site1["site"] != site2["site"]:
        return False
    elif site1["morphology"] != site2["morphology"]:
        return False
    else:
        return True
    
def get_occupied_sites(struct,sites,nslab,allowed_site_dict=dict(),site_bond_cutoff=2.45,
                       site_bond_disruption_cutoff=0.5,cutoff_corrections={"ontop": 0.3}):
    """determine what sites are occupied by what atoms

    Args:
        struct (ase.Atoms): Atoms object for the structure
        sites (list): list of site dictionaries
        nslab (int): number of atoms in the surface
        allowed_site_dict (dict, optional): dictionary mapping atom index to a list of allowed (site,morphology) for that atom
        site_bond_cutoff (float, optional): _description_. Defaults to 2.5.

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    occ_sites = []
    for i in range(nslab,len(struct)):
        pos = struct.positions[i]
        if i in allowed_site_dict.keys():
            allowed_sites = allowed_site_dict[i]
        else:
            allowed_sites = None
        siteout = None
        mindist = None
        for site in sites:
            if allowed_sites and (site["site"],site["morphology"]) not in allowed_sites: #skip disallowed sites
                continue
            v,dist = get_distances([site["position"]], [pos], cell=struct.cell, pbc=struct.pbc)
            if siteout is None:
                siteout = site
                mindist = dist
                n = v
            else:
                if dist < mindist:
                    mindist = dist
                    siteout = site
                    n = v
        
        if mindist is None:
            #print(i)
            #view(struct)
            raise ValueError
        
        if (siteout["site"] not in cutoff_corrections.keys() and mindist < site_bond_cutoff) or (siteout["site"] in cutoff_corrections.keys() and mindist < site_bond_cutoff + cutoff_corrections[siteout["site"]]):
            mindn = None
            for j in range(nslab,len(struct)): #check for site bond disruption by other adsorbed atoms
                if i == j:
                    continue
                else:
                    AB,ABdist = get_distances([struct.positions[j]], [pos], cell=struct.cell, pbc=struct.pbc)
                    x = np.dot(AB[0,0,:],n[0,0,:])
                    if x < 0: #this means our target atom is closer than the other atom, which is okay (for this atom)
                        continue
                    dn = np.linalg.norm(AB[0,0,:] - x/mindist**2*n[0,0,:])
                    if mindn is None:
                        mindn = dn
                    elif dn < mindn:
                        mindn = dn
                    
                
            if mindn is None or mindn > site_bond_disruption_cutoff:
                s = deepcopy(siteout)
                s["bonding_index"] = i
                s["bond_length"] = mindist
                occ_sites.append(s)
    
    return occ_sites

class SiteOccupationException(Exception):
    pass

def get_unique_sym(geoms):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of paths to .xyz or .traj files to compare

    Returns
    ________
    idx_list : list(str)
        a list with prefixes of all symmetrically distinct sites

    '''
    comparator = SymmetryEquivalenceCheck()

    good_adsorbates_atom_obj_list = []
    geos_out = []

    for geom in geoms:
        adsorbate_atom_obj = read(geom)
        adsorbate_atom_obj.set_constraint()
        adsorbate_atom_obj.pbc = True
        adsorbate_atom_obj.set_constraint() # Reset constraints for comparison
        comparision = comparator.compare(
            adsorbate_atom_obj, good_adsorbates_atom_obj_list)

        if comparision is False:
            good_adsorbates_atom_obj_list.append(adsorbate_atom_obj)
            geos_out.append(geom)

    return geos_out

def get_unique_sym_indices(geoms):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of paths to .xyz or .traj files to compare

    '''
    comparator = SymmetryEquivalenceCheck()

    good_adsorbates_atom_obj_list = []
    geos_out = []

    for geom in geoms:
        adsorbate_atom_obj = read(geom)
        adsorbate_atom_obj.set_constraint()
        adsorbate_atom_obj.pbc = True
        adsorbate_atom_obj.set_constraint() # Reset constraints before comparison
        comparision = comparator.compare(
            adsorbate_atom_obj, good_adsorbates_atom_obj_list)

        if comparision is False:
            good_adsorbates_atom_obj_list.append(adsorbate_atom_obj)
            geos_out.append(geom)

    indices = [geoms.index(g) for g in geos_out]

    return indices

def get_unique_sym_structs(geoms):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of Atoms objects to compare


    '''
    comparator = SymmetryEquivalenceCheck()

    geoms_copy = deepcopy(geoms)

    good_adsorbates_atom_obj_list = []
    geos_out = []

    for i,geom in enumerate(geoms_copy):
        adsorbate_atom_obj = geom
        adsorbate_atom_obj.set_constraint()
        adsorbate_atom_obj.pbc = True
        adsorbate_atom_obj.set_constraint() # Reset constraints before comparison
        comparision = comparator.compare(
            adsorbate_atom_obj, good_adsorbates_atom_obj_list)

        if comparision is False:
            good_adsorbates_atom_obj_list.append(adsorbate_atom_obj)
            geos_out.append(geoms[i])

    return geos_out

def get_unique_sym_struct_indices(geoms):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of Atoms objects to compare


    '''
    comparator = SymmetryEquivalenceCheck()

    geoms_copy = deepcopy(geoms)

    good_adsorbates_atom_obj_list = []
    indices = []

    for i,geom in enumerate(geoms_copy):
        adsorbate_atom_obj = geom
        adsorbate_atom_obj.set_constraint()
        adsorbate_atom_obj.pbc = True
        adsorbate_atom_obj.set_constraint() # Reset constraints before comparison
        comparision = comparator.compare(
            adsorbate_atom_obj, good_adsorbates_atom_obj_list)

        if comparision is False:
            good_adsorbates_atom_obj_list.append(adsorbate_atom_obj)
            indices.append(i)

    return indices

def get_unique_sym_struct_index_clusters(geoms):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of Atoms objects to compare


    '''
    comparator = SymmetryEquivalenceCheck()

    geoms_copy = deepcopy(geoms)

    good_adsorbates_atom_obj_list = []
    indices = []

    for i,geom in enumerate(geoms_copy):
        adsorbate_atom_obj = geom
        adsorbate_atom_obj.set_constraint()
        adsorbate_atom_obj.pbc = True
        adsorbate_atom_obj.set_constraint() # Reset constraints before comparison
        comparison = None
        for j,adlist in enumerate(good_adsorbates_atom_obj_list):
            comparison = comparator.compare(adsorbate_atom_obj, [adlist[0]])
            ind = j
            if comparison:
                break
        else:
            comparison = False

        if comparison is False:
            good_adsorbates_atom_obj_list.append([adsorbate_atom_obj])
            indices.append([i])
        else:
            indices[ind].append(i)

    return indices

def filter_nonunique_TS_guess_indices(geoms,Es):
    ''' Check for the symmetry equivalent structures in the given files

    Parameters
    ___________
    geoms: list of paths to .xyz or .traj files to compare

    '''
    comparator = SymmetryEquivalenceCheck()

    good_adsorbates_atom_obj_list = []
    geos_out = []
    Esout = []

    for j,geom in enumerate(geoms):
        adsorbate_atom_obj = read(geom)
        adsorbate_atom_obj.set_constraint()
        adsorbate_atom_obj.pbc = True
        adsorbate_atom_obj.set_constraint() # Reset constraints before comparison
        for i,good_adsorbate in enumerate(good_adsorbates_atom_obj_list):
            comparison = comparator.compare(adsorbate_atom_obj,good_adsorbate)
            if comparison and Es[j] < Esout[i]:
                geos_out[i] = geom
                good_adsorbates_atom_obj_list[i] = adsorbate_atom_obj
                Esout[i] = Es[j]
                break
            elif comparison:
                break
        else:
            good_adsorbates_atom_obj_list.append(adsorbate_atom_obj)
            geos_out.append(geom)
            Esout.append(Es[j])

    return geos_out,Esout

def get_fmax(at):
    return np.max(np.abs([np.linalg.norm(at.get_forces()[i,:]) for i in range(at.get_forces().shape[0])]))

def name_to_ase_software(software_name,module_name=None):
    """
    go from software_name to the associated
    ASE calculator constructor
    """
    if isinstance(software_name,list):
        return [name_to_ase_software(x) for x in software_name]
    
    if module_name is None:
        module_name = software_name
        
    if software_name == "TBLite" or software_name == "XTB":
        module = import_module("tblite.ase")
        return getattr(module, "TBLite")
    elif software_name == "MACECalculator":
        module = import_module("mace.calculators.mace")
        return getattr(module, software_name)
    elif software_name == "TorchDFTD3Calculator":
        module = import_module("torch_dftd.torch_dftd3_calculator")
        return getattr(module, software_name)
    else:
        module = import_module("ase.calculators."+module_name.lower())
        return getattr(module, software_name)

def to_ase_software(object_name,software_kwargs,module_name=None):
    if module_name is None:
        module_name = object_name
    
    software_kwargs = copy.deepcopy(software_kwargs)
    software = name_to_ase_software(object_name,module_name=module_name)
    
    if not isinstance(software,list):
        for k,v in software_kwargs.items():
            if isinstance(v,dict) and "type" in v.keys():
                typ = v["type"]
                dv = {k:v for k,v in v.items() if k != "type"}
                software_kwargs[k] = to_ase_software(typ,dv,module_name=module_name)
        
        return software(**software_kwargs) 
    else:
       return SumCalculator([to_ase_software(object_name[i],software_kwargs[i],module_name=module_name[i]) if module_name is not None else to_ase_software(object_name[i],software_kwargs[i]) for i in range(len(object_name))])

def name_to_ase_opt(opt_name):
    """
    go from the optimizer name to the
    ASE optimizer
    """
    module = import_module("ase.optimize")
    return getattr(module, opt_name)

def clean_pynta_path(path,save_initial_guess=True):
    if save_initial_guess:
        for p in os.listdir(path):
            if p[:2] == "TS": #delete TSs
                shutil.rmtree(os.path.join(path,p))
            elif p == "Adsorbates":
                for ad in os.listdir(os.path.join(path,p)):
                    if ad == ".DS_Store":
                        os.remove(os.path.join(path,p,ad))
                        continue
                    for ind in os.listdir(os.path.join(path,p,ad)):
                        if ind == "info.json":
                            continue
                        elif ind.isdigit():
                            for file in os.listdir(os.path.join(path,p,ad,ind)):
                                if not "_init.xyz" in file:
                                    pa = os.path.join(path,p,ad,ind,file)
                                    if os.path.isdir(pa):
                                        shutil.rmtree(pa)
                                    else:
                                        os.remove(pa)
                        else:
                            os.remove(os.path.join(path,p,ad,ind))
    else:
         for p in os.listdir(path):
            if p[:2] == "TS" or p == "Adsorbates": #delete TSs
                shutil.rmtree(os.path.join(path,p))
    
    if os.path.exists(os.path.join(path,"reaction_library")):
        shutil.rmtree(os.path.join(path,"reaction_library"))
    if os.path.exists(os.path.join(path,"thermo_library.py")):
        os.remove(os.path.join(path,"thermo_library.py"))

def save_slab_sites(path, sites, site_adjacency):
    """Persist the ACAT-derived sites + site_adjacency for a slab to JSON in ``path`` (sites.json,
    site_adjacency.json), so coverage dependence and postprocessing reuse the EXACT same set (and
    ordering) instead of re-deriving them with ACAT -- which can reorder sites or fail get_labels on
    open/stepped surfaces. numpy arrays/tuples are made JSON-safe; load with load_slab_sites."""
    import json
    def _jsonify(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, dict):
            return {k: _jsonify(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [_jsonify(v) for v in o]
        return o
    with open(os.path.join(path, "sites.json"), "w") as f:
        json.dump(_jsonify(sites), f)
    with open(os.path.join(path, "site_adjacency.json"), "w") as f:
        json.dump({int(k): [int(x) for x in v] for k, v in site_adjacency.items()}, f)

def load_slab_sites(path):
    """Load sites + site_adjacency saved by save_slab_sites (falling back to the legacy
    single_sites_lists.json / neighbor_site_list.json names if present), converting position/normal
    back to numpy arrays and indices to tuples. Use this everywhere instead of re-running
    SlabAdsorptionSites so the whole pipeline shares one site definition. Returns (sites, site_adjacency)."""
    import json
    sites_file = os.path.join(path, "sites.json")
    if not os.path.exists(sites_file) and os.path.exists(os.path.join(path, "single_sites_lists.json")):
        sites_file = os.path.join(path, "single_sites_lists.json")
    adj_file = os.path.join(path, "site_adjacency.json")
    if not os.path.exists(adj_file) and os.path.exists(os.path.join(path, "neighbor_site_list.json")):
        adj_file = os.path.join(path, "neighbor_site_list.json")
    if not os.path.exists(sites_file) or not os.path.exists(adj_file):
        raise FileNotFoundError(
            "no saved sites in {} (looked for sites.json/site_adjacency.json and the legacy "
            "single_sites_lists.json/neighbor_site_list.json). Pass sites/site_adjacency explicitly, "
            "or re-run the isolated workflow with current pynta so it saves them.".format(path))
    with open(sites_file) as f:
        sites = json.load(f)
    for s in sites:
        s["position"] = np.array(s["position"])
        if s.get("normal") is not None:
            s["normal"] = np.array(s["normal"])
        if s.get("indices") is not None:
            s["indices"] = tuple(s["indices"])
    with open(adj_file) as f:
        site_adjacency = {int(k): v for k, v in json.load(f).items()}
    return sites, site_adjacency

def write_sites_xyz(slab, sites, out_xyz="sites.xyz", marker_height=0.0, by="site"):
    """Write the slab + one marker atom per adsorption site to an xyz so the sites can be viewed
    spatially (e.g. to check whether nearby sites are genuinely distinct or near-duplicates). The
    marker element encodes the site category (by="site" | "morphology" | "label") so a viewer colors
    them; a legend (marker element -> category) is printed. marker_height raises markers in z.
    Returns (atoms, legend)."""
    marker_pool = ["H","He","Li","Be","B","C","N","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc"]
    def _cat(s):
        if by == "morphology":
            return s.get("morphology", "?")
        if by == "label":
            return str(s.get("label", "?"))
        return s.get("site", "?")
    cats = [_cat(s) for s in sites]
    legend = {c: marker_pool[i % len(marker_pool)] for i, c in enumerate(sorted(set(cats)))}
    symbols = [legend[c] for c in cats]
    positions = []
    for s in sites:
        p = np.array(s["position"], dtype=float).copy()
        p[2] += marker_height
        positions.append(p)
    atoms = slab.copy()
    if positions:
        atoms += Atoms(symbols=symbols, positions=positions)
    write(out_xyz, atoms)
    print("wrote {} (slab + {} site markers). legend (marker element -> {} category):".format(out_xyz, len(sites), by))
    for c, e in sorted(legend.items(), key=lambda kv: kv[1]):
        print("  {:3s} -> {}".format(e, c))
    return atoms, legend

def _interaction_terms(admol, tree, sites):
    """Decompose a 2D config admol into its SIDT interaction terms for visualization. Returns
    (ad_xy, terms): ad_xy maps a surface-bonded adsorbate atom -> its site xy; terms is a list of
    (atom_i, atom_j, contribution, node_id), one per pair-centered subgraph in the SAME order as the
    tree's decomposition/trace (contribution = matched node value, J/mol). Works directly on the 2D
    admol (no 3D round-trip), so adsorbates close in space are not merged. Assumes the admol was built
    with max_dist=inf, where generate_adsorbate_molecule keeps ninds=range(len(sites)), so the admol's
    k-th surface-site atom maps to sites[k]."""
    import logging
    from pynta.mol import find_adsorbate_atoms_surface_sites

    # max_dist=inf admols: k-th surface-site atom <-> sites[k]
    site_atoms = [a for a in admol.atoms if a.is_surface_site()]
    if len(site_atoms) != len(sites):
        logging.warning("interaction_terms: admol has %d surface sites but sites has %d; xy mapping "
                        "assumes the admol was built with max_dist=inf", len(site_atoms), len(sites))
    site_xy = {id(a): np.asarray(sites[k]["position"])[:2]
               for k, a in enumerate(site_atoms) if k < len(sites)}

    def ad_xy(atom):
        ps = [site_xy[id(s)] for s in atom.bonds.keys() if s.is_surface_site() and id(s) in site_xy]
        return np.mean(ps, axis=0) if ps else None

    # pairs in the SAME order as adsorbate_interaction_decomposition (i>j, skip same-adsorbate)
    surf_bonded = [a for a in admol.atoms if a.is_bonded_to_surface() and not a.is_surface_site()]
    sb_inds = [admol.atoms.index(a) for a in surf_bonded]
    pairs = []
    for i, indi in enumerate(sb_inds):
        for j, indj in enumerate(sb_inds):
            if i > j:
                ad_atoms_j, _ = find_adsorbate_atoms_surface_sites(admol.atoms[indj], admol)
                if admol.atoms[indi] in ad_atoms_j:
                    continue
                pairs.append((indi, indj))

    E, sigma, trace = tree.evaluate(admol, trace=True, estimate_uncertainty=True)
    if len(trace) != len(pairs):
        logging.warning("interaction_terms: trace length %d != pair count %d; edge mapping may be off",
                        len(trace), len(pairs))
    terms = []
    for (indi, indj), g in zip(pairs, trace):
        try:
            c = float(tree.nodes[g].rule.value)
        except Exception:
            c = float("nan")
        terms.append((admol.atoms[indi], admol.atoms[indj], c, g))
    return ad_xy, terms

def plot_interaction_graph(admol, tree, sites, site_adjacency, slab,
                           out_png=None, ax=None, cmap="coolwarm", label_values=False, title=None):
    """Overlay the SIDT interaction-energy decomposition of a 2D config admol on the xy site
    projection. The underlying site graph is drawn light gray; markers = surface-bonded adsorbate
    atoms; one edge per pair-centered interaction term, width proportional to |contribution| and
    color by sign (diverging cmap). Works directly from the 2D admol (no 3D->2D round-trip, so
    adsorbates close in space are not merged). Axes span the full slab cell. Contributions in eV.
    Returns (fig, ax).

    admol: 2D config admol (e.g. from build_admol_from_occ); tree: trained interaction regressor;
    slab: ase.Atoms for the slab (cell extent + PBC for the site-graph underlay).
    """
    import matplotlib.pyplot as plt
    from ase.geometry import get_distances

    ad_xy, terms = _interaction_terms(admol, tree, sites)
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))
    else:
        fig = ax.figure

    # --- underlying site graph (light gray), minimum-image segments to avoid cross-cell lines ---
    site_pos = [np.asarray(s["position"]) for s in sites]
    for i, neighbors in site_adjacency.items():
        pi = site_pos[i][:2]
        for j in neighbors:
            v, _ = get_distances([site_pos[i]], [site_pos[j]], cell=slab.cell, pbc=slab.pbc)
            end = pi + v[0][0][:2]
            ax.plot([pi[0], end[0]], [pi[1], end[1]], color="0.85", lw=0.6, zorder=0)
    ax.scatter([p[0] for p in site_pos], [p[1] for p in site_pos], s=8, c="0.8", zorder=1)

    # --- interaction terms (eV) ---
    EV = 1.0 / (96.48530749925793 * 1000.0)  # J/mol -> eV
    contribs = [t[2] * EV for t in terms if not np.isnan(t[2])]
    maxc = max((abs(c) for c in contribs), default=1.0) or 1.0
    norm = plt.Normalize(-maxc, maxc)
    cm = plt.get_cmap(cmap)

    for ai, aj, c, g in terms:
        if np.isnan(c):
            continue
        c = c * EV
        p, q = ad_xy(ai), ad_xy(aj)
        if p is None or q is None:
            continue
        ax.plot([p[0], q[0]], [p[1], q[1]], color=cm(norm(c)),
                linewidth=0.5 + 5.0 * abs(c) / maxc, alpha=0.85, zorder=2)
        if label_values:
            mid = (p + q) / 2.0
            ax.annotate("{:.3f}".format(c), mid, fontsize=6, zorder=5)

    for a in admol.atoms:
        if a.is_bonded_to_surface() and not a.is_surface_site():
            xy = ad_xy(a)
            if xy is not None:
                ax.scatter(xy[0], xy[1], s=120, c="k", zorder=3)
                ax.annotate(a.element.symbol, xy, color="white", ha="center", va="center",
                            fontsize=7, zorder=4)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="pair interaction contribution (eV)")

    # axes span the full slab cell (in-plane bounding box)
    a_xy, b_xy = np.asarray(slab.cell[0])[:2], np.asarray(slab.cell[1])[:2]
    corners = np.array([[0.0, 0.0], a_xy, b_xy, a_xy + b_xy])
    mx = 0.05 * (corners[:, 0].max() - corners[:, 0].min() + 1e-9)
    my = 0.05 * (corners[:, 1].max() - corners[:, 1].min() + 1e-9)
    ax.set_xlim(corners[:, 0].min() - mx, corners[:, 0].max() + mx)
    ax.set_ylim(corners[:, 1].min() - my, corners[:, 1].max() + my)
    ax.set_aspect("equal")
    ax.set_xlabel("x [A]")
    ax.set_ylabel("y [A]")
    total = sum(contribs)
    ax.set_title(title or "interaction E = {:.3f} eV ({} terms)".format(total, len(terms)))
    if out_png:
        fig.savefig(out_png, dpi=150, bbox_inches="tight")
    return fig, ax

def construct_constraint(d):
    """
    construct a constrain from a dictionary that is the input to the constraint
    constructor plus an additional "type" key that indices the name of the constraint
    returns the constraint
    """
    constraint_dict = copy.deepcopy(d)
    constructor = getattr(ase.constraints,constraint_dict["type"])
    del constraint_dict["type"]
    return constructor(**constraint_dict)

def to_dict(self):
    out_dict = dict()
    attrs = [attr for attr in dir(self) if not attr.startswith("_")]
    for attr in attrs:
        val = getattr(self,attr)
    
        if not isinstance(val,list) and not isinstance(val,np.ndarray) and not isinstance(val,dict) and callable(val):
            continue
        
        try:
            json.dumps(val)
            out_dict[attr] = val 
            assert not isinstance(val,np.ndarray)
        except:
            if isinstance(val, ScalarQuantity):
                out_dict[attr] = {
                "class": val.__class__.__name__,
                "value": val.value,
                "units": val.units,
                "uncertainty": val.uncertainty,
                "uncertainty_type": val.uncertainty_type,
            }
            elif isinstance(val, Atoms):
                new_dict = val.todict()
                new_dict["class"] = val.__class__.__name__
                out_dict[attr] = new_dict
            else:
                out_dict[attr] = to_dict
                
    return out_dict