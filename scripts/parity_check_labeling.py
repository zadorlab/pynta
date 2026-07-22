"""
Parity check + speedup measurement for the positional site-labeling change in
`adsorbate_interaction_decomposition` + `get_adsorbed_atom_pairs`
(pynta/coveragedependence.py).

The labeling tags each surface site by its shortest site-hop distance to the nearest
adsorbate endpoint (endpoints stay "*"). It is meant to be a *speed* change only --
same tree, same predictions. Its one real behavior change is that the labels reject
"winding" (non-geodesic) subgraph matches, making separation-binning geodesic-exact.
This script measures the training speedup AND checks whether predictions moved.

It trains the tree twice in one process, on the NEW code (no git toggling):
  * "new"  = current labeled behavior
  * "old"  = same code with the numbered "*1"/"*2"/... labels stripped back off
             (monkeypatched) -> reproduces the pre-change unlabeled-site behavior
then compares node counts, training time, and per-datum predictions.

--------------------------------------------------------------------------------
STANDALONE USAGE (recommended) -- point it at a covdep run folder on the cluster:

    conda activate pynta_env
    python scripts/parity_check_labeling.py /home/jzador/Pt/hox_mc_6 0

  arg1 = run path (the folder containing pairs_datums.json)
  arg2 = iteration to reproduce (default 0; iter 0 is where you saw no speedup)
  arg3 = tol in J/mol for flagging changed predictions (default 1.0)

It loads pairs_datums.json (+ Iterations/<iter>/cumulative_sample_datums.json for
iter>0), rebuilds the datums exactly as tasks.py does, and runs the comparison.

IN-DRIVER USAGE -- if you'd rather call it where the datums are already in scope:

    import sys; sys.path.insert(0, "scripts")
    from parity_check_labeling import run
    run(pairs_datums, sampling_datums, node_fract_training=0.7)
--------------------------------------------------------------------------------

Interpreting the printout:
  * similar node count + tiny prediction MAD + speedup > 1  -> pure speedup, keep it.
  * large node/prediction shift                            -> winding-path matches
                                                              were load-bearing; inspect.
"""
import sys
import os
import json
import time
import re

_LABELED_SITE = re.compile(r"^\*\d+$")  # matches "*1", "*2", ... but NOT bare "*"


def load_datums(run_path, iter=0):
    """Rebuild (pairs_datums, sampling_datums) from a run folder, matching tasks.py."""
    from molecule.molecule import Molecule
    from pysidt.sidt import Datum

    def _load(path):
        with open(path) as f:
            return [Datum(mol=Molecule().from_adjacency_list(d["mol"], check_consistency=False),
                          value=d["value"]) for d in json.load(f)]

    pairs = _load(os.path.join(run_path, "pairs_datums.json"))
    if iter == 0:
        sampling = []
    else:
        sampling = _load(os.path.join(run_path, "Iterations", str(iter),
                                      "cumulative_sample_datums.json"))
    return pairs, sampling


def _strip_numbered_labels(atoms):
    """Clear "*1"/"*2"/... labels (leave bare "*" endpoints): emulates pre-change behavior."""
    for at in atoms:
        lbl = getattr(at, "label", "")
        if lbl and _LABELED_SITE.match(lbl):
            at.label = ""


def _n_surface_sites(mol):
    try:
        return sum(1 for a in mol.atoms if a.is_surface_site())
    except Exception:
        return -1


def _predict(tree, mol):
    try:
        return float(tree.evaluate(mol, estimate_uncertainty=False))
    except TypeError:
        return float(tree.evaluate(mol))


def _train_and_eval(pairs_datums, sampling_datums, train_kwargs, emulate_old):
    """Train one tree; if emulate_old, monkeypatch the two functions to strip site labels."""
    import pynta.coveragedependence as cd

    real_pairs = cd.get_adsorbed_atom_pairs
    real_decomp = cd.adsorbate_interaction_decomposition

    if emulate_old:
        def old_pairs(*a, **k):
            groups = real_pairs(*a, **k)
            for g in groups:
                _strip_numbered_labels(g.atoms)
            return groups

        def old_decomp(mol, *a, **k):
            structs = real_decomp(mol, *a, **k)
            for st in structs:
                _strip_numbered_labels(st.atoms)
            return structs

        cd.get_adsorbed_atom_pairs = old_pairs
        cd.adsorbate_interaction_decomposition = old_decomp

    try:
        t0 = time.perf_counter()
        tree = cd.train_sidt_cov_dep_regressor(pairs_datums, sampling_datums, **train_kwargs)
        seconds = time.perf_counter() - t0
    finally:
        cd.get_adsorbed_atom_pairs = real_pairs
        cd.adsorbate_interaction_decomposition = real_decomp

    datums = list(pairs_datums) + list(sampling_datums)
    preds = []
    for d in datums:
        try:
            preds.append(_predict(tree, d.mol))
        except Exception:
            preds.append(None)
    return tree, seconds, preds, datums


def run_new_only(pairs_datums, sampling_datums, **train_kwargs):
    """Train ONLY the new (labeled) tree, timed, with a prediction sanity check.

    Skips the slow old-emulated baseline -- use when you just want the speed number
    (compare `train s` to your known old iter0 time) and a quick did-it-break check.
    """
    print("[parity] training NEW only (labeled) ...")
    new_tree, new_s, new_p, datums = _train_and_eval(pairs_datums, sampling_datums, train_kwargs, False)

    ok = [p for p in new_p if p is not None]
    failed = sum(1 for p in new_p if p is None)
    print("=" * 62)
    print(f"nodes    : {len(new_tree.nodes)}")
    print(f"train s  : {new_s:.1f}   <- compare to your known OLD iter0 time")
    print(f"datums   : {len(datums)} total, {failed} failed to evaluate")
    if ok:
        print(f"pred J/mol: min {min(ok):.1f}  mean {sum(ok)/len(ok):.1f}  max {max(ok):.1f}")
    print("=" * 62)
    if failed:
        print("WARNING: some datums failed to evaluate -- possible root/label mismatch; "
              "run the full parity (drop --new-only) to inspect.")
    elif not ok:
        print("WARNING: NO datums evaluated -- datums likely detached from root "
              "(unlabeled root_pair not matching labeled sites). Needs fixing.")
    else:
        print("sanity OK: tree built and all datums evaluate. For semantic parity vs old, "
              "run without --new-only when you can spare the time.")
    return {"new_nodes": len(new_tree.nodes), "new_s": new_s, "failed": failed}


def run(pairs_datums, sampling_datums, tol=1.0, **train_kwargs):
    """Train old-emulated vs new in one process and print the parity + speedup report."""
    print("[parity] training OLD-emulated (labels stripped) ...")
    old_tree, old_s, old_p, datums = _train_and_eval(pairs_datums, sampling_datums, train_kwargs, True)
    print(f"[parity]   old: {len(old_tree.nodes)} nodes, {old_s:.1f}s")

    print("[parity] training NEW (labeled) ...")
    new_tree, new_s, new_p, _ = _train_and_eval(pairs_datums, sampling_datums, train_kwargs, False)
    print(f"[parity]   new: {len(new_tree.nodes)} nodes, {new_s:.1f}s")

    diffs = []
    for i, (a, b, d) in enumerate(zip(old_p, new_p, datums)):
        if a is None or b is None:
            continue
        diffs.append((abs(b - a), i, a, b, float(d.value), _n_surface_sites(d.mol)))
    diffs.sort(reverse=True)

    n = len(diffs)
    mad = sum(x[0] for x in diffs) / n if n else float("nan")
    changed = sum(1 for x in diffs if x[0] > tol)
    spd = old_s / new_s if new_s else float("nan")

    print("=" * 62)
    print(f"nodes    : old {len(old_tree.nodes):>6}   new {len(new_tree.nodes):>6}   "
          f"(delta {len(new_tree.nodes) - len(old_tree.nodes):+d})")
    print(f"train s  : old {old_s:>6.1f}   new {new_s:>6.1f}   (speedup {spd:.2f}x)")
    print(f"datums   : {n} comparable")
    print(f"pred MAD : {mad:.3f} J/mol      changed >{tol} J/mol: {changed} / {n}")
    print("worst per-datum prediction shifts (J/mol):")
    for dd, i, a, b, act, ns in diffs[:10]:
        print(f"  datum {i:>5}  |dpred| {dd:10.2f}   old {a:10.2f}  new {b:10.2f}  "
              f"actual {act:10.2f}  n_sites {ns}")
    node_ok = abs(len(new_tree.nodes) - len(old_tree.nodes)) <= max(2, 0.02 * len(old_tree.nodes))
    print("=" * 62)
    print("verdict:", "PURE SPEEDUP (predictions ~unchanged)" if (mad <= tol and node_ok)
          else "BEHAVIOR CHANGED -- inspect the winding-path datums above before keeping")
    return {"old_nodes": len(old_tree.nodes), "new_nodes": len(new_tree.nodes),
            "old_s": old_s, "new_s": new_s, "speedup": spd, "mad": mad, "changed": changed}


if __name__ == "__main__":
    args = sys.argv[1:]
    new_only = "--new-only" in args
    args = [a for a in args if a != "--new-only"]
    if not args:
        print(__doc__)
        print("usage: python scripts/parity_check_labeling.py <run_path> [iter=0] [tol_Jmol=1.0] [--new-only]")
        print("  --new-only : train just the fast new tree (skips the slow baseline); "
              "gives the speed number + a did-it-break check in ~one training.")
        sys.exit(0)
    run_path = args[0]
    it = int(args[1]) if len(args) > 1 else 0
    tol = float(args[2]) if len(args) > 2 else 1.0
    pairs, sampling = load_datums(run_path, it)
    print(f"[parity] loaded {len(pairs)} pairs + {len(sampling)} sampling datums "
          f"from {run_path} (iter {it})")
    if new_only:
        run_new_only(pairs, sampling, node_fract_training=0.7)
    else:
        run(pairs, sampling, tol=tol, node_fract_training=0.7)
