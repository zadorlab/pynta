"""
FAST Check A: load an already-trained covdep regressor (Iterations/<iter>/regressor.json)
and evaluate it on the training pairs -- NO retraining. Tells us whether the trained tree
predicts a real spread or a flat constant, i.e. whether the observed flat interaction
comes from the trained tree itself (TRAINING/FIT problem) or only appears later
(EVAL/POSTPROCESSING problem).

Loads the tree exactly the way postprocessing does:
    nodes = read_nodes(regressor.json)
    tree  = MultiEvalSubgraphIsomorphicDecisionTreeRegressor([adsorbate_interaction_decomposition], nodes=nodes)
    tree.evaluate(mol)

Usage (cluster, pynta_env):
    python scripts/eval_saved_tree.py /home/jzador/Pt/hox_mc_5 0
    (arg1 = run path, arg2 = iteration whose regressor.json to load, default 0)
"""
import sys
import os
import json
import statistics


def main(run_path, iteration=0):
    from pysidt.sidt import read_nodes, MultiEvalSubgraphIsomorphicDecisionTreeRegressor, Datum
    from molecule.molecule import Molecule
    from pynta.coveragedependence import adsorbate_interaction_decomposition

    reg = os.path.join(run_path, "Iterations", str(iteration), "regressor.json")
    nodes = read_nodes(reg)
    tree = MultiEvalSubgraphIsomorphicDecisionTreeRegressor(
        [adsorbate_interaction_decomposition], nodes=nodes)
    print(f"loaded {reg}: {len(nodes)} nodes")

    with open(os.path.join(run_path, "pairs_datums.json")) as f:
        pairs = [Datum(mol=Molecule().from_adjacency_list(d["mol"], check_consistency=False),
                       value=d["value"]) for d in json.load(f)]

    preds, acts = [], []
    fail = 0
    for d in pairs:
        try:
            try:
                p = float(tree.evaluate(d.mol, estimate_uncertainty=False))
            except TypeError:
                p = float(tree.evaluate(d.mol))
        except Exception:
            fail += 1
            continue
        preds.append(p)
        acts.append(float(d.value))

    print(f"evaluated {len(preds)}/{len(pairs)} training pairs ({fail} failed)")
    if not preds:
        print("no predictions -- evaluate failed on everything")
        return

    ps = statistics.pstdev(preds) if len(preds) > 1 else 0.0
    as_ = statistics.pstdev(acts) if len(acts) > 1 else 0.0
    distinct = len({round(p, 1) for p in preds})
    print(f"PRED   J/mol: min {min(preds):10.1f}  mean {statistics.mean(preds):10.1f}  "
          f"max {max(preds):10.1f}  std {ps:10.1f}")
    print(f"ACTUAL J/mol: min {min(acts):10.1f}  mean {statistics.mean(acts):10.1f}  "
          f"max {max(acts):10.1f}  std {as_:10.1f}")
    print(f"distinct predicted values (rounded 0.1 J/mol): {distinct}")
    print("=" * 60)
    if ps < 1000.0 or distinct <= 5:
        print("VERDICT: FLAT -- the trained tree predicts ~constant on its OWN training")
        print("  pairs. The collapse is in TRAINING / the FIT (look at the Lasso), not eval.")
    else:
        print("VERDICT: SPREAD -- the trained tree differentiates its training pairs.")
        print("  The flat [OH][Pt] plot is therefore introduced in EVAL / POSTPROCESSING,")
        print("  not in the trained tree. Hunt the postprocessing evaluation path next.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        print("usage: python scripts/eval_saved_tree.py <run_path> [iteration=0]")
        sys.exit(0)
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 0)
