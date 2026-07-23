"""
Reproduce the flat-interaction collapse on COVERAGE configs (not training pairs).

eval_saved_tree.py showed the trained tree predicts a healthy SPREAD on its training
pairs. The flat [OH][Pt] plot is built from coverage configs, whose predicted energies
were computed during the run by tree.evaluate on multi-adsorbate configs. This script
loads the same tree and the run's stored coverage configs (Ncoad_config_*.json) and
re-evaluates them NOW, reporting:
  - spread of predicted interaction across configs (flat here == reproduced the bug),
  - for the first few configs: how many interaction subgraphs the decomposition yields,
    and how many of those match >=1 tree node (gim=False) -- i.e. whether coverage-config
    subgraphs fail to descend past the root.

If configs come back FLAT with subs that match 0 nodes -> the coverage-config decomposition
produces labeled subgraphs that don't match the (also labeled) tree nodes, even though the
training-pair subgraphs did. If configs come back with SPREAD -> the stored plot values are
stale/from a different computation.

Usage: python scripts/eval_coverage_configs.py /home/jzador/Pt/hox_mc_7 0 [n_configs=40]
"""
import sys
import os
import json
import glob
import statistics


def main(run_path, iteration=0, ncfg=40):
    from pysidt.sidt import read_nodes, MultiEvalSubgraphIsomorphicDecisionTreeRegressor
    from molecule.molecule import Molecule
    from pynta.coveragedependence import adsorbate_interaction_decomposition

    iter_path = os.path.join(run_path, "Iterations", str(iteration))
    nodes = read_nodes(os.path.join(iter_path, "regressor.json"))
    tree = MultiEvalSubgraphIsomorphicDecisionTreeRegressor(
        [adsorbate_interaction_decomposition], nodes=nodes)
    non_root = [n.group for n in nodes.values() if n.parent is not None and n.group is not None]
    print(f"loaded tree: {len(nodes)} nodes ({len(non_root)} non-root groups)")

    # collect coverage configs from Ncoad_config_*.json
    configs = []
    for f in sorted(glob.glob(os.path.join(iter_path, "Ncoad_config_*.json"))):
        with open(f) as fh:
            d = json.load(fh)
        for v in d.values():
            for adj in (v if isinstance(v, list) else [v]):
                try:
                    configs.append(Molecule().from_adjacency_list(adj, check_consistency=False))
                except Exception:
                    pass
    print(f"loaded {len(configs)} coverage configs from Ncoad_config_*.json")
    if not configs:
        print("no coverage configs found -- check the iteration path / file names")
        return

    preds = []
    for m in configs[:ncfg]:
        try:
            try:
                p = float(tree.evaluate(m, estimate_uncertainty=False))
            except TypeError:
                p = float(tree.evaluate(m))
            preds.append(p)
        except Exception as e:
            print("  eval failed:", repr(e))

    if preds:
        ps = statistics.pstdev(preds) if len(preds) > 1 else 0.0
        print(f"predicted interaction J/mol over {len(preds)} configs: "
              f"min {min(preds):.1f}  mean {statistics.mean(preds):.1f}  max {max(preds):.1f}  std {ps:.1f}")
        print("=> FLAT (bug reproduced)" if ps < 1000 else "=> SPREAD (stored plot may be stale)")

    # decomposition / match detail on the first few configs
    print("\nper-config decomposition (n_subgraphs, n_subs matching >=1 node):")
    for m in configs[:8]:
        try:
            subs = adsorbate_interaction_decomposition(m)
        except Exception as e:
            print(f"  decomposition RAISED: {e!r}")
            continue
        matched = 0
        for s in subs:
            if any(s.is_subgraph_isomorphic(g, generate_initial_map=False, save_order=True)
                   for g in non_root):
                matched += 1
        print(f"  n_subs={len(subs):>3}   matching_nonroot={matched:>3}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        print("usage: python scripts/eval_coverage_configs.py <run_path> [iteration=0] [n_configs=40]")
        sys.exit(0)
    main(sys.argv[1],
         int(sys.argv[2]) if len(sys.argv) > 2 else 0,
         int(sys.argv[3]) if len(sys.argv) > 3 else 40)
