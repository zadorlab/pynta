"""
Diagnostic: WHY did the distance-labeled interaction subgraphs stop matching the
tree's child node groups, collapsing the model to its root constant?

For a sample of datums it decomposes each with adsorbate_interaction_decomposition
(these labeled sub-structures are exactly what descend the tree), and for each
sub-structure asks whether it subgraph-matches ANY of the get_adsorbed_atom_pairs
child groups under three regimes:

  (1) labeled group,   generate_initial_map=False   <- what runs now (suspected broken)
  (2) labeled group,   generate_initial_map=True     <- semantically complete engine
  (3) unlabeled group, generate_initial_map=False    <- original pre-change behavior
                                                        (numbered *N labels stripped
                                                         from BOTH sides; bare "*" kept)

Verdict logic:
  (3) high, (1) low, (2) HIGH   -> label_isomorphism / generate_initial_map=False ENGINE
                                    bug: the labels are fine, the fast matcher returns
                                    false negatives on multi-duplicate-label patterns.
                                    Fix = generate_initial_map=True on these calls (slow)
                                    or report the molecule branch upstream.
  (3) high, (1) low, (2) LOW    -> the LABELING itself is over-restrictive; idea is flawed
                                    as implemented, drop it.
  (1) ~= (3), both high         -> matching is NOT the problem; look at the fit/Lasso.

Usage (cluster, pynta_env):
    python scripts/diagnose_labeling_match.py /home/jzador/Pt/hox_mc_5 [n_datums=20]
"""
import sys
import os
import json
import re

_NUM = re.compile(r"^\*\d+$")  # "*1","*2",... but not bare "*"


def load_pairs(run_path):
    from molecule.molecule import Molecule
    from pysidt.sidt import Datum
    with open(os.path.join(run_path, "pairs_datums.json")) as f:
        return [Datum(mol=Molecule().from_adjacency_list(d["mol"], check_consistency=False),
                      value=d["value"]) for d in json.load(f)]


def strip_numbered(struct):
    """Copy with numbered *N labels cleared (keep bare '*'), emulating pre-change labeling."""
    s = struct.copy(deep=True)
    for a in s.atoms:
        if getattr(a, "label", "") and _NUM.match(a.label):
            a.label = ""
    return s


def matches_any(sub, groups, gim):
    c = 0
    for g in groups:
        try:
            if sub.is_subgraph_isomorphic(g, generate_initial_map=gim, save_order=True):
                c += 1
        except Exception:
            pass
    return c


def main(run_path, n=20):
    import pynta.coveragedependence as cd

    datums = load_pairs(run_path)[:n]
    groups = cd.get_adsorbed_atom_pairs(length=7, r_bonds=[1, 2, 3, 0.05, 0])
    groups_stripped = [strip_numbered(g) for g in groups]

    n_sub = 0
    hit1 = hit2 = hit3 = 0
    examples = []
    for d in datums:
        try:
            subs = cd.adsorbate_interaction_decomposition(d.mol)
        except Exception as e:
            print(f"  (decomposition failed on a datum: {e!r})")
            continue
        for sub in subs:
            n_sub += 1
            m1 = matches_any(sub, groups, False)
            m2 = matches_any(sub, groups, True)
            sub_stripped = strip_numbered(sub)
            m3 = matches_any(sub_stripped, groups_stripped, False)
            hit1 += (m1 > 0)
            hit2 += (m2 > 0)
            hit3 += (m3 > 0)
            if len(examples) < 8:
                examples.append((m1, m2, m3))

    print("=" * 60)
    print(f"datums sampled       : {len(datums)}")
    print(f"interaction subgraphs: {n_sub}")
    if n_sub == 0:
        print("no sub-structures produced -- try more datums or check the run folder")
        return
    print(f"(1) labeled + gim=False (current) : {hit1:>4}/{n_sub} match >=1 child group")
    print(f"(2) labeled + gim=True  (complete): {hit2:>4}/{n_sub}")
    print(f"(3) unlabeled (original behavior) : {hit3:>4}/{n_sub}")
    print("per-substructure #child-groups matched (first few): (gimF, gimT, unlabeled)")
    for e in examples:
        print("   ", e)

    print("=" * 60)
    f1, f2, f3 = hit1 / n_sub, hit2 / n_sub, hit3 / n_sub
    if f3 > 0.5 and f1 < 0.2 and f2 > 0.5:
        print("VERDICT: ENGINE BUG in label_isomorphism / generate_initial_map=False.")
        print("  Labels are semantically fine (gim=True matches), the fast matcher gives")
        print("  false negatives. Fix = gim=True on these calls, or report molecule upstream.")
    elif f3 > 0.5 and f1 < 0.2 and f2 < 0.5:
        print("VERDICT: the LABELING is over-restrictive (even gim=True fails). Idea is")
        print("  flawed as implemented -- revert the labeling.")
    elif abs(f1 - f3) < 0.15:
        print("VERDICT: matching is NOT the problem (labeled ~ unlabeled). Look at the")
        print("  fit/Lasso for the constant-prediction collapse instead.")
    else:
        print("VERDICT: mixed -- inspect the per-substructure numbers above.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        print("usage: python scripts/diagnose_labeling_match.py <run_path> [n_datums=20]")
        sys.exit(0)
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 20)
