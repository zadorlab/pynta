#!/usr/bin/env python
"""Double-diffusion (concerted two-species hop) saddle search, v0 hack.

Hypothesis: pairing two known single-diffusion saddles (TS_A + TS_B placed near each other)
seeds guesses for the CONCERTED saddle where both species change sites through one first-order
saddle. A paired-saddle start is ~index-2; Sella (order=1) will kill one uphill direction, so
each result is classified afterwards by its imaginary-mode composition:
  concerted    : exactly 1 imaginary mode with significant amplitude on BOTH species
  single+spect : imaginary mode lives on one species (B slid into a well) -- still a useful
                 "TS_A in the presence of B" datum, but not the double diffusion
  other        : wrong number of imaginary modes / not converged

Usage (run in the isolated pynta run directory, pynta_env):
  python double_diffusion.py generate            # writes doublediff/<A>-<B>/<k>/init.xyz + meta.json
  python double_diffusion.py opt <placement_dir> # one Sella order-1 opt + vibrations (sbatch this)
  python double_diffusion.py analyze             # collate results.json -> table + CSV

generate prints an sbatch loop template. No pynta-workflow plumbing is touched; this reuses the
run's sites/neighbors files and finished TS guess geometries only (nothing is recomputed).
"""
import os
import sys
import json
import glob
import numpy as np
from ase.io import read, write
from ase.geometry import get_distances

# ---------------------------------------------------------------- configuration
RUN_DIR   = "."                                   # the isolated pynta run (TS*/ dirs, slab.xyz)
SITES     = "sites.json"                          # the run's saved site list
MODEL     = "/home/jzador/macemodel/mace-mh-1.model"
MACE_HEAD = "omat_pbe"

# diffusion TS dirs -> species label; guess index None = auto-pick lowest-energy guess with files
TS_SPECIES = {"TS1": "NH3", "TS2": "O", "TS3": "NH2", "TS4": "OH"}
TS_GUESSES = {"TS1": "18", "TS2": None, "TS3": None, "TS4": None}

MAX_SEP       = 6.0    # A: max separation (min A-atom..B-atom distance) of a placement
MIN_CLEARANCE = 1.6    # A: reject placements with any A..B atom pair closer than this
N_PER_COMBO   = 8      # keep the N closest placements per species pair
FMAX          = 0.05
OPT_STEPS     = 200
IMAG_THRESH   = 25.0   # cm^-1: |Im(freq)| above this counts as an imaginary mode
CONCERTED_MIN = 0.25   # min mode-0 amplitude fraction on EACH species to call it concerted
OUT_ROOT      = "doublediff"
# -------------------------------------------------------------------------------


def load_sites():
    with open(os.path.join(RUN_DIR, SITES)) as f:
        sites = json.load(f)
    for s in sites:
        s["position"] = np.array(s["position"])
    return sites


def pick_guess(tsdir):
    """Configured guess, or the lowest-energy guess dir that has opt.xyz + vib.0.traj."""
    g = TS_GUESSES.get(tsdir)
    if g is not None:
        return str(g)
    best = None
    for d in sorted(os.listdir(os.path.join(RUN_DIR, tsdir))):
        p = os.path.join(RUN_DIR, tsdir, d)
        if not (os.path.isdir(p) and os.path.exists(os.path.join(p, "opt.xyz"))
                and os.path.exists(os.path.join(p, "vib.0.traj"))):
            continue
        try:
            E = read(os.path.join(p, "opt.xyz")).get_potential_energy()
        except Exception:
            continue
        if best is None or E < best[1]:
            best = (d, E)
    if best is None:
        raise RuntimeError(f"{tsdir}: no guess dir with opt.xyz + vib.0.traj and an energy; "
                           f"set TS_GUESSES['{tsdir}'] explicitly")
    return best[0]


def load_saddle(tsdir, guess, nslab):
    atoms = read(os.path.join(RUN_DIR, tsdir, guess, "opt.xyz"))
    try:
        E = atoms.get_potential_energy()
    except Exception:
        E = None
    return atoms, E


def anchor_site(atoms, nslab, sites):
    """(site_dict, binding_atom_index) for the fragment: heavy (non-H if possible) adsorbate atom
    closest to the surface, anchored to its nearest site (xy distance)."""
    syms = atoms.get_chemical_symbols()
    ads = list(range(nslab, len(atoms)))
    heavy = [i for i in ads if syms[i] != "H"] or ads
    b = min(heavy, key=lambda i: atoms.positions[i, 2])
    pos = atoms.positions[b]
    best = min(sites, key=lambda s: np.hypot(*(s["position"][:2] - pos[:2])))
    return best, b


def mic_min_dist(p1, p2, cell):
    v, d = get_distances(p1, p2, cell=cell, pbc=[True, True, False])
    return float(d.min())


def generate():
    slab = read(os.path.join(RUN_DIR, "slab.xyz"))
    nslab = len(slab)
    sites = load_sites()

    saddles = {}
    for tsdir, name in TS_SPECIES.items():
        guess = pick_guess(tsdir)
        atoms, E = load_saddle(tsdir, guess, nslab)
        site, b = anchor_site(atoms, nslab, sites)
        saddles[tsdir] = dict(name=name, guess=guess, atoms=atoms, E=E,
                              anchor=site, anchor_pos=atoms.positions[b].copy())
        print(f"{tsdir} ({name}): guess {guess}, E={E}, anchor {site['site']}/{site['morphology']}")

    combos = []
    keys = list(TS_SPECIES.keys())
    for i, a in enumerate(keys):
        for b in keys[i:]:
            combos.append((a, b))

    all_dirs = []
    for A, B in combos:
        sa, sb = saddles[A], saddles[B]
        atomsA, atomsB = sa["atoms"], sb["atoms"]
        fragB = atomsB[nslab:]
        cell = atomsA.cell
        aid = (sb["anchor"]["site"], sb["anchor"]["morphology"])
        # candidate translations: B's anchor moved onto every congruent site (same type/morph)
        cands = []
        for s in sites:
            if (s["site"], s["morphology"]) != aid:
                continue
            t = np.array([s["position"][0] - sb["anchor_pos"][0],
                          s["position"][1] - sb["anchor_pos"][1], 0.0])
            posB = fragB.positions + t
            sep = mic_min_dist(atomsA.positions[nslab:], posB, cell)
            if sep < MIN_CLEARANCE or sep > MAX_SEP:
                continue
            cands.append((sep, t, s))
        cands.sort(key=lambda c: c[0])
        # dedup near-identical translations (MIC-equivalent placements)
        seen = []
        kept = []
        for sep, t, s in cands:
            key = tuple(np.round(t[:2], 1))
            if key in seen:
                continue
            seen.append(key)
            kept.append((sep, t, s))
            if len(kept) >= N_PER_COMBO:
                break

        label = f"{sa['name']}-{sb['name']}"
        for k, (sep, t, s) in enumerate(kept):
            d = os.path.join(RUN_DIR, OUT_ROOT, label, str(k))
            os.makedirs(d, exist_ok=True)
            combined = atomsA.copy()          # A's saddle incl. its (relaxed) slab
            fB = fragB.copy()
            fB.translate(t)
            combined += fB
            combined.wrap()
            write(os.path.join(d, "init.xyz"), combined)
            meta = dict(A=dict(tsdir=A, name=sa["name"], guess=sa["guess"], E=sa["E"]),
                        B=dict(tsdir=B, name=sb["name"], guess=sb["guess"], E=sb["E"]),
                        nslab=nslab, nA=len(atomsA) - nslab, nB=len(fragB),
                        target_site=dict(site=s["site"], morphology=s["morphology"],
                                         position=list(map(float, s["position"]))),
                        separation=float(sep))
            with open(os.path.join(d, "meta.json"), "w") as f:
                json.dump(meta, f, indent=1)
            all_dirs.append(d)
        print(f"{label}: {len(kept)} placements (of {len(cands)} candidates)")

    listfile = os.path.join(RUN_DIR, OUT_ROOT, "placements.txt")
    with open(listfile, "w") as f:
        f.write("\n".join(all_dirs) + "\n")
    print(f"\n{len(all_dirs)} placements -> {listfile}")
    print("\nsbatch template (1 core each, adjust partition/account):\n")
    print("  while read d; do")
    print("    sbatch -n1 --mem-per-cpu=8000 -t 24:00:00 \\")
    print(f"      --wrap \"python {os.path.abspath(__file__)} opt $d\"")
    print(f"  done < {listfile}")


def make_calc():
    from mace.calculators import MACECalculator
    from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
    from ase.calculators.mixing import SumCalculator
    # kwargs mirror the covdep run's software_kwargs exactly (model_path singular -- the spelling
    # the installed mace version demonstrably accepts)
    mace = MACECalculator(model_path=MODEL, head=MACE_HEAD, device="cpu", default_dtype="float64")
    d3 = TorchDFTD3Calculator(device="cpu", damping="bj", xc="pbe")
    return SumCalculator([mace, d3])


def opt(d):
    from ase.constraints import FixAtoms
    from ase.vibrations import Vibrations
    from sella import Sella

    with open(os.path.join(d, "meta.json")) as f:
        meta = json.load(f)
    nslab, nA, nB = meta["nslab"], meta["nA"], meta["nB"]
    atoms = read(os.path.join(d, "init.xyz"))
    atoms.pbc = [True, True, False]
    atoms.set_constraint(FixAtoms(indices=list(range(nslab))))
    atoms.calc = make_calc()

    res = dict(converged=False)
    try:
        dyn = Sella(atoms, order=1, trajectory=os.path.join(d, "opt.traj"))
        res["converged"] = bool(dyn.run(fmax=FMAX, steps=OPT_STEPS))
        write(os.path.join(d, "opt.xyz"), atoms)
        res["energy"] = float(atoms.get_potential_energy())

        free = list(range(nslab, len(atoms)))
        vib = Vibrations(atoms, indices=free, name=os.path.join(d, "vib"))
        vib.run()
        freqs = vib.get_frequencies()          # complex, cm^-1
        imag = [f for f in freqs if abs(f.imag) > IMAG_THRESH]
        res["freqs_imag"] = [float(f.imag) for f in imag]
        res["n_imag"] = len(imag)
        mode0 = vib.get_mode(0)                # (natoms, 3), zeros on frozen atoms
        aA = float(np.linalg.norm(mode0[nslab:nslab + nA]) ** 2)
        aB = float(np.linalg.norm(mode0[nslab + nA:nslab + nA + nB]) ** 2)
        tot = aA + aB
        res["frac_A"] = aA / tot if tot else 0.0
        res["frac_B"] = aB / tot if tot else 0.0
    except Exception as e:
        res["error"] = f"{type(e).__name__}: {e}"
    with open(os.path.join(d, "results.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


def analyze():
    try:
        Eslab = read(os.path.join(RUN_DIR, "slab.xyz")).get_potential_energy()
    except Exception:
        Eslab = None
    rows = []
    for d in sorted(glob.glob(os.path.join(RUN_DIR, OUT_ROOT, "*", "*"))):
        rp = os.path.join(d, "results.json")
        mp = os.path.join(d, "meta.json")
        if not (os.path.exists(rp) and os.path.exists(mp)):
            continue
        res = json.load(open(rp))
        meta = json.load(open(mp))
        cls = "error" if "error" in res else "unconverged" if not res.get("converged") else \
              "not_saddle" if res.get("n_imag") != 1 else \
              "concerted" if min(res.get("frac_A", 0), res.get("frac_B", 0)) >= CONCERTED_MIN else \
              "single+spectator"
        # TS-TS interaction: E(A+B saddle) - E(A saddle) - E(B saddle) + E(slab)
        inter = ""
        if (Eslab is not None and res.get("energy") is not None
                and meta["A"]["E"] is not None and meta["B"]["E"] is not None):
            inter = round(res["energy"] - meta["A"]["E"] - meta["B"]["E"] + Eslab, 3)
        rows.append([os.path.relpath(d, RUN_DIR), meta["A"]["name"], meta["B"]["name"],
                     round(meta["separation"], 2), cls, res.get("n_imag", ""),
                     round(res.get("frac_A", 0), 2), round(res.get("frac_B", 0), 2), inter])
    hdr = ["dir", "A", "B", "sep_A", "class", "n_imag", "frac_A", "frac_B", "Einter_eV"]
    wid = [max(len(str(r[i])) for r in rows + [hdr]) for i in range(len(hdr))] if rows else []
    for r in [hdr] + rows:
        print("  ".join(str(x).ljust(w) for x, w in zip(r, wid)))
    out = os.path.join(RUN_DIR, OUT_ROOT, "double_diffusion_results.csv")
    with open(out, "w") as f:
        f.write(",".join(hdr) + "\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")
    print(f"\nwrote {out}  ({sum(1 for r in rows if r[4] == 'concerted')} concerted "
          f"of {len(rows)} finished placements)")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "generate":
        generate()
    elif cmd == "opt" and len(sys.argv) > 2:
        opt(sys.argv[2])
    elif cmd == "analyze":
        analyze()
    else:
        print(__doc__)
