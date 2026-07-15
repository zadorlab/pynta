"""
Canonical Metropolis Monte Carlo sampler over coverage-dependence configurations.

This is the MC front-end that replaces the enumerative configuration listing
(``get_configurations`` + the per-coverage minimum scan in ``get_cov_energies*``) for
regimes where enumeration is impossible (higher coverage on moderate slabs). It works
purely on the 2D adsorbate-site graph (``admol``):

  * a configuration is a base structure (one valid central adsorbate/TS site arrangement)
    decorated with a fixed number ``Ncoad`` of copies of a single coadsorbate on stable sites;
  * the central species is given as a LIST of valid base arrangements (all stable adsorbate
    geometries / all valid saddles). The MC explores the joint space of (which central
    arrangement, which sites hold a coadsorbate), so the central species can "hop" between
    arrangements rather than each chain being stuck on one;
  * the MC state is ``(base_idx, set of occupied coadsorbate-site ranks)`` -- hashable, cheap to
    dedup, and independent of atom-identity bookkeeping. Site "ranks" are the canonical order of
    surface-site atoms, shared across all base arrangements (same slab/site graph);
  * energies/uncertainties come from the trained SIDT regressor exactly as in
    ``get_cov_energies`` (E = atom-centered + interaction, both in J/mol);
  * two move types, canonical in #coadsorbates: a coadsorbate relocation (move one coad
    occupied->empty, mostly nearby with an occasional global jump) and, with probability
    ``p_central``, a central hop (swap to another valid central arrangement, keeping the coads);
  * acceptance is Metropolis on the uncertainty-corrected energy ``E_corr = E - lambda*sigma``
    with a within-chain annealed temperature ``T`` and exploration weight ``lambda``.

The visited (validly-evaluated) configurations are returned as a pool of
``(admol, E, sigma, trace)`` tuples, which the existing ``get_configs_for_calculation``
selection consumes unchanged. The 2D->3D geometry generation and calculation back-end is
untouched and only runs on the downselected configs.

See the enumerative counterpart in ``coveragedependence.py``.
"""

import os
import json
import numpy as np
import random
import logging
from collections import deque
from ase.io import read, write
import ase.units

from pynta.coveragedependence import (add_ad_to_site, configuration_is_valid,
                                       get_atom_centered_correction, get_central_templates,
                                       remove_slab, mol_to_atoms,
                                       generate_allowed_structure_site_structures)

R = ase.units._k * ase.units._Nav        # gas constant, J/mol/K (Boltzmann * Avogadro)
EV_TO_JMOL = ase.units.mol / ase.units.J  # eV -> J/mol (Faraday); matches get_cov_energies


def _coad_stable_sites_set(coad_stable_sites):
    return set(tuple(x) if not isinstance(x, str) else x for x in coad_stable_sites)


def get_stable_empty_site_inds(base_admol, coad_stable_sites):
    """Indices (into base_admol.atoms) of stable coadsorbate sites that are EMPTY in the base
    structure, i.e. the sites a coadsorbate may occupy. Mirrors the empty-site detection in
    ``get_configurations``."""
    s = _coad_stable_sites_set(coad_stable_sites)
    return [i for i, a in enumerate(base_admol.atoms)
            if a.is_surface_site() and (a.site, a.morphology) in s
            and not any(not a2.is_surface_site() for a2 in a.bonds.keys())]


class CanonicalCoverageMC:
    """Canonical (fixed-Ncoad) Metropolis MC over coadsorbate placements, with central hopping.

    Args:
        base_admols: list of valid central arrangements (all stable adsorbate geometries / all
            valid saddles). A single Molecule is accepted and wrapped in a list. All must be built
            on the same slab/site ordering (checked).
        coad: the (slab-free) coadsorbate Molecule to place
        coad_stable_sites: stable coadsorbate (site, morphology) list
        tree_interaction_regressor: trained SIDT interaction regressor; queried with
            ``evaluate(m, trace=True, estimate_uncertainty=True) -> (E, sigma, trace)``
        tree_atom_regressor / coadmol_E_dict: atom-centered (1-body) energy source (one required)
        unstable_pairs: list of unstable-pair Groups; proposals matching one are rejected. Falsy
            disables stability filtering.
        is_ts: whether the base structures are transition states (affects validity splitting)
        local_radius: site-graph hop radius defining a "nearby" empty site for local moves
        seed: optional base RNG seed
    """

    def __init__(self, base_admols, coad, coad_stable_sites, tree_interaction_regressor,
                 tree_atom_regressor=None, coadmol_E_dict=None, unstable_pairs=None,
                 is_ts=False, local_radius=2, seed=None):
        if tree_atom_regressor is None and coadmol_E_dict is None:
            raise ValueError("provide either tree_atom_regressor or coadmol_E_dict")
        if not isinstance(base_admols, (list, tuple)):
            base_admols = [base_admols]
        if len(base_admols) == 0:
            raise ValueError("base_admols is empty")
        self.base_admols = list(base_admols)
        self.coad = coad
        self.tree_interaction_regressor = tree_interaction_regressor
        self.tree_atom_regressor = tree_atom_regressor
        self.coadmol_E_dict = coadmol_E_dict
        self.unstable_pairs = unstable_pairs
        self.is_ts = is_ts
        self.local_radius = local_radius
        self.seed = seed
        self._coad_set = _coad_stable_sites_set(coad_stable_sites)

        # canonical site ranks: the ordered surface-site atoms, assumed consistent across bases
        # (same slab/sites -> same site atoms in the same order, with adsorbate atoms appended)
        self._site_inds = [[i for i, a in enumerate(b.atoms) if a.is_surface_site()]
                           for b in self.base_admols]
        self.nsites = len(self._site_inds[0])
        for b, si in enumerate(self._site_inds):
            if len(si) != self.nsites:
                raise ValueError("base structures have inconsistent surface-site counts "
                                 "({} vs {} for base {})".format(len(si), self.nsites, b))
        self._site_sm = [(self.base_admols[0].atoms[self._site_inds[0][k]].site,
                          self.base_admols[0].atoms[self._site_inds[0][k]].morphology)
                         for k in range(self.nsites)]
        for b, si in enumerate(self.base_admols):
            for k in range(self.nsites):
                a = self.base_admols[b].atoms[self._site_inds[b][k]]
                if (a.site, a.morphology) != self._site_sm[k]:
                    raise ValueError("site ordering differs across base structures at rank "
                                     "{} (base {})".format(k, b))

        # per-base: ranks occupied by the central species, and stable-empty ranks for coads
        self._central_ranks = []
        self._stable_empty = []
        for b, base in enumerate(self.base_admols):
            occ_ranks = set()
            for k, idx in enumerate(self._site_inds[b]):
                if any(not a2.is_surface_site() for a2 in base.atoms[idx].bonds.keys()):
                    occ_ranks.add(k)
            self._central_ranks.append(occ_ranks)
            self._stable_empty.append(set(k for k in range(self.nsites)
                                          if self._site_sm[k] in self._coad_set
                                          and k not in occ_ranks))
        if all(len(se) == 0 for se in self._stable_empty):
            raise ValueError("no empty stable coadsorbate sites on any base structure")

        self._dist = self._build_site_distances()

    # ------------------------------------------------------------------ graph helpers
    def _build_site_distances(self):
        """BFS site-rank hop distances on the lattice site graph (base-independent; uses base 0).
        Returns {rank: {rank: hops}} from every stable-empty rank."""
        base = self.base_admols[0]
        si = self._site_inds[0]
        atom_to_idx = {base.atoms[i]: i for i in si}
        idx_to_rank = {si[k]: k for k in range(self.nsites)}
        nbrs = {k: [idx_to_rank[atom_to_idx[a2]] for a2 in base.atoms[si[k]].bonds.keys()
                    if a2.is_surface_site()] for k in range(self.nsites)}
        stable_ranks = set().union(*self._stable_empty)
        dist = {}
        for s in stable_ranks:
            d = {s: 0}
            q = deque([s])
            while q:
                u = q.popleft()
                for v in nbrs[u]:
                    if v not in d:
                        d[v] = d[u] + 1
                        q.append(v)
            dist[s] = d
        return dist

    def _empty(self, base_idx, occ):
        return [k for k in self._stable_empty[base_idx] if k not in occ]

    def _nearby_empty(self, base_idx, o_rank, occ):
        d = self._dist.get(o_rank, {})
        return [k for k in self._stable_empty[base_idx]
                if k not in occ and d.get(k, np.inf) <= self.local_radius]

    # ------------------------------------------------------------------ config build/eval
    def build_config(self, base_idx, occ):
        """Build the admol for (base_idx, occ) by adding a coadsorbate at each occupied rank of the
        chosen base. Site ranks map to base-atom indices via self._site_inds[base_idx]; those
        indices survive add_ad_to_site (coad atoms are appended). May raise on a bad placement;
        callers treat that as an invalid proposal."""
        m = self.base_admols[base_idx]
        si = self._site_inds[base_idx]
        for rank in occ:
            m = add_ad_to_site(m, self.coad, m.atoms[si[rank]])
        return m

    def energy(self, m):
        """(E, sigma, trace) for a config, identical decomposition to get_cov_energies (J/mol)."""
        Einteraction, sigma, tr = self.tree_interaction_regressor.evaluate(
            m, trace=True, estimate_uncertainty=True)
        if self.tree_atom_regressor is not None:
            E = self.tree_atom_regressor.evaluate(m) + Einteraction
        else:
            E = get_atom_centered_correction(m, self.coadmol_E_dict) * EV_TO_JMOL + Einteraction
        return E, sigma, tr

    def is_valid(self, m):
        if not self.unstable_pairs:
            return True
        return configuration_is_valid(m, self.base_admols[0], self.is_ts, self.unstable_pairs)

    # ------------------------------------------------------------------ initialization
    def _greedy_fill(self, base_idx, Ncoad, rng):
        """Greedily place Ncoad valid coadsorbates on base_idx in random site order.
        Returns (occ, admol) or (None, None) if it gets stuck."""
        occ = set()
        m = self.base_admols[base_idx]
        while len(occ) < Ncoad:
            cands = list(self._stable_empty[base_idx] - occ)
            rng.shuffle(cands)
            placed = False
            for e in cands:
                cand = occ | {e}
                try:
                    mc = self.build_config(base_idx, cand)
                except Exception:
                    continue
                if self.is_valid(mc):
                    occ, m, placed = cand, mc, True
                    break
            if not placed:
                return None, None
        return occ, m

    def _greedy_valid_state(self, Ncoad, rng, max_restarts=200):
        """Greedily build a valid (base_idx, occ) with len(occ)==Ncoad. Bases are tried in index
        order; since get_central_templates sorts arrangements by central energy (lowest first), this
        starts the chain from the central's lowest-energy geometry -- usually the best site even under
        coverage, which helps learning. The chain can still leave it via central hops (p_central).
        Escalates to a higher-energy arrangement only if the lower one cannot host Ncoad coadsorbates
        validly. Returns (base_idx, occ, admol)."""
        viable = [b for b in range(len(self.base_admols)) if len(self._stable_empty[b]) >= Ncoad]
        if not viable:
            raise RuntimeError("no base structure can host Ncoad={} coadsorbates "
                               "(max stable sites {})".format(Ncoad, max(len(se) for se in self._stable_empty)))
        per_base = max(10, max_restarts // len(viable))
        for b in viable:
            for _ in range(per_base):
                occ, m = self._greedy_fill(b, Ncoad, rng)
                if occ is not None:
                    return b, occ, m
        raise RuntimeError("could not build a valid starting configuration with Ncoad={}".format(Ncoad))

    # ------------------------------------------------------------------ the sampler
    def run(self, Ncoad, n_steps=10000, T=(2000.0, 300.0), lam=(1.0, 0.0),
            p_local=0.8, p_central=0.1, seed=None, record_trajectory=False):
        """Run canonical Metropolis MC at fixed Ncoad with within-chain annealing.

        Annealing happens WITHIN the chain (T cooled T_start->T_end, lam start->end); it is not
        cooled between AL rounds.

        Args:
            Ncoad: number of coadsorbates (fixed throughout the chain)
            n_steps: number of MC steps
            T: temperature in K; (T_start, T_end) tuple geometrically annealed, scalar held constant
            lam: weight in E_corr = E - lam*sigma; (start, end) tuple linearly annealed, scalar constant
            p_local: probability a coad relocation targets a nearby (within local_radius) empty site
            p_central: probability a step is a central hop (swap to another valid central
                arrangement) rather than a coad relocation; only active with >1 base structure
            seed: RNG seed for this chain (defaults to the sampler's seed)
            record_trajectory: if True, also return "trajectory": the ORDERED chain as a list of
                (base_idx, sorted occupied-site-rank tuple, E_cur) after each evaluated step
                (state repeats on a rejected move). Rebuild a frame with build_config(base_idx, occ)
                and mol_to_atoms for 3D visualization. Off by default (memory).

        Returns a dict with: Ncoad, best=(admol, E), pool=[(admol, E, sigma, trace), ...] (unique
        validly-evaluated configs), diagnostics (n_steps, n_accept, accept_frac, n_central,
        n_central_accept, n_unique, Emin), and trajectory (None unless record_trajectory).
        """
        rng = random.Random(self.seed if seed is None else seed)
        T0, T1 = (T, T) if np.isscalar(T) else T
        l0, l1 = (lam, lam) if np.isscalar(lam) else lam
        nbases = len(self.base_admols)

        b, occ, m_cur = self._greedy_valid_state(Ncoad, rng)
        E_cur, s_cur, tr_cur = self.energy(m_cur)

        visited = {}  # (base_idx, frozenset(occ)) -> (admol, E, sigma, trace)

        def record(b_, occ_, m, E, s, tr):
            k = (b_, frozenset(occ_))
            if k not in visited:
                visited[k] = (m, E, s, tr)

        record(b, occ, m_cur, E_cur, s_cur, tr_cur)

        traj = [(b, tuple(sorted(occ)), E_cur)] if record_trajectory else None

        n_accept = 0
        n_central = 0
        n_central_accept = 0
        for step in range(n_steps):
            f = step / max(n_steps - 1, 1)
            Tk = T0 * (T1 / T0) ** f if (T0 > 0 and T1 > 0) else (T0 + (T1 - T0) * f)
            lam_s = l0 + (l1 - l0) * f

            do_central = (nbases > 1 and Ncoad > 0 and rng.random() < p_central)
            if do_central:
                # central hop: swap to another valid central arrangement, keeping the coads
                b_new = rng.choice([x for x in range(nbases) if x != b])
                if not all(k in self._stable_empty[b_new] for k in occ):
                    continue  # coad placement collides with the new central arrangement
                occ_new = set(occ)
            else:
                # coad relocation (canonical: preserves #coads)
                o = rng.choice(tuple(occ))
                if rng.random() < p_local:
                    cands = self._nearby_empty(b, o, occ) or self._empty(b, occ)
                else:
                    cands = self._empty(b, occ)
                if not cands:
                    continue
                e = rng.choice(cands)
                b_new = b
                occ_new = (occ - {o}) | {e}

            try:
                m_new = self.build_config(b_new, occ_new)
            except Exception:
                continue
            if not self.is_valid(m_new):
                continue

            E_new, s_new, tr_new = self.energy(m_new)
            record(b_new, occ_new, m_new, E_new, s_new, tr_new)
            if do_central:
                n_central += 1

            dEcorr = (E_new - lam_s * s_new) - (E_cur - lam_s * s_cur)
            if dEcorr <= 0 or rng.random() < np.exp(-dEcorr / (R * Tk)):
                b, occ, m_cur, E_cur, s_cur, tr_cur = b_new, occ_new, m_new, E_new, s_new, tr_new
                n_accept += 1
                if do_central:
                    n_central_accept += 1

            if record_trajectory:
                traj.append((b, tuple(sorted(occ)), E_cur))

        pool = list(visited.values())
        best = min(pool, key=lambda v: v[1])
        return {"Ncoad": Ncoad, "best": (best[0], best[1]), "pool": pool,
                "n_steps": n_steps, "n_accept": n_accept,
                "accept_frac": n_accept / n_steps if n_steps else 0.0,
                "n_central": n_central, "n_central_accept": n_central_accept,
                "n_unique": len(pool), "Emin": best[1], "trajectory": traj}

    def scan_coverages(self, Ncoads=None, max_coadsorbates=None, n_jobs=1, **run_kwargs):
        """Run an MC chain at each coverage level (one sampler explores all base arrangements via
        central hops) and assemble per-coverage minima plus a combined candidate pool.

        Args:
            Ncoads: coverage levels (number of coadsorbates) to sample; defaults to 1..max hostable
            max_coadsorbates: cap on Ncoad
            n_jobs: joblib workers for the independent per-coverage chains
            **run_kwargs: forwarded to run() (n_steps, T, lam, p_local, p_central)

        Returns a dict with Ncoad_energy_dict / Ncoad_config_dict (per-coverage best energy and
        config adjacency list), a combined ``pool``, and per-Ncoad ``results``.
        """
        max_possible = max(len(se) for se in self._stable_empty)
        if Ncoads is None:
            top = max_possible
            if max_coadsorbates is not None:
                top = min(top, max_coadsorbates)
            Ncoads = list(range(1, top + 1))
        else:
            Ncoads = [N for N in Ncoads
                      if (max_coadsorbates is None or N <= max_coadsorbates) and N <= max_possible]

        seeds = [None if self.seed is None else self.seed + N for N in Ncoads]

        if n_jobs == 1:
            runs = [self.run(N, seed=s, **run_kwargs) for N, s in zip(Ncoads, seeds)]
        else:
            from joblib import Parallel, delayed
            runs = Parallel(n_jobs=n_jobs)(
                delayed(self.run)(N, seed=s, **run_kwargs) for N, s in zip(Ncoads, seeds))

        results = {}
        Ncoad_energy_dict = {}
        Ncoad_config_dict = {}
        pool = []
        for N, r in zip(Ncoads, runs):
            results[N] = r
            Ncoad_energy_dict[N] = r["best"][1]
            Ncoad_config_dict[N] = r["best"][0].to_adjacency_list()
            pool.extend(r["pool"])
            logging.info("coverage_mc: Ncoad=%d Emin=%.1f J/mol unique=%d accept=%.2f central=%d/%d",
                         N, r["Emin"], r["n_unique"], r["accept_frac"],
                         r["n_central_accept"], r["n_central"])
        return {"results": results, "Ncoad_energy_dict": Ncoad_energy_dict,
                "Ncoad_config_dict": Ncoad_config_dict, "pool": pool}


def mc_cov_energies_configs_concern(base_admols, coad, coad_stable_sites, tree_interaction_regressor,
                                    Nocc_isolated=None, concern_energy_tol=None, coadmol_E_dict=None,
                                    tree_atom_regressor=None, unstable_pairs=None, is_ts=False,
                                    max_coadsorbates=None, Ncoads=None, n_jobs=1, seed=None, local_radius=2,
                                    chain_max_frames=500, **mc_run_kwargs):
    """MC analogue of ``get_cov_energies_configs_concern_tree`` (coveragedependence.py).

    Drop-in replacement: returns ``(Ncoad_energy_dict, Ncoad_config_dict, configs_of_concern)`` in
    the SAME shapes, so the enumerative pieces downstream (CalculateConfigurationEnergiesTask
    serialization, SelectCalculationsTask, get_configs_for_calculation) work unchanged. Runs ONE
    canonical-MC sampler that explores all base central arrangements (via central hops) and pools
    visited configurations per coverage.

    Args:
        base_admols: list of valid central arrangements (all stable adsorbate geometries / all
            valid saddles) -- the central species can hop among these during sampling
        coad: the (slab-free) coadsorbate Molecule
        coad_stable_sites: stable coadsorbate (site, morphology) list
        tree_interaction_regressor: trained SIDT interaction regressor
        Nocc_isolated: accepted for API compatibility but unused -- coverage is the #coadsorbates
            placed (len(occ)), which is robust to central arrangements occupying different site counts
        concern_energy_tol: configs within this many J/mol of their coverage's MC minimum are flagged
            "of concern" (the candidate set for selection); None keeps all visited
        coadmol_E_dict / tree_atom_regressor: atom-centered (1-body) energy source (one required)
        unstable_pairs: unstable-pair Groups for validity filtering
        is_ts: whether the base structures are transition states
        max_coadsorbates: cap on coverage
        n_jobs: joblib workers for the per-coverage chains
        seed: base RNG seed
        local_radius: site-graph hop radius for local moves (hops, not Angstrom; default 2)
        **mc_run_kwargs: forwarded to run (n_steps, T, lam, p_local, p_central)

    Returns:
        Ncoad_energy_dict: {Ncoad: min energy [J/mol]}
        Ncoad_config_dict: {Ncoad: adjacency list of the min-energy config}
        configs_of_concern: {i: (admol, E, trace, sigma)} (same value order as the enumerative
            function, which CalculateConfigurationEnergiesTask serializes as [adjlist, E, trace, sigma])
        diagnostics: {Ncoad: {Emin, n_steps, n_accept, accept_frac, n_unique, n_central,
            n_central_accept}} -- per-coverage MC health, persisted for inspection
        chains: {Ncoad: [[base_idx, [occ ranks], E], ...]} -- the strided ordered MC walk per
            coverage (capped at chain_max_frames), for visualization via write_mc_chain_xyz
        (the enumerative function returns only the first three values)
    """
    mc = CanonicalCoverageMC(base_admols, coad, coad_stable_sites, tree_interaction_regressor,
                             tree_atom_regressor=tree_atom_regressor, coadmol_E_dict=coadmol_E_dict,
                             unstable_pairs=unstable_pairs, is_ts=is_ts, local_radius=local_radius,
                             seed=seed)
    res = mc.scan_coverages(Ncoads=Ncoads, max_coadsorbates=max_coadsorbates, n_jobs=n_jobs,
                            record_trajectory=True, **mc_run_kwargs)

    Ncoad_energy_dict = res["Ncoad_energy_dict"]
    Ncoad_config_dict = res["Ncoad_config_dict"]
    configs_of_concern = {}
    diagnostics = {}
    chains = {}
    i = 0
    for N, r in res["results"].items():
        Emin = Ncoad_energy_dict[N]
        diagnostics[N] = {"Emin": r["Emin"], "n_steps": r["n_steps"], "n_accept": r["n_accept"],
                          "accept_frac": r["accept_frac"], "n_unique": r["n_unique"],
                          "n_central": r["n_central"], "n_central_accept": r["n_central_accept"]}
        # strided, capped chain for visualization: [(base_idx, [occ ranks], E), ...]
        traj = r.get("trajectory") or []
        if chain_max_frames and len(traj) > chain_max_frames:
            stride = (len(traj) + chain_max_frames - 1) // chain_max_frames
            traj = traj[::stride]
        chains[N] = [[b, list(occ), E] for (b, occ, E) in traj]
        for (m, E, sigma, tr) in r["pool"]:
            if concern_energy_tol is None or Emin + concern_energy_tol > E:
                configs_of_concern[i] = (m, E, tr, sigma)
                i += 1

    return Ncoad_energy_dict, Ncoad_config_dict, configs_of_concern, diagnostics, chains


def build_admol_from_occ(base_admol, coad, occ_ranks):
    """Rebuild a config admol from a base arrangement + occupied site ranks (the chain's compact
    representation). Mirrors CanonicalCoverageMC.build_config but standalone, so rendering a saved
    chain needs no sampler/MC state."""
    site_inds = [i for i, a in enumerate(base_admol.atoms) if a.is_surface_site()]
    m = base_admol
    for rank in occ_ranks:
        m = add_ad_to_site(m, coad, m.atoms[site_inds[rank]])
    return m


def write_mc_chain_xyz(path, pynta_run_directory, central, coadname, Ncoad, sites, site_adjacency,
                       metal, facet, iteration=None, adsorbate_site_energy_cutoff=None,
                       out_xyz=None, slab=None, allowed_structure_site_structures=None):
    """Render a SAVED MC chain to a multi-frame .xyz (one frame per recorded step) for visualization.

    Needs only RUN INPUTS (the same things the covdep run / postprocess was given) -- no MC
    hyperparameters (T/lam/p_local/p_central/local_radius): the chain is already determined, this
    just turns each 2D config into 3D. Reads Iterations/<iter>/mc_chain_<central>_<coad>.json,
    rebuilds the central base arrangements (get_central_templates -> must match the run, so pass the
    run's adsorbate_site_energy_cutoff), and places coadsorbates via mol_to_atoms.

    Args:
        path: covdep run directory
        pynta_run_directory: the isolated pynta run (slab, Adsorbates/, TS*/)
        central: central species admol_name (e.g. "TS0_21" or an adsorbate SMILES)
        coadname: coadsorbate name
        Ncoad: which coverage's chain to render
        sites, site_adjacency, metal, facet: run inputs (facet is the ACAT surface type the run
            used, e.g. "fcc111" -- same variable the covdep postprocess already defines)
        iteration: which Iterations/<i>/; defaults to the latest present
        adsorbate_site_energy_cutoff: must match the run (defines the base-arrangement set/order)
        out_xyz: output path; defaults to mc_chain_<central>_<coad>_N<Ncoad>_iter<i>.xyz
        slab, allowed_structure_site_structures: optional precomputed inputs

    Returns (out_xyz, frames, energies): the written path (auto-named with the resolved iteration
    unless out_xyz was given), the list of ase.Atoms, and the per-frame energies [J/mol].
    """
    if iteration is None:
        iters = [int(d) for d in os.listdir(os.path.join(path, "Iterations")) if d.isdigit()]
        iteration = max(iters)
    chain_file = os.path.join(path, "Iterations", str(iteration), "mc_chain_" + central + "_" + coadname + ".json")
    with open(chain_file) as f:
        chains = json.load(f)
    if str(Ncoad) not in chains:
        raise KeyError("no chain for Ncoad={} in {} (have {})".format(Ncoad, chain_file, list(chains.keys())))
    chain = chains[str(Ncoad)]

    if slab is None:
        slab = read(os.path.join(pynta_run_directory, "slab.xyz"))
    nslab = len(slab)
    if allowed_structure_site_structures is None:
        allowed_structure_site_structures = generate_allowed_structure_site_structures(
            os.path.join(pynta_run_directory, "Adsorbates"), sites, site_adjacency, nslab, max_dist=np.inf)

    is_ts = central.startswith("TS")
    templates = get_central_templates(central, is_ts, pynta_run_directory, metal, facet, sites,
                                      site_adjacency, nslab,
                                      allowed_structure_site_structures=allowed_structure_site_structures,
                                      energy_cutoff=adsorbate_site_energy_cutoff)
    coad_templates = get_central_templates(coadname, False, pynta_run_directory, metal, facet, sites,
                                           site_adjacency, nslab,
                                           allowed_structure_site_structures=allowed_structure_site_structures,
                                           energy_cutoff=None)
    coad_simple = remove_slab(coad_templates[0][1])

    frames = []
    energies = []
    for entry in chain:
        b, occ, E = entry[0], entry[1], entry[2]
        base_admol = templates[b][1]
        nbase = len(base_admol.atoms)
        admol = build_admol_from_occ(base_admol, coad_simple, occ)
        # central atoms are the base's atoms (build_admol_from_occ appends coadsorbates), so we can
        # tell mol_to_atoms directly and skip the subgraph isomorphism
        central_atoms = [a for a in admol.atoms[:nbase] if not a.is_surface_site()]
        # random_offset=False -> the same config renders identically, so converged (rejected) runs
        # show as still frames and only real site-hops change the picture
        atoms = mol_to_atoms(admol, slab, sites, metal,
                             partial_atoms=templates[b][0], atoms_in_partial=central_atoms,
                             random_offset=False)
        frames.append(atoms)
        energies.append(E)

    if out_xyz is None:
        out_xyz = "mc_chain_{}_{}_N{}_iter{}.xyz".format(central, coadname, Ncoad, iteration)
    write(out_xyz, frames)
    return out_xyz, frames, energies


def plot_mc_frame_interaction_graph(path, pynta_run_directory, central, coadname, Ncoad, frame,
                                    sites, site_adjacency, metal, facet, iteration=None,
                                    adsorbate_site_energy_cutoff=None, slab=None,
                                    allowed_structure_site_structures=None, out_png=None, **plot_kwargs):
    """Plot the SIDT interaction-energy decomposition of ONE saved MC chain frame (no 3D round-trip).

    Thin wrapper around pynta.utils.plot_interaction_graph: takes the same run inputs as
    write_mc_chain_xyz plus ``frame`` (index into this Ncoad's chain). Loads the chain entry
    (base_idx, occ), rebuilds the 2D config admol via build_admol_from_occ (so close adsorbates are
    not merged), loads the iteration's interaction regressor, and plots. Returns (fig, ax).
    """
    from pynta.utils import plot_interaction_graph
    from pysidt.sidt import read_nodes, MultiEvalSubgraphIsomorphicDecisionTreeRegressor
    from pynta.coveragedependence import adsorbate_interaction_decomposition

    if iteration is None:
        iters = [int(d) for d in os.listdir(os.path.join(path, "Iterations")) if d.isdigit()]
        iteration = max(iters)
    chain_file = os.path.join(path, "Iterations", str(iteration), "mc_chain_" + central + "_" + coadname + ".json")
    with open(chain_file) as f:
        chains = json.load(f)
    if str(Ncoad) not in chains:
        raise KeyError("no chain for Ncoad={} in {} (have {})".format(Ncoad, chain_file, list(chains.keys())))
    entry = chains[str(Ncoad)][frame]
    b, occ, E = entry[0], entry[1], entry[2]

    if slab is None:
        slab = read(os.path.join(pynta_run_directory, "slab.xyz"))
    nslab = len(slab)
    if allowed_structure_site_structures is None:
        allowed_structure_site_structures = generate_allowed_structure_site_structures(
            os.path.join(pynta_run_directory, "Adsorbates"), sites, site_adjacency, nslab, max_dist=np.inf)

    is_ts = central.startswith("TS")
    templates = get_central_templates(central, is_ts, pynta_run_directory, metal, facet, sites,
                                      site_adjacency, nslab,
                                      allowed_structure_site_structures=allowed_structure_site_structures,
                                      energy_cutoff=adsorbate_site_energy_cutoff)
    coad_templates = get_central_templates(coadname, False, pynta_run_directory, metal, facet, sites,
                                           site_adjacency, nslab,
                                           allowed_structure_site_structures=allowed_structure_site_structures,
                                           energy_cutoff=None)
    coad_simple = remove_slab(coad_templates[0][1])
    admol = build_admol_from_occ(templates[b][1], coad_simple, occ)

    nodes = read_nodes(os.path.join(path, "Iterations", str(iteration), "regressor.json"))
    tree = MultiEvalSubgraphIsomorphicDecisionTreeRegressor([adsorbate_interaction_decomposition], nodes=nodes)

    title = plot_kwargs.pop("title", "MC {} +{}x{} iter{} frame{} (E={:.3f} eV)".format(
        central, Ncoad, coadname, iteration, frame, E / EV_TO_JMOL))
    return plot_interaction_graph(admol, tree, sites, site_adjacency, slab,
                                  out_png=out_png, title=title, **plot_kwargs)
