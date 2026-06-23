"""
Canonical Metropolis Monte Carlo sampler over coverage-dependence configurations.

This is the MC front-end that replaces the enumerative configuration listing
(``get_configurations`` + the per-coverage minimum scan in ``get_cov_energies*``) for
regimes where enumeration is impossible (higher coverage on moderate slabs). It works
purely on the 2D adsorbate-site graph (``admol``):

  * a configuration is the base structure (central adsorbate/TS) decorated with a fixed
    number ``Ncoad`` of copies of a single coadsorbate on stable sites;
  * the MC state is just the SET of base-graph site indices currently holding a coad, so
    it is hashable (cheap dedup/memoization) and independent of atom-identity bookkeeping;
  * energies/uncertainties come from the trained SIDT regressor exactly as in
    ``get_cov_energies`` (E = atom-centered + interaction, both in J/mol);
  * moves are canonical relocations (move one coad occupied->empty, preserving Ncoad),
    mostly to a nearby empty site with an occasional global jump;
  * acceptance is Metropolis on the uncertainty-corrected energy ``E_corr = E - lambda*sigma``
    with an annealed temperature ``T`` and exploration weight ``lambda``.

The visited (validly-evaluated) configurations are returned as a pool of
``(admol, E, sigma, trace)`` tuples, which the existing ``get_configs_for_calculation``
selection consumes unchanged. The 2D->3D geometry generation and calculation back-end is
untouched and only runs on the downselected configs.

See the enumerative counterpart in ``coveragedependence.py``.
"""

import numpy as np
import random
import logging
from collections import deque

from pynta.coveragedependence import (add_ad_to_site, configuration_is_valid,
                                       get_atom_centered_correction)

R = 8.314  # J/mol/K, matches the Boltzmann weighting used elsewhere in covdep
EV_TO_JMOL = 96.48530749925793 * 1000.0  # eV -> J/mol, matches get_cov_energies


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
    """Canonical (fixed-Ncoad) Metropolis MC over coadsorbate placements on a fixed base admol.

    Args:
        base_admol: the base 2D structure (central adsorbate/TS on the slab, no coadsorbates)
        coad: the coadsorbate Molecule to place (single species per sampler)
        coad_stable_sites: iterable of (site, morphology) tuples (or strings) that are stable
            coadsorbate sites
        tree_interaction_regressor: trained SIDT interaction regressor; queried with
            ``evaluate(m, trace=True, estimate_uncertainty=True) -> (E, sigma, trace)``
        tree_atom_regressor: optional atom-centered (1-body) regressor; if given E uses it
        coadmol_E_dict: optional atom-centered correction dict; used if tree_atom_regressor is None
        unstable_pairs: list of unstable-pair Groups (from the pair calculations); proposals
            matching one are rejected. If falsy, no stability filtering is applied.
        is_ts: whether the base structure is a transition state (affects validity splitting)
        local_radius: site-graph hop radius defining a "nearby" empty site for local moves
        seed: optional RNG seed for reproducibility
    """

    def __init__(self, base_admol, coad, coad_stable_sites, tree_interaction_regressor,
                 tree_atom_regressor=None, coadmol_E_dict=None, unstable_pairs=None,
                 is_ts=False, local_radius=2, seed=None):
        if tree_atom_regressor is None and coadmol_E_dict is None:
            raise ValueError("provide either tree_atom_regressor or coadmol_E_dict")
        self.base_admol = base_admol
        self.coad = coad
        self.coad_stable_sites = coad_stable_sites
        self.tree_interaction_regressor = tree_interaction_regressor
        self.tree_atom_regressor = tree_atom_regressor
        self.coadmol_E_dict = coadmol_E_dict
        self.unstable_pairs = unstable_pairs
        self.is_ts = is_ts
        self.local_radius = local_radius
        if seed is not None:
            random.seed(seed)

        self.stable_inds = get_stable_empty_site_inds(base_admol, coad_stable_sites)
        if len(self.stable_inds) == 0:
            raise ValueError("no empty stable coadsorbate sites found on the base structure")
        self._dist = self._build_site_distances()

    # ------------------------------------------------------------------ graph helpers
    def _build_site_distances(self):
        """BFS site-to-site hop distances on the base lattice graph, from every stable site to
        every site. Returns {stable_ind: {site_ind: hops}}. Paths may pass through occupied or
        central-species sites, which is the correct lattice distance."""
        base = self.base_admol
        site_inds = [i for i, a in enumerate(base.atoms) if a.is_surface_site()]
        atom_to_ind = {base.atoms[i]: i for i in site_inds}
        nbrs = {i: [atom_to_ind[a2] for a2 in base.atoms[i].bonds.keys()
                    if a2.is_surface_site()] for i in site_inds}
        dist = {}
        for s in self.stable_inds:
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

    def _empty(self, occ):
        return [i for i in self.stable_inds if i not in occ]

    def _nearby_empty(self, sind, occ):
        d = self._dist[sind]
        return [i for i in self.stable_inds
                if i not in occ and d.get(i, np.inf) <= self.local_radius]

    # ------------------------------------------------------------------ config build/eval
    def build_config(self, occ):
        """Build the admol for an occupancy set by adding a coadsorbate to each occupied site.
        Site indices are preserved across add_ad_to_site (coad atoms are appended), so
        base_admol.atoms[sind] remains the correct site throughout. May raise on a placement
        that fails to form a valid molecule; callers treat that as an invalid proposal."""
        m = self.base_admol
        for sind in occ:
            m = add_ad_to_site(m, self.coad, m.atoms[sind])
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
        return configuration_is_valid(m, self.base_admol, self.is_ts, self.unstable_pairs)

    # ------------------------------------------------------------------ initialization
    def _greedy_valid_state(self, Ncoad, max_restarts=200):
        """Greedily place Ncoad coadsorbates on empty sites, keeping the config valid at each
        step. Restarts on a dead end. Returns (occ_set, admol)."""
        for _ in range(max_restarts):
            occ = set()
            m = self.base_admol
            stuck = False
            while len(occ) < Ncoad:
                empties = self._empty(occ)
                random.shuffle(empties)
                placed = False
                for e in empties:
                    cand = occ | {e}
                    try:
                        mc = self.build_config(cand)
                    except Exception:
                        continue
                    if self.is_valid(mc):
                        occ, m, placed = cand, mc, True
                        break
                if not placed:
                    stuck = True
                    break
            if not stuck and len(occ) == Ncoad:
                return occ, m
        raise RuntimeError("could not build a valid starting configuration with Ncoad={} "
                           "(stable sites available: {})".format(Ncoad, len(self.stable_inds)))

    # ------------------------------------------------------------------ the sampler
    def run(self, Ncoad, n_steps=10000, T=(2000.0, 300.0), lam=(1.0, 0.0), p_local=0.8):
        """Run canonical Metropolis MC at fixed Ncoad with annealed T and lambda.

        Args:
            Ncoad: number of coadsorbates (fixed throughout the chain)
            n_steps: number of MC steps
            T: temperature in K; a (T_start, T_end) tuple is geometrically annealed, a scalar
                is held constant
            lam: uncertainty weight in E_corr = E - lam*sigma; a (lam_start, lam_end) tuple is
                linearly annealed, a scalar is held constant
            p_local: probability a relocation targets a nearby (within local_radius) empty site
                rather than any empty site

        Returns a dict with: Ncoad, best=(admol, E), pool=[(admol, E, sigma, trace), ...]
        (unique validly-evaluated configs, the candidate set for selection), and run
        diagnostics (n_steps, n_accept, accept_frac, n_unique, Emin).
        """
        if Ncoad > len(self.stable_inds):
            raise ValueError("Ncoad={} exceeds available stable sites ({})".format(
                Ncoad, len(self.stable_inds)))
        T0, T1 = (T, T) if np.isscalar(T) else T
        l0, l1 = (lam, lam) if np.isscalar(lam) else lam

        occ, m_cur = self._greedy_valid_state(Ncoad)
        E_cur, s_cur, tr_cur = self.energy(m_cur)

        visited = {}  # frozenset(occ) -> (admol, E, sigma, trace)

        def record(occ_set, m, E, s, tr):
            k = frozenset(occ_set)
            if k not in visited:
                visited[k] = (m, E, s, tr)

        record(occ, m_cur, E_cur, s_cur, tr_cur)

        n_accept = 0
        for step in range(n_steps):
            f = step / max(n_steps - 1, 1)
            Tk = T0 * (T1 / T0) ** f if (T0 > 0 and T1 > 0) else (T0 + (T1 - T0) * f)
            lam_s = l0 + (l1 - l0) * f

            # propose a canonical relocation: move one coad from an occupied to an empty site
            o = random.choice(tuple(occ))
            if random.random() < p_local:
                cands = self._nearby_empty(o, occ) or self._empty(occ)
            else:
                cands = self._empty(occ)
            if not cands:
                continue
            e = random.choice(cands)
            new_occ = (occ - {o}) | {e}

            try:
                m_new = self.build_config(new_occ)
            except Exception:
                continue
            if not self.is_valid(m_new):
                continue

            E_new, s_new, tr_new = self.energy(m_new)
            record(new_occ, m_new, E_new, s_new, tr_new)  # keep every validly-evaluated proposal

            dEcorr = (E_new - lam_s * s_new) - (E_cur - lam_s * s_cur)
            if dEcorr <= 0 or random.random() < np.exp(-dEcorr / (R * Tk)):
                occ, m_cur, E_cur, s_cur, tr_cur = new_occ, m_new, E_new, s_new, tr_new
                n_accept += 1

        pool = list(visited.values())
        best = min(pool, key=lambda v: v[1])
        return {"Ncoad": Ncoad, "best": (best[0], best[1]), "pool": pool,
                "n_steps": n_steps, "n_accept": n_accept,
                "accept_frac": n_accept / n_steps if n_steps else 0.0,
                "n_unique": len(pool), "Emin": best[1]}

    def scan_coverages(self, Ncoads=None, **run_kwargs):
        """Run an MC chain at each coverage level and assemble per-coverage minima plus a
        combined candidate pool. ``Ncoads`` defaults to 1..len(stable_inds).

        Returns a dict with Ncoad_energy_dict / Ncoad_config_dict (per-coverage best energy and
        config adjacency list), a combined ``pool`` for selection, and per-Ncoad ``results``.
        """
        if Ncoads is None:
            Ncoads = range(1, len(self.stable_inds) + 1)
        results = {}
        Ncoad_energy_dict = {}
        Ncoad_config_dict = {}
        pool = []
        for N in Ncoads:
            r = self.run(N, **run_kwargs)
            results[N] = r
            Ncoad_energy_dict[N] = r["best"][1]
            Ncoad_config_dict[N] = r["best"][0].to_adjacency_list()
            pool.extend(r["pool"])
            logging.info("coverage_mc: Ncoad=%d Emin=%.1f J/mol unique=%d accept=%.2f",
                         N, r["Emin"], r["n_unique"], r["accept_frac"])
        return {"results": results, "Ncoad_energy_dict": Ncoad_energy_dict,
                "Ncoad_config_dict": Ncoad_config_dict, "pool": pool}
