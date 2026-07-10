"""
pes_diagram.py  —  single-file potential-energy-surface (PES) plotter
=====================================================================
Draw reaction-coordinate diagrams from a simple ordered list of states.

Two ways to reference energies:
  * per-state, by giving one state's label to `reference` (it becomes 0 eV); or
  * per-PATH, by giving each path a `gas_energy` (the energy of that branch's
    gas-phase reference species on a COMMON scale). The plotter then sets the
    LOWEST gas-phase species to 0 and shifts every higher branch up by the gap,
    so overlaid branches share one absolute zero and their heights are directly
    comparable.

A path is an ordered list of states, each one of:
    ("label", energy)            -> a minimum (reactant / intermediate / product)
    ("label", energy, "ts")      -> a transition state (barrier peak)

Run this file (`python pes_diagram.py`) to produce the CH4-vs-NH3 overlay.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence, Union

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ============================================================== data model ===

@dataclass
class _State:
    label: str
    energy: float
    kind: str = "min"           # "min" or "ts"


@dataclass
class _Path:
    name: str
    states: list                # list[_State]
    color: str = "black"
    linewidth: float = 4.0
    label_side: str = "below"   # "below" or "above" — where names/energies go
    gas_energy: float = 0.0     # energy of this branch's gas-phase reference,
                                # on a scale common to all paths


@dataclass
class _BranchLevel:
    attach: Union[int, str]
    energy: float
    label: str = ""


@dataclass
class _Branch:
    name: str
    levels: list
    color: str = "green"
    linewidth: float = 3.0


def _coerce_state(item) -> _State:
    """Accept ('label', E) or ('label', E, 'ts')."""
    if isinstance(item, _State):
        return item
    if len(item) == 2:
        label, energy = item
        return _State(label, float(energy), "min")
    label, energy, kind = item
    kind = "ts" if str(kind).lower() in ("ts", "tst", "barrier", "saddle") else "min"
    return _State(label, float(energy), kind)


# ================================================================ helper ======

def build_ladder(start_label, steps, adsorbate_label=None, E_ads=0.0):
    """
    Chain elementary steps into a cumulative profile, reactant at 0 eV.

    steps : list of (Ea, dE, product_label) — forward barrier and step energy.

    If `adsorbate_label` is given, the profile starts at the gas-phase species
    (`start_label`, 0 eV), drops to the adsorbed reactant at -|E_ads|, then
    chains the steps from there. Otherwise it starts at `start_label` = 0.

        E(TS_k)    = E(state_{k-1}) + Ea_k
        E(state_k) = E(state_{k-1}) + dE_k
    """
    if adsorbate_label is not None:
        states = [(start_label, 0.0), (adsorbate_label, -abs(E_ads))]
        E = -abs(E_ads)
    else:
        states = [(start_label, 0.0)]
        E = 0.0
    for Ea, dE, product in steps:
        states.append(("TS", E + Ea, "ts"))
        E += dE
        states.append((product, E))
    return states


# ============================================================== main class ===

class PESDiagram:
    def __init__(self, ylabel="Energy (eV)", reference=None,
                 bar_halfwidth=0.34, slot_width=1.0, x_start=0.6):
        self.ylabel = ylabel
        self.reference = reference
        self.hw = bar_halfwidth
        self.dx = slot_width
        self.x0 = x_start
        self._paths: list[_Path] = []
        self._branches: list[_Branch] = []

    # ---- building ----------------------------------------------------------

    def add_path(self, name, states, color="black", linewidth=4.0,
                 label_side="below", gas_energy=0.0):
        """gas_energy: energy of this branch's gas-phase reference species on a
        scale shared with the other paths. Lowest across all paths -> 0 eV."""
        self._paths.append(_Path(name, [_coerce_state(s) for s in states],
                                 color, linewidth, label_side, float(gas_energy)))
        return self

    def add_branch(self, name, levels, color="green", linewidth=3.0):
        parsed = []
        for lv in levels:
            if len(lv) == 2:
                attach, energy = lv; label = ""
            else:
                attach, energy, label = lv
            parsed.append(_BranchLevel(attach, float(energy), label))
        self._branches.append(_Branch(name, parsed, color, linewidth))
        return self

    # ---- internals ---------------------------------------------------------

    def _shift(self):
        if self.reference is None:
            return 0.0
        for p in self._paths:
            for s in p.states:
                if s.label == self.reference:
                    return -s.energy
        raise ValueError(f"reference label {self.reference!r} not found in any path.")

    def _xcenter(self, i):
        return self.x0 + i * self.dx

    def _resolve_attach(self, main, attach):
        if isinstance(attach, int):
            return attach
        for i, s in enumerate(main.states):
            if s.label == attach:
                return i
        raise ValueError(f"branch attach label {attach!r} not found on main path.")

    # ---- plotting ----------------------------------------------------------

    def plot(self, figsize=(11, 6), save=None, dpi=200,
             show_reference_line=True, zero_line=False, zero_label=None,
             show_gas_gap=True, show_ts_labels=False, show_branch_labels=False,
             label_fontsize=9, energy_fontsize=9, ypad=0.18, ax=None):

        shift = self._shift()
        gas_min = min((p.gas_energy for p in self._paths), default=0.0)

        def poff(p):                       # total vertical offset for a path
            return shift + (p.gas_energy - gas_min)

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        # collect all plotted energies to size the axis / label offsets
        all_e = [s.energy + poff(p) for p in self._paths for s in p.states]
        main = self._paths[0] if self._paths else None
        if main is not None:
            all_e += [lv.energy + poff(main) for b in self._branches for lv in b.levels]
        ymin, ymax = min(all_e), max(all_e)
        yr = (ymax - ymin) or 1.0
        name_dy = 0.045 * yr
        en_dy_min = 0.105 * yr
        en_dy_ts = 0.055 * yr

        # --- main path(s): bars + connectors ---
        for p in self._paths:
            off = poff(p)
            below = (p.label_side != "above")
            va = "top" if below else "bottom"
            sgn = -1 if below else 1

            xs_end = []
            for i, s in enumerate(p.states):
                xc = self._xcenter(i)
                y = s.energy + off
                ax.plot([xc - self.hw, xc + self.hw], [y, y], color=p.color,
                        lw=p.linewidth, solid_capstyle="round", zorder=3)
                xs_end.append((xc, y))

                if s.kind == "min":
                    ax.annotate(s.label, (xc, y + sgn * name_dy), ha="center",
                                va=va, fontsize=label_fontsize, fontweight="bold",
                                color=p.color, zorder=4)
                    ax.annotate(f"{y:.2f}", (xc, y + sgn * (name_dy + en_dy_min)),
                                ha="center", va=va, fontsize=energy_fontsize,
                                color=p.color, zorder=4)
                else:
                    if show_ts_labels and s.label:
                        ax.annotate(s.label, (xc, y + sgn * name_dy), ha="center",
                                    va=va, fontsize=label_fontsize, color=p.color,
                                    zorder=4)
                    ax.annotate(f"{y:.2f}", (xc, y + sgn * en_dy_ts), ha="center",
                                va=va, fontsize=energy_fontsize, color=p.color,
                                zorder=4)

            for (xc, y), (xc2, y2) in zip(xs_end, xs_end[1:]):
                ax.plot([xc + self.hw, xc2 - self.hw], [y, y2],
                        color=p.color, ls="--", lw=1.0, alpha=0.45, zorder=1)

        # --- branches attached to the main path ---
        if main is not None:
            off = poff(main)
            for b in self._branches:
                for lv in b.levels:
                    idx = self._resolve_attach(main, lv.attach)
                    xc = self._xcenter(idx)
                    y_branch = lv.energy + off
                    y_main = main.states[idx].energy + off
                    ax.plot([xc - self.hw, xc + self.hw], [y_branch, y_branch],
                            color=b.color, lw=b.linewidth, solid_capstyle="round",
                            zorder=3)
                    for xend in (xc - self.hw, xc + self.hw):
                        ax.plot([xend, xend], [y_main, y_branch],
                                color=b.color, ls="--", lw=1.0, zorder=1)
                    dn = y_branch < y_main
                    ydir = -1 if dn else 1
                    va = "top" if dn else "bottom"
                    ax.annotate(f"{y_branch:.2f}", (xc, y_branch + ydir * en_dy_ts),
                                ha="center", va=va, fontsize=energy_fontsize,
                                color=b.color, zorder=4)
                    if show_branch_labels and lv.label:
                        ax.annotate(lv.label, (xc, y_branch + ydir * (en_dy_ts + name_dy)),
                                    ha="center", va=va, fontsize=label_fontsize,
                                    color=b.color, zorder=4)

        # --- limits ---
        n = max((len(p.states) for p in self._paths), default=1)
        xmax = self._xcenter(n - 1) + self.x0 + self.hw
        ax.set_xlim(0, xmax)
        ax.set_ylim(ymin - ypad * yr - en_dy_min, ymax + ypad * yr)

        # --- gas-gap arrow between the two gas-phase references (slot 0) ---
        if show_gas_gap and len(self._paths) == 2:
            p_lo = min(self._paths, key=lambda p: p.gas_energy)
            p_hi = max(self._paths, key=lambda p: p.gas_energy)
            gap = p_hi.gas_energy - p_lo.gas_energy
            if abs(gap) > 1e-9:
                x = self._xcenter(0) - self.hw - 0.12
                y0 = p_lo.states[0].energy + poff(p_lo)
                y1 = p_hi.states[0].energy + poff(p_hi)
                ax.annotate("", xy=(x, y1), xytext=(x, y0),
                            arrowprops=dict(arrowstyle="<->", color="0.35", lw=1.2))
                ax.annotate(f"$\\Delta_{{gas}}$ = {gap:.2f} eV",
                            (x, 0.5 * (y0 + y1)), xytext=(-6, 0),
                            textcoords="offset points", ha="right", va="center",
                            rotation=90, fontsize=8, color="0.35")

        # --- reference / zero line ---
        if (show_reference_line and self.reference is not None) or zero_line \
                or any(p.gas_energy != gas_min for p in self._paths):
            ax.axhline(0.0, color="0.5", ls=":", lw=1.2, zorder=0)
            cap = zero_label or (f"{self.reference} = 0 eV" if self.reference
                                 else "lower gas-phase reference = 0 eV")
            ax.annotate(cap, (xmax, 0.0), xytext=(-4, 3),
                        textcoords="offset points", ha="right", va="bottom",
                        fontsize=8, color="0.4")

        ax.set_ylabel(self.ylabel, fontsize=14)
        ax.set_xlabel("Reaction coordinate", fontsize=13)
        ax.get_xaxis().set_ticks([])
        ax.grid(False)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

        handles = [mpatches.Patch(color=p.color, label=p.name) for p in self._paths]
        handles += [mpatches.Patch(color=b.color, label=b.name) for b in self._branches]
        if handles:
            ax.legend(handles=handles, loc="upper left", frameon=False)

        fig.tight_layout()
        if save:
            fig.savefig(save, dpi=dpi, bbox_inches="tight")
        return fig, ax


# ============================================================== examples ======

def hcn_overlay():
    """CH4 vs NH3, on a COMMON gas-phase reference (lower species = 0)."""
    # ----- USER INPUTS -------------------------------------------------------
    # Relative energies of the two gas-phase reference species, on whatever
    # common scale your DFT uses (e.g. E[X(g)+clean slab]). The lower one is set
    # to 0 automatically; the higher branch is shifted up by the difference.
    E_gas_CH4 = 0.00        # CH4(g) + *   (lower here -> becomes 0)
    E_gas_NH3 = 0.50        # NH3(g) + *   <-- PLACEHOLDER: put your value here
    # Adsorption energies (positive = binding strength):
    Eads_CH4 = 0.00         # methane: no chemisorption well (RPBE)
    Eads_NH3 = 0.60         # ammonia: ~0.6 eV (paper)
    # -------------------------------------------------------------------------

    # Dehydrogenation steps (Ea, dE, product) — DFT Table 1, JPCC 2011, 115, 5667
    ch4_steps = [
        (0.73,  0.04, "CH$_3$* + H*"),
        (0.89,  0.30, "CH$_2$* + 2H*"),
        (0.38, -0.43, "CH* + 3H*"),
        (1.51,  0.94, "C* + 4H*"),
    ]
    nh3_steps = [
        (1.39,  0.45, "NH$_2$* + H*"),
        (1.30,  0.19, "NH* + 2H*"),
        (1.32,  0.54, "N* + 3H*"),
    ]

    ch4 = build_ladder("CH$_4$(g)+*", ch4_steps, "CH$_4$*", Eads_CH4)
    nh3 = build_ladder("NH$_3$(g)+*", nh3_steps, "NH$_3$*", Eads_NH3)

    d = PESDiagram(reference=None)
    d.add_path("C-H activation (CH$_4$)", ch4, color="#1f3b73",
               label_side="below", gas_energy=E_gas_CH4)
    d.add_path("N-H activation (NH$_3$)", nh3, color="#b22222",
               label_side="above", gas_energy=E_gas_NH3)
    d.plot(save="hcn_overlay.png", zero_line=True, show_gas_gap=True)
    print("wrote hcn_overlay.png")


if __name__ == "__main__":
    hcn_overlay()
