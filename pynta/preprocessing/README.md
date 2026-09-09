# Preprocessing

Tools for preparing DFT inputs and analyzing surface sites before running Pynta.

## Files

| File | Purpose |
|------|---------|
| `preprocessing.py` | All classes, workflow steps, alloy helpers, and plot functions |
| `runprep.py` | Driver script — configure and run `Prep` here |
| `run.prep.sh` | SLURM batch job script that calls `runprep.py` |
| `Pynta preprocessing.ipynb` | Postprocessing notebook — reads `prep_convergence.json` and plots results |
| `site_analysis.py` | Python module for adsorption-site analysis, imported by `preprocessing.py` and `site_analysis.ipynb` |
| `site_analysis.ipynb` | Site-analysis notebook — edit Cell 2, then Run All |
| `custom_surfaces/` | Example slab `.xyz` files used by `site_analysis.ipynb` |

---

## Part 1 — DFT Convergence Workflow

DFT preprocessing workflow for Pynta: lattice constant optimization, slab generation, and convergence testing (ecut and k-points). Supports pure metals and binary alloys (fcc/bcc).

### Quick Start

Edit `runprep.py` with your system parameters and run:

```bash
python runprep.py
# or on a cluster:
sbatch run.prep.sh
```

The workflow is **restart-safe**: completed steps are skipped on re-run by reading `prep_convergence.json`.

### How It Works

`Prep` uses two decorators to define which steps run automatically:

- `@requires("keyword")` — a step only runs if that keyword was passed to `Prep()`
- `@step(n)` — sets the execution order

Pass a keyword to opt into a step:

```python
prep = Prep(
    optimize_lattice=True,       # activates calculate_lattice_parameters
    generate_slab=True,          # activates generate_slab
    ecut_values=[40, 60, 80],    # activates run_convergence (with kmesh_values)
    kmesh_values=[(3,3,1), (5,5,1)],
    ...
)
prep.run_selected()
```

### `Prep` Parameters

#### System

| Parameter | Default | Description |
|-----------|---------|-------------|
| `metal` | `"Pt"` | Element symbol |
| `surface_type` | `"fcc111"` | ASE surface builder name (`fcc111`, `bcc110`, `hcp0001`, …) |
| `a0` | `3.96` | Initial lattice constant guess (Å) |
| `c` | `None` | c/a ratio for hcp; omit for cubic |
| `repeats` | `(3,3,4)` | Slab supercell size `(nx, ny, nz)` |
| `vacuum` | `10` | Vacuum thickness (Å) |
| `pbc` | `(True,True,False)` | Periodic boundary conditions |
| `fmax` | `0.05` | Force convergence threshold for slab relaxation (eV/Å) |
| `frozen_layers` | `None` | Number of bottom slab layers to freeze |

#### Calculator

| Parameter | Default | Description |
|-----------|---------|-------------|
| `software` | `"Espresso"` | ASE calculator name passed to `name_to_ase_software()` |
| `software_kwargs` | see below | Calculator kwargs used for slab and convergence calculations |
| `lattice_opt_software_kwargs` | see below | Calculator kwargs used only for lattice optimization (denser k-mesh, higher ecut) |

**Default `software_kwargs`** (Quantum ESPRESSO / BEEF-vdW):
```python
{
    "kpts": (3, 3, 1),
    "ecutwfc": 40,
    "occupations": "smearing",
    "smearing": "marzari-vanderbilt",
    "degauss": 0.01,
    "input_dft": "BEEF-VDW",
    "nosym": True,
    "conv_thr": 1e-6,
    "pseudopotentials": {...},
}
```

**Default `lattice_opt_software_kwargs`**:
```python
{
    "kpts": (25, 25, 25),
    "ecutwfc": 70,
    "ecutrho": 280,
    "degauss": 0.02,
    "input_dft": "BEEF-VDW",
    "mixing_mode": "plain",
    "pseudopotentials": <inherited from software_kwargs>,
}
```

#### Workflow flags

Pass these to activate the corresponding step:

| Keyword | Activates |
|---------|-----------|
| `optimize_lattice=True` | `calculate_lattice_parameters` |
| `optimize_alloy_lattice=True` | `calculate_alloy_lattice_parameters` |
| `generate_slab=True` | `generate_slab` |
| `ecut_values=[...]` + `kmesh_values=[...]` | `run_convergence` |
| `run_adsorbate_convergence=True` + `ecut_values` + `kmesh_values` | `run_convergence_adsorbate` |

#### Alloy parameters

| Parameter | Description |
|-----------|-------------|
| `alloy_path` | Path to a pre-built alloy structure file (xyz/cif/…). If provided, skips the ordered-cell builder. |
| `alloy_A` | Majority element symbol (e.g. `"Pt"`) |
| `alloy_B` | Minority element symbol (e.g. `"Au"`) |
| `alloy_xA` | Fraction of A (e.g. `0.75` for 75% Pt) |
| `alloy_maxrep` | Max cubic repeat to match composition exactly (default `6`) |
| `alloy_scan_step` | Lattice scan step in Å (default `0.01`) |

### Workflow Steps

#### 1. `calculate_alloy_lattice_parameters` — alloy lattice constant
Activated by `optimize_alloy_lattice=True`.

Builds an ordered alloy supercell (smallest cubic repeat that exactly matches `xA`) or reads from `alloy_path`, then scans total energy vs. lattice constant *a* and fits a quadratic to find the minimum. Sets `self.a`. Supports fcc and bcc only.

#### 2. `calculate_lattice_parameters` — pure metal lattice constant
Activated by `optimize_lattice=True`.

Scans a cubic bulk cell over a range of *a* values, fits a quadratic, and sets `self.a`.

#### 3. `generate_slab`
Activated by `generate_slab=True`.

Builds a slab using `self.a` and `surface_type` via ASE's surface builders, relaxes it with BFGS, and writes `slab.xyz` and `slab.traj`.

#### 4. `run_convergence`
Activated by passing both `ecut_values` and `kmesh_values`.

Loops over all `(ecut, kpts)` combinations, computes total energy of the slab, and records results.

#### 5. `run_convergence_adsorbate`
Activated by `run_adsorbate_convergence=True` + `ecut_values` + `kmesh_values`.

Same loop as `run_convergence`, but computes adsorption energy:
```
E_ads = E(slab+adsorbate) - E(slab) - E(reference)
```
Defaults to H with H₂ as reference. Pass `adsorbate="CO"` etc. for other species.

### Output

All results are written to **`prep_convergence.json`** after each step. This file drives restart behaviour and is read by the postprocessing notebook.

```
prep_convergence.json
├── created          — timestamp
├── user_arguments   — parameters passed to Prep()
└── steps[]
    ├── step         — step name
    ├── time         — completion timestamp
    ├── status       — "completed" or "failed"
    ├── inputs       — parameters used
    └── outputs      — results (lattice_constant, scan_avals_A, scan_Evals_eV, results, …)
```

Log messages are written to `preprocessing.log` and stdout simultaneously.

### Postprocessing

Open `Pynta preprocessing.ipynb`. It reads `prep_convergence.json` and calls three functions from `preprocessing.py`:

```python
from preprocessing import plot_lattice_constant, plot_ecut_convergence, plot_kpoint_convergence

a_opt            = plot_lattice_constant(lout)
ecuts_all, kpts_all = plot_ecut_convergence(results)
                     plot_kpoint_convergence(results)
```

Each function produces a dual-panel plot (absolute energy + |ΔE| on log scale) and prints a per-point table.

### Alloy Example

```python
prep = Prep(
    metal="Pt",
    surface_type="fcc111",
    software="Espresso",
    a0=3.95,
    optimize_alloy_lattice=True,
    alloy_A="Pt",
    alloy_B="Au",
    alloy_xA=0.75,
    software_kwargs={...},
    lattice_opt_software_kwargs={...},
)
prep.run_selected()
```

Or with a pre-built structure:
```python
prep = Prep(
    ...
    optimize_alloy_lattice=True,
    alloy_path="my_PtAu_alloy.xyz",
    a0=3.95,
)
```

---

## Part 2 — Custom Surface Site Analysis

Tools for identifying symmetry-unique adsorption sites on non-standard surfaces using graph isomorphism. Works with stepped metal surfaces, surfaces with vacancies, oxide overlayers, and any geometry representable as an ASE `Atoms` object.

### Example slabs (in `custom_surfaces/`)

| File | Surface type |
|---|---|
| `pd332.xyz` | Pd(332) stepped metal surface |
| `pt111_no_defect.xyz` | Pt(111) flat terrace |
| `pt111_defect_middle.xyz` | Pt(111) with a single vacancy |
| `rutile_110.xyz` | TiO₂ rutile(110) oxide surface |
| `alpha_PtO2_layer.xyz` | α-PtO₂ oxide layer |
| `fe2o3.xyz` | Fe₂O₃ oxide surface |
| `single_atom_PtAu.xyz` | Single-atom PtAu alloy slab |

### Quick start

Open `site_analysis.ipynb` and edit **Cell 2 only**:

```python
xyz_path = "custom_surfaces/pd332.xyz"   # path to your slab
n_layers = 10                            # number of atomic layers (use suggest_n_layers helper if unsure)
facet    = None                          # "fcc111", "fcc332", etc. — or None for CustomSurface
```

Then **Run All**. The workflow auto-detects whether the slab has a vacancy and dispatches accordingly. If a vacancy is detected, the "Defect Site Analysis" cells run automatically and produce `panel1_void_ring.traj`, `panel2_drop_steps.traj`, `panel3_void_sites.traj` and matching PNG figures.

### Workflow

`workflow_auto` has two branches:

**No defect** (`workflow_no_defect_unique_sites`)
1. Enumerate all ACAT adsorption sites
2. Tag each site with a noble-gas atom at `adsorbate_height`
3. Build RMG graphs, cluster by graph isomorphism + (site type, morphology)
4. Assign labels (`ontop0`, `bridge0`, `3fold0`, …) and reduce to one representative per label

**Defect / vacancy** (`workflow_defect_vacancy_drop`)
1. Detect vacancy positions from coordination number analysis
2. For each terrace site: drop a tagged atom toward the vacancy, recording the trajectory
3. Find the position of maximum coordination number as the "drop site"
4. Apply graph-isomorphism labeling to all found sites

### Output files

Both workflows write output files to the **current working directory** (wherever the notebook is run from):

| File | Content |
|---|---|
| `sites.json` | All raw ACAT sites |
| `labeled_sites.json` | Representative sites with graph-isomorphism labels |
| `neighbor_site_list.json` | Neighbor relationships between sites |
| `unique_sites.traj` | One geometry per unique site *(no-defect path)* |
| `geom_all.traj` | All site geometries *(defect path)* |
| `drop_steps.traj` | Full drop trajectories *(defect path)* |
| `maxcn_geoms.xyz` | Max-coordination-number positions per defect *(defect path)* |
| `panel1_void_ring.traj`, `panel2_drop_steps.traj`, `panel3_void_sites.traj` + matching `.png` | Defect Site Analysis panels |

### Key parameters

| Parameter | Default | Description |
|---|---|---|
| `n_layers` | 4 | Number of atomic layers for `CustomSurface`; use `sa.suggest_n_layers(slab)` to check |
| `facet` | `None` | ACAT facet string (`"fcc111"`, `"fcc332"`, …); `None` uses `CustomSurface` |
| `adsorbate_height` | 1.0 Å | Height above surface to place the tag atom |
| `site_bond_cutoff` | 1.5 Å | Distance cutoff for site–slab bond detection |
| `tag_symbol` | `"Ne"` | Noble-gas symbol used to tag sites |
| `dz` | 0.1 Å | Step size for the vacancy-drop algorithm |
| `stable_steps` | 3 | Patience: stop dropping after this many steps with no CN increase |
| `max_drop` | 25.0 Å | Defect analysis only — must exceed void depth + 3 Å starting offset |
| `min_clearance` | 1.3 Å | Defect analysis only — stop the probe drop this close to any slab atom |
| `margin` | 0.25 Å | Defect analysis only — added to the nearest-neighbor distance to set the CN cutoff radius |
| `uniq_tol` | 0.25 Å | Defect analysis only — deduplication tolerance for unique probe positions |

### Dependencies

- [ASE](https://wiki.fysik.dtu.dk/ase/)
- [ACAT](https://github.com/HongxuLiu/acat) — `CustomSurface`, `SlabAdsorptionSites`
- [RMG-Py](https://github.com/ReactionMechanismGenerator/RMG-Py) — graph isomorphism via `Molecule`
- NetworkX, NumPy, Matplotlib
