# Custom Surface Site Analysis

Tools for identifying symmetry-unique adsorption sites on non-standard surfaces using graph isomorphism. Works with stepped metal surfaces, surfaces with vacancies, oxide overlayers, and any geometry representable as an ASE `Atoms` object.

## Files

| File | Description |
|---|---|
| `site_analysis.ipynb` | Main notebook — edit Cell 2, then Run All |
| `site_analysis.py` | Python module imported by the notebook |
| `pd332.xyz` | Pd(332) stepped metal surface |
| `pt111_no_defect.xyz` | Pt(111) flat terrace |
| `pt111_defect_middle.xyz` | Pt(111) with a single vacancy |
| `rutile_110.xyz` | TiO₂ rutile(110) oxide surface |
| `alpha_PtO2_layer.xyz` | α-PtO₂ oxide layer |
| `fe2o3.xyz` | Fe₂O₃ oxide surface |

## Quick start

Open `site_analysis.ipynb` and edit **Cell 2 only**:

```python
xyz_path = "pd332.xyz"   # path to your slab
n_layers = 10            # number of atomic layers (use suggest_n_layers helper if unsure)
facet    = None          # "fcc111", "fcc332", etc. — or None for CustomSurface
```

Then **Run All**. The workflow auto-detects whether the slab has a vacancy and dispatches accordingly.

## Workflow

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

## Output files

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

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `n_layers` | 4 | Number of atomic layers for `CustomSurface`; use `sa.suggest_n_layers(slab)` to check |
| `facet` | `None` | ACAT facet string (`"fcc111"`, `"fcc332"`, …); `None` uses `CustomSurface` |
| `adsorbate_height` | 1.0 Å | Height above surface to place the tag atom |
| `site_bond_cutoff` | 1.5 Å | Distance cutoff for site–slab bond detection |
| `tag_symbol` | `"Ne"` | Noble-gas symbol used to tag sites |
| `dz` | 0.1 Å | Step size for the vacancy-drop algorithm |
| `stable_steps` | 3 | Patience: stop dropping after this many steps with no CN increase |

## Dependencies

- [ASE](https://wiki.fysik.dtu.dk/ase/)
- [ACAT](https://github.com/HongxuLiu/acat) — `CustomSurface`, `SlabAdsorptionSites`
- [RMG-Py](https://github.com/ReactionMechanismGenerator/RMG-Py) — graph isomorphism via `Molecule`
- NetworkX, NumPy, Matplotlib
