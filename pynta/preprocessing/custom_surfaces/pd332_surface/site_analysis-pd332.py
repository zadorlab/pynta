# Custom surface site analysis
# Edit **Cell 2 only**, then **Run All**.

# Imports (DO NOT EDIT)
import sys, os
sys.path.append("/home/shikim/pynta")

from ase.io import read, write, Trajectory
from ase import Atoms
from acat.settings import CustomSurface
from acat.adsorption_sites import SlabAdsorptionSites
from ase.build import surface
from ase.visualize import view

import site_analysis as sa
from site_analysis import get_occupied_sites, sites_match
import numpy as np

# USER INPUT
name = '/home/shikim/pynta/pynta/preprocessing/custom_surfaces/pd332_surface/pd332.xyz'
n_layers = 10
#atoms = surface('Pd',(3,3,2),n_layers) #lattice constant
#atoms.center(vacuum=10, axis=2)


# Load slab
slab = read(name)
surface = CustomSurface(slab, n_layers=n_layers)
nslab = len(slab)


adsorbate_height = 2
site_bond_cutoff = 2

#slabrat.rattle(stdev=0.3)

# Generate symmetry-unique site geometries
#cas = SlabAdsorptionSites(slab, "fcc332", composition_effect=True)  # ACAT has surface, custom does not find them all!
cas = SlabAdsorptionSites(slab, surface, composition_effect=True, allow_6fold=True)  # ACAT has surface, custom does not find them all!
#cas = SlabAdsorptionSites(slabrat, surface, composition_effect=True)  # ACAT has surface, custom does not find them all!

sites = cas.get_sites()
occ = get_occupied_sites(slab, sites, nslab, site_bond_cutoff)
unocc = [s for s in sites if not any(sites_match(s, o, slab) for o in occ)]

print(f"Total sites detected: {len(sites)}")
print(f"Occupied sites (within {site_bond_cutoff} Å of atoms): {len(occ)}")
print(f"Unoccupied sites: {len(unocc)}")

single_geoms, single_sites_lists = sa.generate_unique_sites(
    slab,
    cas.get_sites(),
    nslab,
    site_bond_cutoff,
    adsorbate_height
)

print(f"Symmetry-unique sites (after clustering): {len(single_sites_lists)}")

print(f'There are {len(single_sites_lists)} unique sites out of {len(cas.get_sites())}.')
#print(cas.get_sites())

#traj = Trajectory("unique_sites.traj", "w")
#for g in single_geoms:
#    traj.write(g)
#traj.close()

#do not rattle the surface. it distorts the surface

# ... your code that makes single_sites_lists ...

sa.save_sites_to_json(single_sites_lists, filename="single_sites_lists.json")

# Load the trajectory file

# Load the trajectory file
#unique_sites = read('unique_sites.traj', index=1)

#test
#cas = SlabAdsorptionSites(slab, surface, composition_effect=True)
#single_geoms, single_sites_lists = sa.generate_unique_sites(
#    slab,
#    cas.get_sites(),
#    nslab,
#    site_bond_cutoff,
#    adsorbate_height
#)

#print(f'There are {len(single_sites_lists)} unique sites out of {len(cas.get_sites())}.')
##print(cas.get_sites())
##sa.save_sites_to_json(single_sites_lists, filename="single_sites_lists_customsurf.json")
##
##traj = Trajectory("unique_sites.traj", "w")
##for g in single_geoms:
##    traj.write(g)
##traj.close()

#print(cas.get_sites())
#len(cas.get_sites())

#[cas.get_sites()[i]['morphology'] for i in range(len(cas.get_sites()))]

# Extract 3-fold site graphs
admols, threefold_geom_indices = sa.classify_all_sites(
    single_geoms, single_sites_lists
)

# Graph isomorphism clustering
iso_mat, clusters = sa.cluster_isomorphic_graphs(admols)

# Update 3-fold site labels
type_map = sa.update_threefold_site_labels(
    single_sites_lists,
    clusters,
    threefold_geom_indices
)

# Write 3-fold-only trajectory
traj3 = Trajectory("threefold_sites.traj", "w")
for i in threefold_geom_indices:
    traj3.write(single_geoms[i])
traj3.close()

# Pairwise strict isomorphism (PRINT)
print("Pairwise strict isomorphism:")
for i in range(len(admols)):
    for j in range(i + 1, len(admols)):
        print(f"{i} vs {j} =", iso_mat[i, j])

# Distinct 3-fold site types (PRINT)
print("Number of distinct 3-fold site types:", len(clusters))
for members in clusters.values():
    print("3-fold site type:", members)

# Updated 3-fold site labels (PRINT)
print("Updated 3-fold site labels per geometry:")
for geom_idx, label in type_map.items():
    print(f"Geometry {geom_idx} -> {label}")

# Site yaml file generated
print("All sites for the custom surfaces are saved in site.yaml")
#sa.write_sites_yaml(single_sites_lists, clusters)

# Compare with "fcc332"
print("\n--- Comparison with 'fcc332' ---")
cas_fcc = SlabAdsorptionSites(slab, "fcc332", composition_effect=True)
sites_fcc = cas_fcc.get_sites()
single_geoms_fcc, single_sites_lists_fcc = sa.generate_unique_sites(
    slab,
    sites_fcc,
    nslab,
    site_bond_cutoff,
    adsorbate_height
)

print(f"CustomSurface: {len(single_sites_lists)} unique sites out of {len(cas.get_sites())} total")
print(f"'fcc332': {len(single_sites_lists_fcc)} unique sites out of {len(sites_fcc)} total")

# Identify overlapping sites
overlapping_sites = []
tolerance = 0.1  # Angstroms

# Flatten the lists with geom indices
custom_sites_flat = []
custom_geom_indices = []
for idx, sublist in enumerate(single_sites_lists):
    for site in sublist:
        custom_sites_flat.append(site)
        custom_geom_indices.append(idx)

fcc_sites_flat = []
fcc_geom_indices = []
for idx, sublist in enumerate(single_sites_lists_fcc):
    for site in sublist:
        fcc_sites_flat.append(site)
        fcc_geom_indices.append(idx)

# Get excluded sites (symmetry duplicates)
from pynta.utils import get_unique_sym_struct_index_clusters
from pynta.mol import add_adsorbate_to_site

# Recreate unocc and geoms as in generate_unique_sites
geoms = []
sites_per_geom = []
for site in unocc:
    g = slab.copy()
    add_adsorbate_to_site(
        g,
        adsorbate=Atoms("He"),
        surf_ind=0,
        site=site,
        height=adsorbate_height,
        tilt_angle=25.24,
        offset=False)
    geoms.append(g)
    sites_per_geom.append([site])

clusters = get_unique_sym_struct_index_clusters(geoms)
unique_indices = set(c[0] for c in clusters)
excluded_sites = [sites_per_geom[i][0] for i in range(len(geoms)) if i not in unique_indices]

print(f"Excluded symmetry-duplicate sites: {len(excluded_sites)}")

# Check overlap of excluded sites with fcc332 sites
excluded_overlapping = []
for site_ex in excluded_sites:
    pos_ex = np.array(site_ex['position'])
    for site_fcc in fcc_sites_flat:
        pos_fcc = np.array(site_fcc['position'])
        dist = np.linalg.norm(pos_ex - pos_fcc)
        if dist < tolerance:
            excluded_overlapping.append(site_ex)
            break

print(f"Excluded sites overlapping with 'fcc332': {len(excluded_overlapping)} out of {len(excluded_sites)}")

# Map overlapping custom and fcc unique site geometries
overlapping_custom_indices = set()
overlapping_fcc_indices = set()
for i, site_custom in enumerate(custom_sites_flat):
    pos_custom = np.array(site_custom['position'])
    for j, site_fcc in enumerate(fcc_sites_flat):
        pos_fcc = np.array(site_fcc['position'])
        dist = np.linalg.norm(pos_custom - pos_fcc)
        if dist < tolerance:
            overlapping_sites.append((site_custom, site_fcc))
            overlapping_custom_indices.add(custom_geom_indices[i])
            overlapping_fcc_indices.add(fcc_geom_indices[j])
            break  # Only count once per custom site

print(f"Number of overlapping sites: {len(overlapping_custom_indices)} unique custom sites out of {len(custom_sites_flat)} custom sites")

# Write trajectory files for all identified sets
traj_custom = Trajectory("customsurface_sites.traj", "w")
for geom in single_geoms:
    traj_custom.write(geom)
traj_custom.close()

traj_fcc = Trajectory("fcc332_sites.traj", "w")
for geom in single_geoms_fcc:
    traj_fcc.write(geom)
traj_fcc.close()

traj_overlap = Trajectory("overlapping_sites.traj", "w")
for idx in sorted(overlapping_custom_indices):
    traj_overlap.write(single_geoms[idx])
traj_overlap.close()

print("Wrote trajectories: customsurface_sites.traj, fcc332_sites.traj, overlapping_sites.traj")

# Save fcc332 sites
sa.save_sites_to_json(single_sites_lists_fcc, filename="single_sites_lists_fcc332.json")