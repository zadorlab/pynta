from catkit import Gratoms
from catkit.gen.adsorption import Builder
from catkit.build import molecule

from rmgcat_to_sella.adjacency_to_3d import get_edges, rmgcat_to_gratoms
from rmgcat_to_sella.find_all_nebs import get_all_species
from rmgcat_to_sella.ts import get_index_adatom
from rmgcat_to_sella.ts import get_index_surface_atom
from rmgcat_to_sella.ts import get_bond_dist
from rmgcat_to_sella.ts import get_av_dist

from ase.io import read, write

import os

import shutil

import yaml

def scan_bond_len_all(slab, repeats, yamlfile, facetpath, bondlen_increment, bondlen_cutoff):
    with open(yamlfile, 'r') as f:
        yamltxt = f.read()
    reactions = yaml.safe_load(yamltxt)

    species_unique = dict()
    nslab = len(slab)

    speciesInd = []
    bonds = []

    for rxn in reactions:
        # transforming reactions data to gratom objects
        reactants, rbonds = rmgcat_to_gratoms(rxn['reactant'].split('\n'))
        products, pbonds = rmgcat_to_gratoms(rxn['product'].split('\n'))
        speciesInd += reactants + products
        bonds += rbonds + pbonds

    r_unique = []
    p_unique = []
    symbols_list = []

    for rp, uniquelist in ((reactants, r_unique), (products, p_unique)):
        for species in rp:
            symbols = str(species.symbols)
            symbols_list.append(symbols)
            speciesdir = os.path.join(facetpath, 'minima_unique', symbols)
            if symbols not in species_unique:
                species_unique[symbols] = get_all_species(speciesdir)
            uniquelist.append(species_unique[symbols])


 = '+'.join([str(species.symbols) for species in reactants])
    p_name = '+'.join([str(species.symbols) for species in products])

    r_name_list = [str(species.symbols) for species in reactants]
    p_name_list = [str(species.symbols) for species in products]

    if len(r_name_list) < len(p_name_list):
        print('Reactant structure will be used to estimate to TS')
        rpDir = 'reactants'
    elif len(r_name_list) > len(p_name_list):
        print('Products will be used to estimate to TS')
        rpDir = 'products'

    saveDir = os.path.join(facetpath, 'scan_bond_all')
    if os.path.exists(saveDir):
        shutil.rmtree(saveDir)
        os.makedirs(saveDir)
    else:
        os.makedirs(saveDir)

    if len(r_name_list) < len(p_name_list):
        TS_candidate = molecule(r_name_list[0])[0]
    elif len(r_name_list) > len(p_name_list):
        TS_candidate = molecule(p_name_list[0])[0]

    if len(TS_candidate.get_tags()) < 3:
        TS_candidate.rotate(90, 'y')

    slabedges, tags = get_edges(slab, True)

    grslab = Gratoms(numbers=slab.numbers,
                     positions=slab.positions,
                     cell=slab.cell,
                     pbc=slab.pbc,
                     edges=slabedges)
    grslab.arrays['surface_atoms'] = tags

    # building adsorbtion structures
    ads_builder = Builder(grslab)

    count = 0
    bond_dist = TS_candidate.get_distance(0, 1)
    bond_dist_cutoff = bondlen_cutoff

    while bond_dist <= bond_dist_cutoff:
        TS_candidate.set_distance(0, 1, bond_dist, fix=0)
        
        structs = ads_builder.add_adsorbate(
            TS_candidate, [0], -1, auto_construct=False)
        big_slab = slab * repeats
        nslab = len(slab)
        rounded_bond = bond_dist.round(2)
        for i, struc in enumerate(structs):
            big_slab_ads = big_slab + struc[nslab:]
            write(os.path.join(saveDir, '{}'.format(
                str(i + len(structs) * count).zfill(3)) + '_' + p_name + '_{}'.format(str(rounded_bond).replace('.','_')) + '.png'), big_slab_ads)
            write(os.path.join(saveDir, '{}'.format(
                str(i + len(structs) * count).zfill(3)) + '_' + p_name + '_{}'.format(str(rounded_bond).replace('.','_')) + '.xyz'), big_slab_ads)

        # TS_candidate.rotate(rotAngle, 'z')
        bond_dist += bondlen_increment
        count += 1


def scan_bond(facetpath, position):
    all_geoms_list = []
    scanAllDir = os.path.join(facetpath, 'scan_bond_all')
    scanDir = os.path.join(facetpath, 'scan_bond_' + position)
    os.makedirs(scanDir, exist_ok=True)
    for geom in sorted(os.listdir(scanAllDir), key=str):
        if geom.endswith('.xyz'):
            all_geoms_list.append(geom)
    if position == 'top':
        count = 0
    elif position == 'bridge':
        count = 1
    elif position == 'hollow1':
        count = 2
    elif position == 'hollow2':
        count = 3
    selected_pos = []
    while count < len(all_geoms_list):
        selected_pos.append(all_geoms_list[count])
        count += 4
    for el in selected_pos:
        selected_src_xyz = os.path.join(scanAllDir, el)
        selected_src_png = os.path.join(selected_src_xyz[:-3] + 'png')
        shutil.copy2(selected_src_xyz, scanDir)
    scanGeomList = []
    for scanGeom in sorted(os.listdir(scanDir), key=str):
        if scanGeom.endswith('.xyz'):
            scanGeomList.append(scanGeom)
    for i, el in enumerate(scanGeomList):
        oldFname = os.path.join(scanDir, el)
        newFdir = os.path.join(scanDir, str(i).zfill(2))
        os.makedirs(newFdir, exist_ok=True)
        newFname = os.path.join(newFdir, '{}'.format(str(i).zfill(2)) + el[3:])
        os.rename(oldFname, newFname)
        write(newFname[:-4] + '_initial.png', read(newFname))


def scan_bond_only_scf(position):
    all_geoms_list = []
    scanAllDir = os.path.join(facetpath, 'scan_bond_all')
    scanDir = os.path.join(facetpath, 'scan_bond_' + position + '_only_scf')
    os.makedirs(scanDir, exist_ok=True)
    for geom in sorted(os.listdir(scanAllDir), key=str):
        if geom.endswith('.xyz'):
            all_geoms_list.append(geom)
    if position == 'top':
        count = 0
    elif position == 'bridge':
        count = 1
    elif position == 'hollow1':
        count = 2
    elif position == 'hollow2':
        count = 3
    selected_pos = []
    while count < len(all_geoms_list):
        selected_pos.append(all_geoms_list[count])
        count += 4
    for el in selected_pos:
        selected_src_xyz = os.path.join(scanAllDir, el)
        selected_src_png = os.path.join(selected_src_xyz[:-3] + 'png')
        shutil.copy2(selected_src_xyz, scanDir)
    scanGeomList = []
    for scanGeom in sorted(os.listdir(scanDir), key=str):
        if scanGeom.endswith('.xyz'):
            scanGeomList.append(scanGeom)
    for i, el in enumerate(scanGeomList):
        oldFname = os.path.join(scanDir, el)
        newFdir = os.path.join(scanDir, str(i).zfill(2))
        os.makedirs(newFdir, exist_ok=True)
        newFname = os.path.join(newFdir, '{}'.format(str(i).zfill(2)) + el[3:])
        os.rename(oldFname, newFname)
        write(newFname[:-4] + '_initial.png', read(newFname))


def set_up_scan_bond(path, pytemplate, avPath, species1, species2, repeats):
    with open(pytemplate, 'r') as f:
        pytemplate = f.read()

    for conf in os.listdir(path):
        confpath = os.path.join(path, conf)
        if os.path.isdir(confpath):
            for geom in os.listdir(confpath):
                if geom.endswith('.xyz'):
                    geomPath = os.path.join(confpath, geom)
                    O_index = get_index_adatom('O', geomPath)
                    H_index = get_index_adatom('H', geomPath)
                    Cu_index1 = get_index_surface_atom('O', geomPath)
                    Cu_index2 = get_index_surface_atom('H', geomPath)
                    # dist_Cu_O = get_bond_dist('O', geomPath)
                    # dist_Cu_H = get_bond_dist('H', geomPath)

                    prefix = geom.split('_')
                    calcDir = os.path.join(path, prefix[0])
                    os.makedirs(calcDir, exist_ok=True)

                    geomName = geom[:-4]

                    # List of bonds we want O-Cu; H-Cu
                    bonds = ((O_index, Cu_index1),
                             (H_index, Cu_index2))
                    avDist1 = get_av_dist(avPath, species1)
                    avDist2 = get_av_dist(avPath, species2)
                    avDists = (avDist1, avDist2)
                    trajPath = os.path.join(geom[:-4] + '.traj')
                    fname = os.path.join(calcDir, geom[:-4] + '.py')
                    with open(fname, 'w') as f:
                        f.write(pytemplate.format(geom=geom, bonds=bonds, avDists=avDists,
                                                  trajPath=trajPath, repeats=repeats, prefix=prefix[0], geomName=geomName))
                    f.close()
                f.close()


def set_up_scf_after_xtb(path, pytemplate):
    with open(pytemplate, 'r') as f:
        template = f.read()
        for conf in os.listdir(path):
            confPath = os.path.join(path, conf)
            if os.path.isdir(confPath):
                for geom in os.listdir(confPath):
                    if geom.endswith('final.xyz'):
                        # geomPath = os.path.join(confPath, geom)
                        fname = os.path.join(confPath, geom[:-9] + 'scf.py')
                        prefix = geom[:2]
                        rxn = geom[3:6]
                        with open(fname, 'w') as ff:
                            ff.write(template.format(prefix=prefix, rxn=rxn, TS_xyz=geom))
                        ff.close()

        f.close()


def set_up_scf_only(path, pytemplate):
    with open(pytemplate, 'r') as f:
        template = f.read()
        for conf in os.listdir(path):
            confPath = os.path.join(path, conf)
            if os.path.isdir(confPath):
                for geom in os.listdir(confPath):
                    if geom.endswith('.xyz'):
                        # geomPath = os.path.join(confPath, geom)
                        fname = os.path.join(confPath, geom[:-4] + '_scf.py')
                        prefix = geom[:2]
                        rxn = geom[3:6]
                        with open(fname, 'w') as ff:
                            ff.write(template.format(prefix=prefix, rxn=rxn, TS_xyz=geom))
                        ff.close()

        f.close()