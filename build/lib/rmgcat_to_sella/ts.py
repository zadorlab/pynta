from catkit import Gratoms
from catkit.gen.adsorption import AdsorptionSites, Builder
from catkit.build import molecule
from catkit.gen.molecules import get_3D_positions

from rmgcat_to_sella.adjacency_to_3d import get_edges, rmgcat_to_gratoms
from rmgcat_to_sella.find_all_nebs import get_all_species
from rmgcat_to_sella.graph_utils import node_test
from rmgcat_to_sella.main import CheckIfMinimasAlreadyCalculated

from ase.io import read, write
from ase import Atoms
from ase.optimize import LBFGS
from ase.calculators.emt import EMT
from ase.constraints import FixBondLength
from ase.utils.structure_comparator import SymmetryEquivalenceCheck

import numpy as np

from xtb import GFN0, GFN1

import itertools

from scipy.optimize import minimize

import yaml

import os

import shutil

from statistics import mean

from pathlib import Path

import networkx as nx

# def prepare_ts_estimate()
#     reactName = genTSestimate()
#     filtered_out_equiv_ts_estimate(saveDir, reactName)

def genTSestimate(slab, repeats, yamlfile, facetpath, rotAngle, scfactor):
    # Caclulate how many distinct rotations exists based on the slab symmetry
    # rotAngle = slab.get_cell_lengths_and_angles()[5]
    # possibleRotations = (360 / slab.get_cell_lengths_and_angles()[5]) - 1

    possibleRotations = (360 / rotAngle) - 1

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

    r_name = '+'.join([str(species.symbols) for species in reactants])
    p_name = '+'.join([str(species.symbols) for species in products])

    # rxn_name = r_name + '_' + p_name

    unique_species = []
    unique_bonds = []
    images = []
    # check if any products are the same as any reactants
    for species1, bond in zip(speciesInd, bonds):
        for species2 in unique_species:
            if nx.is_isomorphic(species1.graph, species2.graph, node_test):
                break
        else:
            images.append(get_3D_positions(species1))
            unique_species.append(species1)
            unique_bonds.append(bond)
    # print(unique_species)

    r_name_list = [str(species.symbols) for species in reactants]
    p_name_list = [str(species.symbols) for species in products]

    if len(r_name_list) < len(p_name_list):
        print('Reactant structure will be used to estimate TS')
        rpDir = 'reactants'
    else:
        print('Products will be used to estimate TS')
        rpDir = 'products'

    saveDir = os.path.join(facetpath, 'TS_estimate')
    if os.path.exists(saveDir):
        shutil.rmtree(saveDir)
        os.makedirs(saveDir)
    else:
        os.makedirs(saveDir)

    # ADSORBATES
    '''
    This part here is a bit hard coded, especially when dealing with reactants with >= 3 atoms. The code works for couple of species included in if else statement below. Probably will fail for other species. This is becouse ase.build.molecule function numbers atoms very diferently depending on the molecule specify. To be resolved somehow.
    '''
    atom1 = 0  # atom1 should always have index 0, i.e. the first atom
    atom2 = 1  # well, here is the problem becouse ase.build.molecule counts atoms differently than my input in the yaml file. Probably the problem with the way yamlfile is converted to species - to be resolved
    if len(r_name_list) <= len(p_name_list):
        TS_candidate = molecule(r_name_list[0])[0]
        reactName = p_name
    else:
        # TS_candidate = molecule(p_name_list[0])[0]
        '''images[2] keeps info about product. It's a Gratom object. Cannot use catkit's molecule method becouse it generates gemetry based on wierd order chemical formula string. Sometimes it works to use molecule method but in general it is better to use get_3D_positions() as it converts yaml info into 3d position. In images I keep the oryginal order of elements with connectivity info. In TS_candidate.get_chemical_formula() elements are sorted and grouped '''
        TS_candidate = images[2] 
        reactName = r_name

    ''' Different cases '''
    if TS_candidate.get_chemical_formula() == 'CHO':
        atom2 = 2
        bondedThrough = [0]
    elif TS_candidate.get_chemical_formula() == 'COH':
        atom2 = 2
        bondedThrough = [0]
    elif TS_candidate.get_chemical_formula() == 'CHO2':
        # atom1 = 0
        # atom2 = 1
        bondedThrough = [2]  # connect through oxygen

    elif TS_candidate.get_chemical_formula() == 'C2H5O2':
        atom0 = 0
        atom1 = 1
        atom2 = 5
        bondedThrough = [6]  # connect through oxygen
    else:
        # atom2 = 1
        bondedThrough = [0]
        
    bondlen = TS_candidate.get_distance(atom1, atom2)
    # bondlen = TS_candidate.get_distance(atom1, atom2)
    # angle = TS_candidate.get_angle(0, 1, 5)

    # TS_candidate.set_distance(atom1, atom2, blen * scfactor, fix=0)

    '''Final ridgid rotations to orientade TS_candidate on the surface'''
    if len(TS_candidate.get_tags()) < 3:
        TS_candidate.rotate(90, 'y')
        TS_candidate.set_distance(atom1, atom2, bondlen * scfactor, fix=0)
    elif len(TS_candidate.get_tags()) == 3:
        TS_candidate.rotate(90, 'z')
        TS_candidate.set_distance(atom1, atom2, bondlen * scfactor, fix=0)
    # else:
    elif TS_candidate.get_chemical_formula() == 'C2H5O2':
        # TS_candidate.rotate(60, 'y')
        TS_candidate.rotate(90, 'z')
        ''' Should work for H2CO*OCH3, i.e. COH3+HCOH '''
        TS_candidate.set_angle(atom2, atom1, atom0, -45, indices=[0, 1, 2, 3, 4], add=True)
        TS_candidate.set_distance(atom1, atom2, bondlen * scfactor, fix=1, indices=[0, 1, 2, 3, 4])
        # indices=[0, 1, 2, 3, 4]
    elif TS_candidate.get_chemical_formula() == 'CHO2':
        TS_candidate.rotate(90, 'z')
        TS_candidate.set_distance(atom1, atom2, bondlen * scfactor, fix=0)


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

    while count <= possibleRotations:
        structs = ads_builder.add_adsorbate(
            TS_candidate, bondedThrough, -1, auto_construct=False)  # change to True will make bondedThrough work. Now it uses TS_candidate,rotate... to generate adsorbed strucutres
        big_slab = slab * repeats
        nslab = len(slab)

        for i, struc in enumerate(structs):
            big_slab_ads = big_slab + struc[nslab:]
            write(os.path.join(saveDir, '{}'.format(
                str(i + len(structs) * count).zfill(3)) + '_' + reactName + '.xyz'), big_slab_ads)
        
        TS_candidate.rotate(rotAngle, 'z')
        count += 1
    # return reactName

    '''Filtering out symmetry equivalent sites '''
# def filtered_out_equiv_ts_estimate(saveDir, reactName)
    filtered_equivalen_sites = checkSymmBeforeXTB(saveDir)
    for eqsites in filtered_equivalen_sites:
        try:
            fileToRemove = os.path.join(saveDir, eqsites + '_' + reactName + '.xyz')
            os.remove(fileToRemove)
        except OSError:
            print("Error while deleting file : ", fileToRemove)

    for prefix, noneqsites in enumerate(sorted(os.listdir(saveDir))):
        prefix = str(prefix).zfill(3)
        oldfname = os.path.join(saveDir, noneqsites)
        newfname = os.path.join(saveDir, prefix + noneqsites[3:])
        os.rename(oldfname, newfname)

def gen_xyz_from_traj(avDistPath, species):
    # if species == 'C':
    #     speciesPath = os.path.join(avDistPath, species + 'O')
    # else:
    speciesPath = os.path.join(avDistPath, species)
    for traj in sorted(os.listdir(speciesPath), key=str):
        if traj.endswith('.traj'):
            srcTrajPath = os.path.join(speciesPath, traj)
            desTrajPath = os.path.join(speciesPath, traj[:-5] + '_final.xyz')
            write(desTrajPath, read(srcTrajPath))


def get_av_dist(avDistPath, species, scfactor_surface, scaled=False):
    # nslab = len(read(slab) * repeats)
    surface_atoms_indices = []
    adsorbate_atoms_indices = []
    all_conf_dist = []
    gen_xyz_from_traj(avDistPath, species)
    speciesPath = os.path.join(avDistPath, species)
    if species in ['CH3O', 'CH2O']:
        species = 'O'
    if len(species) > 1:
        species = species[:1]
    
    # print(species)
    for xyz in sorted(os.listdir(speciesPath), key=str):
        if xyz.endswith('_final.xyz'):
            # if xyz.endswith('.traj'):
            xyzPath = os.path.join(speciesPath, xyz)
            # trajPath = os.path.join(speciesPath, xyz[:-10] + '.traj')
            # print(trajPath)
            conf = read(xyzPath)
            # print(conf)
        # if xyz.endswith('*traj'):
            with open(xyzPath, 'r') as f:
                xyzFile = f.readlines()
                for num, line in enumerate(xyzFile):
                    '''For Cu(111) we can put ' 1 ' instead of 'Cu ' to limit calculations to surface Cu atoms only. For Cu(211) ase is not generating tags like this, so full calculation have to be performed, i.e. all surface atoms'''
                    if 'Cu ' in line:
                        surface_atoms_indices.append(num - 2)
                    elif species in line and not 'Cu' in line:
                        ''' We need to have additional statement 'not 'Cu' in line' because for C it does not work without it'''
                        # if not 'Cu' in line:
                        adsorbate_atoms_indices.append(num - 2)
            f.close()
            dist = float(min(conf.get_distances(
                adsorbate_atoms_indices[0], surface_atoms_indices)))
            all_conf_dist.append(dist)
            if scaled:
                meanDist = mean(all_conf_dist) * scfactor_surface
            else:
                meanDist = mean(all_conf_dist)
    # print(meanDist)
    return meanDist


# specify adsorbate atom symbol and bond distance with the closest surface
# metal atom will be calculated
def get_bond_dist(ads_atom, geom):
    surface_atom = []
    adsorbate_atom = []
    struc = read(geom)
    if len(ads_atom) > 1:
        ads_atom = ads_atom[:-1]
    # if ads_atom == 'C':
    #     ads_atom == 'CO'
    with open(geom, 'r') as f:
        xyz_geom_file = f.readlines()
        for num, line in enumerate(xyz_geom_file):
            if ' 1 ' in line:
                # reading each line and adding only the one with tags = 1
                # (surface atoms). It is necessary to subtract 2 as indices in
                # atom object starts with 0 and we also have to take into
                # account that the first line in xyz file contains non xyz
                # information
                surface_atom.append(num - 2)
            elif ads_atom in line:
                if not "Cu" in line:
                    adsorbate_atom.append(num - 2)
    f.close()

    if len(adsorbate_atom) > 1:
        # print(adsorbate_atom[0])
        dist_Cu_adsorbate = min(struc.get_distances(
            adsorbate_atom[0], surface_atom))
        return dist_Cu_adsorbate
    else:
        dist_Cu_adsorbate = min(
            struc.get_distances(adsorbate_atom, surface_atom))
        return dist_Cu_adsorbate


def get_index_adatom(ads_atom, geom):
    adsorbate_atom = []
    if len(ads_atom) > 1:
        ads_atom = ads_atom[:1]
    with open(geom, 'r') as f:
        xyz_geom_file = f.readlines()
        for num, line in enumerate(xyz_geom_file):
            if ads_atom in line:
                if not 'Cu' in line:
                    adsorbate_atom.append(num - 2)
    f.close()
    return adsorbate_atom[0]


def get_index_surface_atom(ads_atom, geom):
    surface_atom = []
    adsorbate_atom = []
    if len(ads_atom) > 1:
        ads_atom = ads_atom[:1]
    # print(ads_atom)
    struc = read(geom)
    with open(geom, 'r') as f:
        xyz_geom_file = f.readlines()
        for num, line in enumerate(xyz_geom_file):
            if 'Cu ' in line:
                surface_atom.append(num - 2)
            elif ads_atom in line:
                if not 'Cu' in line:
                    adsorbate_atom.append(num - 2)
    f.close()

    all_dist_surface_adsorbate = struc.get_distances(
        adsorbate_atom[0], surface_atom)
    min_dist_surface_adsorbate = min(
        struc.get_distances(adsorbate_atom[0], surface_atom))
    # get index of the surface atom for which distance to adsorbate is the
    # lowest
    index = np.where(all_dist_surface_adsorbate == min_dist_surface_adsorbate)

    return surface_atom[index[0][0]]


'''
def set_adsorbate_positions(x):
    # make a copy of adsorbate
    ads = TS_candidate.copy()
    # extract center and axis/angle
    center = x[:3]
    axis = x[3:]
    # angle is encoded as the norm of the axis
    # ASE expects degrees, but we want to use radians so that the center and
    # terms have about the same magnitude
    angle = np.linalg.norm(axis) * 180 / np.pi
    # Shift the adsorbate center to "center"
    ads.positions += center - ads.positions.mean(0)
    # Rotate the adsorbate by angle "angle" through axis "axis"
    ads.rotate(angle, v=axis, center=center)
    return ads


def penalty(x):
    # Temporarily combine slab and rotated/translated adsorbate
    atoms = bigSlab + set_adsorbate_positions(x)
    calc = EMT()
    # calc = GPAW(xc='PBE', mode = 'pw')
    # c = FixBondLength(36, 37)
    atoms.set_calculator(calc)
    # atoms.set_constraint(c)
    # print(atoms.get_potential_energy())
    energy = 0
    for bond, dist in zip(bonds, dists):
        energy += (atoms.get_distance(*bond) - dist)**2 + \
            atoms.get_potential_energy()
    # with open ('results.log', 'a+') as f:
    #     print(energy, '\t', atoms.get_potential_energy(), file=f)
    # f.close()
    print(energy, '\t', atoms.get_potential_energy())
    return energy


def optimizePenalty(path, geom):
    x0 = np.zeros(6)
    x0[:3] = TS_candidate.positions.mean(0)
    x0[3:] = [0, np.pi / 2, 0]  # just an arbitrary initial guess
    res = minimize(penalty, x0, tol=1e-03)
    print(res)
    ads_opt = set_adsorbate_positions(res['x'])
    # view(slab + ads_opt)
    # write('maciek_emt.png', bigSlab + ads_opt, rotation='10z,-80x')

    # def unique_file(basename, ext):
    #     actualname = "%s.%s" % (basename, ext)
    #     c = itertools.count()
    #     while os.path.exists(actualname):
    #         actualname = "%s (%d).%s" % (basename, next(c), ext)
    #     return actualname

    save_png = os.path.join(path, geom + '.png')
    save_xyz = os.path.join(path, geom + '.xyz')
    write(save_png, bigSlab + ads_opt)
    write(save_xyz, bigSlab + ads_opt)
    # write('maciek_emt_2.xyz', bigSlab + ads_opt)
    # write('maciek_emt_3.png', bigSlab + ads_opt, rotation='10z,-80x')


def run_preopt():
    for geom in sorted(os.listdir(saveDir), key=str):
        print('E_pot + penalty', '\t', 'E_pot')

        geomPath = os.path.join(saveDir, geom)
        SaveDir_final = os.path.join(saveDir, 'neb_candidate_xtb')
        os.makedirs(SaveDir_final, exist_ok=True)
        if geom.endswith('.xyz'):
            O_index = get_index_adatom('O', geomPath)
            H_index = get_index_adatom('H', geomPath)
            Cu_index1 = get_index_surface_atom('O', geomPath)
            Cu_index2 = get_index_surface_atom('H', geomPath)
            dist_Cu_O = get_bond_dist('O', geomPath)
            dist_Cu_H = get_bond_dist('H', geomPath)

            adsorbed = read(geomPath)
            bigSlab = slab * repeats
            nbigSlab = len(bigSlab)
            TS_candidate = adsorbed[nbigSlab:]
        # write('maciek_emt_1.png', bigSlab + TS_candidate)

        # List of bonds we want O-Cu; H-Cu
            bonds = ((O_index, Cu_index1),
                     (H_index, Cu_index2))
        # For now, assume the optimal Cu-C and Cu-O bond distances are both 2.0
        # Ang
            avDist1 = get_av_dist('H', 'Cu_111')
            avDist2 = get_av_dist('O', 'Cu_111')
            dists = (avDist1, avDist2)
            geom = geom[:-4]
            # print(bonds)
            optimizePenalty(SaveDir_final, geom)
'''


# def checkSymm(path, TSdir):
#     good_adsorbate = []
#     result_list = []
#     gpath = os.path.join(path, TSdir)
#     for geom in sorted(os.listdir(gpath), key=str):
#         geomDir = os.path.join(gpath, geom)
#         if os.path.isdir(geomDir):
#             for traj in os.listdir(geomDir):
#                 if traj.endswith('.traj'):
#                     adsorbed = read(os.path.join(geomDir, traj))
#                     adsorbed.pbc = True
#                     comparator = SymmetryEquivalenceCheck()
#                     result = comparator.compare(adsorbed, good_adsorbate)
#                     result_list.append(result)
#                     if result is False:
#                         good_adsorbate.append(adsorbed)
#     unique_index = []
#     for num, res in enumerate(result_list):
#         if res is False:
#             unique_index.append(str(num).zfill(3))
#     return unique_index


def checkSymm(path):
    good_adsorbate = []
    result_list = []
    geomlist = sorted(Path(path).glob('**/*.traj'))
    for geom in geomlist:
        adsorbed = read(geom)
        adsorbed.pbc = True
        comparator = SymmetryEquivalenceCheck()
        result = comparator.compare(adsorbed, good_adsorbate)
        result_list.append(result)
        if result is False:
            good_adsorbate.append(adsorbed)
    unique_index = []
    for num, res in enumerate(result_list):
        if res is False:
            unique_index.append(str(num).zfill(3))
    return unique_index


def checkSymmBeforeXTB(path):
    good_adsorbate = []
    result_list = []
    geomlist = sorted(Path(path).glob('*.xyz'))
    for geom in geomlist:
        adsorbed = read(geom)
        adsorbed.pbc = True
        comparator = SymmetryEquivalenceCheck()
        result = comparator.compare(adsorbed, good_adsorbate)
        result_list.append(result)
        if result is False:
            good_adsorbate.append(adsorbed)
    notunique_index = []
    for num, res in enumerate(result_list):
        if res is True:
            ''' Better to have all symmetry equivalent site here in a list. The workflow will remove them in getTSestimate function keeping all symmetry distinct sites '''
            notunique_index.append(str(num).zfill(3))
    return notunique_index

    # for geom in sorted(os.listdir(path), key=str):
    #     geomDir = os.path.join(path, geom)
    #     for traj in os.listdir(geomDir):
    #         if traj.endswith('.xyz'):
    #             adsorbed = read(os.path.join(geomDir, traj))
    #             adsorbed.pbc = True
    #             comparator = SymmetryEquivalenceCheck()
    #             result = comparator.compare(adsorbed, good_adsorbate)
    #             result_list.append(result)
    #             if result is False:
    #                 good_adsorbate.append(adsorbed)
    # unique_index = []
    # for num, res in enumerate(result_list):
    #     if res is False:
    #         unique_index.append(str(num).zfill(3))
    # return unique_index


def create_unique_TS(facetpath, TSdir):
    # checkSymmPath = os.path.join(path, 'TS_estimate')

    # gd_ads_index = checkSymm(facetpath, TSdir)
    gd_ads_index = checkSymm(os.path.join(facetpath, TSdir)) # new checksym with globbing
    for i, index in enumerate(gd_ads_index):
        uniqueTSdir = os.path.join(
            facetpath, TSdir + '_unique', str(i).zfill(2))
        if os.path.isdir(uniqueTSdir):
            shutil.rmtree(uniqueTSdir)
            os.makedirs(uniqueTSdir, exist_ok=True)
        else:
            os.makedirs(uniqueTSdir, exist_ok=True)
        gpath = os.path.join(facetpath, TSdir)
        for geom in os.listdir(gpath):
            geomDir = os.path.join(gpath, geom)
            if os.path.isdir(geomDir):
                for traj in sorted(os.listdir(geomDir)):
                    if traj.startswith(gd_ads_index[i]) and traj.endswith('.traj'):
                        srcFile = os.path.join(geomDir, traj)
                        destFile = os.path.join(uniqueTSdir, traj[:-5] + '_ts')
                        write(destFile + '.xyz', read(srcFile))
                        write(destFile + '.png', read(srcFile))
        #         shutil.copy2(os.path.join(path, geom), uniqueTSdir)
        for ts in os.listdir(uniqueTSdir):
            # if neb.endswith('.xyz'):
            TS_xyz_Dir = os.path.join(uniqueTSdir, ts)
            newTS_xyz_Dir = os.path.join(
                uniqueTSdir, str(i).zfill(2) + ts[3:])
            # NEB_png_Dir = os.path.join(uniqueTSdir, str(
            #     i).zfill(3) + neb[2:][:-4] + '.png')
            os.rename(TS_xyz_Dir, newTS_xyz_Dir)
            # write(NEB_png_Dir, read(newNEB_xyz_Dir))


# uniqueNEB = '/Users/mgierad/00_SANDIA_WORK/05_rmgcat_to_stella/test/rmgcat_to_sella/Cu_111_tests/slab_optimized/Cu_111/neb_estimate/OH_O+H/neb_candidate_unique/'

# pytemplate = '/Users/mgierad/00_SANDIA_WORK/05_rmgcat_to_stella/test/rmgcat_to_sella/Cu_111_tests/slab_optimized/pytemplate_neb.py'


def create_all_TS(facetpath):
    all_TS_candidate_path = os.path.join(facetpath, 'TS_candidate')
    os.makedirs(all_TS_candidate_path, exist_ok=True)
    gpath = os.path.join(facetpath, 'TS_estimate')
    for geom in os.listdir(gpath):
        geomDir = os.path.join(gpath, geom)
        if os.path.isdir(geomDir):
            for traj in sorted(os.listdir(geomDir)):
                if traj.endswith('.traj'):
                    # print(str(num/2).zfill(3))
                    xyz_dir_path = os.path.join(
                        all_TS_candidate_path, traj[:3])
                    os.makedirs(xyz_dir_path, exist_ok=True)
                    srcFile = os.path.join(geomDir, traj)
                    destFile = os.path.join(xyz_dir_path, traj[:-5])
                    write(destFile + '.xyz', read(srcFile))
                    write(destFile + '_initial.png', read(srcFile))


def create_all_TS_job_files(facetpath, pytemplate):
    all_TS_candidate_path = os.path.join(facetpath, 'TS_candidate')
    with open(pytemplate, 'r') as f:
        pytemplate = f.read()

    for struc in os.listdir(all_TS_candidate_path):
        TSdir = os.path.join(all_TS_candidate_path, struc)
        if os.path.isdir(TSdir):
            TSpath = os.path.join(all_TS_candidate_path, struc)
            for file in os.listdir(TSpath):
                if file.endswith('.xyz'):
                    fname = os.path.join(
                        all_TS_candidate_path, f'{struc}_{file[:-4][4:]}_relax.py')
                    with open(fname, 'w') as f:
                        f.write(pytemplate.format(TS=os.path.join(
                            struc, file), rxn=file[:-4][4:], prefix=file[:3]))
                    f.close()
    f.close()


def create_TS_unique_job_files(facetpath, TSdir, pytemplate):
    unique_TS_candidate_path = os.path.join(facetpath, TSdir + '_unique')
    with open(pytemplate, 'r') as f:
        pytemplate = f.read()

    for struc in os.listdir(unique_TS_candidate_path):
        TSdir = os.path.join(unique_TS_candidate_path, struc)
        if os.path.isdir(TSdir):
            TSpath = os.path.join(unique_TS_candidate_path, struc)
            for fl in os.listdir(TSpath):
                if fl.endswith('.xyz'):
                    fname = os.path.join(
                        unique_TS_candidate_path, fl[:-4] + '.py')
                    with open(fname, 'w') as f:
                        f.write(pytemplate.format(TS=os.path.join(
                            struc, fl), rxn=fl[3:-7], prefix=fl[:2]))
                    f.close()
    f.close()


def set_up_penalty_xtb(path, pytemplate, repeats, slabopt, species_list, scfactor_surface, scaled1, scaled2):
    '''Species3 is/should be optional '''
    # print(species3)
    fPath = os.path.split(path)
    avPath = os.path.join(fPath[0], 'minima')
    avDist1 = get_av_dist(avPath, species_list[0], scfactor_surface, scaled1)
    avDist2 = get_av_dist(avPath, species_list[1], scfactor_surface, scaled2)
    '''the code belowe does not work in the loop, probably two differet avPath variable needed to accouut for the poscible scenaario that one species was already calculated (other set of calculations - different reactions, whereas the second in calculated here for the first time) '''
    # checkMinimaPath = os.path.dirname(os.getcwd())
    # spList = [species1, species2]

    # for species in spList:
    #     isItCalculated = CheckIfMinimasAlreadyCalculated(
    #         checkMinimaPath, species)
    #     if isItCalculated is False:
    #         fPath = os.path.split(path)
    #         avPath = os.path.join(fPath[0], 'minima')
    #     else:
    #         avPath = isItCalculated[1]

    # species_list = [species1, species2, species3]

    with open(pytemplate, 'r') as f:
        pytemplate = f.read()
    
    for geom in sorted(os.listdir(path)):
        if geom.endswith('.xyz'):
            bonds = []
            geomPath = os.path.join(path, geom)
            
            for species in species_list:
                sp_index = get_index_adatom(species, geomPath)
                Cu_index = get_index_surface_atom(species, geomPath)
                # print(sp_index, Cu_index)
                bonds.append((sp_index, Cu_index))

            # sp1_index = get_index_adatom(species1, geomPath)
            # sp2_index = get_index_adatom(species2, geomPath)
            # # sp3_index = get_index_adatom(species3, geomPath)
            # Cu_index1 = get_index_surface_atom(species1, geomPath)
            # Cu_index2 = get_index_surface_atom(species2, geomPath)
            # # Cu_index3 = get_index_surface_atom(species3, geomPath)

            prefix = geom.split('_')
            calcDir = os.path.join(path, prefix[0])
            os.makedirs(calcDir, exist_ok=True)

            geomName = geom[:-4]

            # List of bonds we want O-Cu; H-Cu
            # bonds = ((sp1_index, Cu_index1),
            #          (sp2_index, Cu_index2))
            # print(type(bonds))
            avDists = (avDist1, avDist2)
            trajPath = os.path.join(geom[:-4] + '.traj')
            init_png = os.path.join(calcDir, geom[:-4] + '_initial.png')
            write(init_png, read(geomPath))
            fname = os.path.join(calcDir, geom[:-4] + '.py')
            with open(fname, 'w') as f:
                f.write(pytemplate.format(geom=geom, bonds=bonds,
                                          avDists=avDists, trajPath=trajPath,
                                          repeats=repeats, prefix=prefix[0],
                                          geomName=geomName, slabopt=slabopt))
            f.close()
            shutil.move(geomPath, calcDir)
    f.close()

    try:
        rmpath = os.path.join(path, 'initial_png')
        shutil.rmtree(rmpath)
    except FileNotFoundError:
        pass
        # print('No files to delete')

    


def copyMinimasPrevCalculated(checkMinimaDir, sp1, sp2, dstDir, slabname):
    '''
    Copy directories
    '''
    speciesList = [sp1, sp2]
    for species in speciesList:
        isItCalculated = CheckIfMinimasAlreadyCalculated(
            checkMinimaDir, species, slabname)
        if isItCalculated is False:
            pass
        else:
            try:
                copyPath = isItCalculated[1]
                dstDir = os.path.join(dstDir, species)
                shutil.copytree(copyPath, dstDir)
                dstDir = os.path.split(dstDir)[0]
            except FileExistsError:
                dstDir = os.path.split(dstDir)[0]
                print('All required files already exist.')
                pass
