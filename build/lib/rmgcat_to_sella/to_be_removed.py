
'''TS'''
# def optimize(path):
#     for num, geom in enumerate(sorted(os.listdir(path), key=str)):
#         # for num, geom in os.listdir(path) if not endswith('.png'):
#         if geom.endswith('.xyz'):
#             print(geom)
#             job = read(os.path.join(path, geom))
#             calc = EMT()
#             job.set_calculator(calc)
#             trajPath = os.path.join(path, geom + '.traj')
#             opt = LBFGS(job, trajectory=trajPath)
#             # opt.calc(BFGS)
#             opt.run(fmax=0.05)

#             writeDir = os.path.join(saveDir, 'opt')
#             os.makedirs(writeDir, exist_ok=True)

#             write(os.path.join(writeDir, '{}'.format(
#                 str(num).zfill(3)) + 'result_relax.png'), read(trajPath))


# saveDir = '/Users/mgierad/00_SANDIA_WORK/05_rmgcat_to_stella/test/rmgcat_to_sella/Cu_111_tests/slab_optimized/Cu_111/neb_estimate/OH_O+H/'
# optimize(saveDir)
# geom = os.path.join(saveDir, '10_O+H.xyz')


# def optimize_fix_bond(geom):
#     job = read(geom)
#     trajPath = os.path.join(saveDir, '10.traj')

#     fix_O = []
#     fix_H = []

#     with open(geom, 'r') as f:
#         xyz_geom_file = f.readlines()
#         for num, line in enumerate(xyz_geom_file):
#             if 'H' in line:
#                 fix_H.append(num - 2)  # reading each line and adding only the one with tags = 1 (surface atoms). It is necessary to subtract 2 as indices in atom object starts with 0 and we also have to take into account that the first line in xyz file contains non xyz information
#             elif 'O' in line:
#                 fix_O.append(num - 2)
#     f.close()

#     # print(ads_O, ads_H)

#     calc = EMT()
#     job.set_calculator(calc)
#     c = FixBondLength(fix_H[0], fix_O[0])
#     job.set_constraint(c)
#     opt = LBFGS(job, trajectory=trajPath)
#     # opt.calc(BFGS)
#     opt.run(fmax=0.05)

#     writeDir = os.path.join(saveDir, 'opt')
#     os.makedirs(writeDir, exist_ok=True)

#     write(os.path.join(writeDir, 'result_fix_bond.png'), read(trajPath))

# optimize_fix_bond(geom)


# avDist_Cu_O = (1.75 + 1.90 + 1.90 + 1.90) / 4
# avDist_Cu_H = (1.53 + 1.73 + 1.75 + 1.73) / 4
# print('Average Cu-O distance: ', avDist_Cu_O)
# print('Average Cu-H distance: ', avDist_Cu_H)
# print()

# path_minima = os.path.join(facetpath, 'minima')

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
    # new checksym with globbing
    gd_ads_index = checkSymm(os.path.join(facetpath, TSdir))
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