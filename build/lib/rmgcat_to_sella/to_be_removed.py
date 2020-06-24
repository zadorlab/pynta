
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

