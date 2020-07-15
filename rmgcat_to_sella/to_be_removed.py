
# '''TS'''
# # def optimize(path):
# #     for num, geom in enumerate(sorted(os.listdir(path), key=str)):
# #         # for num, geom in os.listdir(path) if not endswith('.png'):
# #         if geom.endswith('.xyz'):
# #             print(geom)
# #             job = read(os.path.join(path, geom))
# #             calc = EMT()
# #             job.set_calculator(calc)
# #             trajPath = os.path.join(path, geom + '.traj')
# #             opt = LBFGS(job, trajectory=trajPath)
# #             # opt.calc(BFGS)
# #             opt.run(fmax=0.05)

# #             writeDir = os.path.join(saveDir, 'opt')
# #             os.makedirs(writeDir, exist_ok=True)

# #             write(os.path.join(writeDir, '{}'.format(
# #                 str(num).zfill(3)) + 'result_relax.png'), read(trajPath))


# # saveDir = '/Users/mgierad/00_SANDIA_WORK/05_rmgcat_to_stella/test/rmgcat_to_sella/Cu_111_tests/slab_optimized/Cu_111/neb_estimate/OH_O+H/'
# # optimize(saveDir)
# # geom = os.path.join(saveDir, '10_O+H.xyz')


# # def optimize_fix_bond(geom):
# #     job = read(geom)
# #     trajPath = os.path.join(saveDir, '10.traj')

# #     fix_O = []
# #     fix_H = []

# #     with open(geom, 'r') as f:
# #         xyz_geom_file = f.readlines()
# #         for num, line in enumerate(xyz_geom_file):
# #             if 'H' in line:
# #                 fix_H.append(num - 2)  # reading each line and adding only the one with tags = 1 (surface atoms). It is necessary to subtract 2 as indices in atom object starts with 0 and we also have to take into account that the first line in xyz file contains non xyz information
# #             elif 'O' in line:
# #                 fix_O.append(num - 2)
# #     f.close()

# #     # print(ads_O, ads_H)

# #     calc = EMT()
# #     job.set_calculator(calc)
# #     c = FixBondLength(fix_H[0], fix_O[0])
# #     job.set_constraint(c)
# #     opt = LBFGS(job, trajectory=trajPath)
# #     # opt.calc(BFGS)
# #     opt.run(fmax=0.05)

# #     writeDir = os.path.join(saveDir, 'opt')
# #     os.makedirs(writeDir, exist_ok=True)

# #     write(os.path.join(writeDir, 'result_fix_bond.png'), read(trajPath))

# # optimize_fix_bond(geom)


# # avDist_Cu_O = (1.75 + 1.90 + 1.90 + 1.90) / 4
# # avDist_Cu_H = (1.53 + 1.73 + 1.75 + 1.73) / 4
# # print('Average Cu-O distance: ', avDist_Cu_O)
# # print('Average Cu-H distance: ', avDist_Cu_H)
# # print()

# # path_minima = os.path.join(facetpath, 'minima')

# '''
# def set_adsorbate_positions(x):
#     # make a copy of adsorbate
#     ads = TS_candidate.copy()
#     # extract center and axis/angle
#     center = x[:3]
#     axis = x[3:]
#     # angle is encoded as the norm of the axis
#     # ASE expects degrees, but we want to use radians so that the center and
#     # terms have about the same magnitude
#     angle = np.linalg.norm(axis) * 180 / np.pi
#     # Shift the adsorbate center to "center"
#     ads.positions += center - ads.positions.mean(0)
#     # Rotate the adsorbate by angle "angle" through axis "axis"
#     ads.rotate(angle, v=axis, center=center)
#     return ads


# def penalty(x):
#     # Temporarily combine slab and rotated/translated adsorbate
#     atoms = bigSlab + set_adsorbate_positions(x)
#     calc = EMT()
#     # calc = GPAW(xc='PBE', mode = 'pw')
#     # c = FixBondLength(36, 37)
#     atoms.set_calculator(calc)
#     # atoms.set_constraint(c)
#     # print(atoms.get_potential_energy())
#     energy = 0
#     for bond, dist in zip(bonds, dists):
#         energy += (atoms.get_distance(*bond) - dist)**2 + \
#             atoms.get_potential_energy()
#     # with open ('results.log', 'a+') as f:
#     #     print(energy, '\t', atoms.get_potential_energy(), file=f)
#     # f.close()
#     print(energy, '\t', atoms.get_potential_energy())
#     return energy


# def optimizePenalty(path, geom):
#     x0 = np.zeros(6)
#     x0[:3] = TS_candidate.positions.mean(0)
#     x0[3:] = [0, np.pi / 2, 0]  # just an arbitrary initial guess
#     res = minimize(penalty, x0, tol=1e-03)
#     print(res)
#     ads_opt = set_adsorbate_positions(res['x'])
#     # view(slab + ads_opt)
#     # write('maciek_emt.png', bigSlab + ads_opt, rotation='10z,-80x')

#     # def unique_file(basename, ext):
#     #     actualname = "%s.%s" % (basename, ext)
#     #     c = itertools.count()
#     #     while os.path.exists(actualname):
#     #         actualname = "%s (%d).%s" % (basename, next(c), ext)
#     #     return actualname

#     save_png = os.path.join(path, geom + '.png')
#     save_xyz = os.path.join(path, geom + '.xyz')
#     write(save_png, bigSlab + ads_opt)
#     write(save_xyz, bigSlab + ads_opt)
#     # write('maciek_emt_2.xyz', bigSlab + ads_opt)
#     # write('maciek_emt_3.png', bigSlab + ads_opt, rotation='10z,-80x')


# def run_preopt():
#     for geom in sorted(os.listdir(saveDir), key=str):
#         print('E_pot + penalty', '\t', 'E_pot')

#         geomPath = os.path.join(saveDir, geom)
#         SaveDir_final = os.path.join(saveDir, 'neb_candidate_xtb')
#         os.makedirs(SaveDir_final, exist_ok=True)
#         if geom.endswith('.xyz'):
#             O_index = get_index_adatom('O', geomPath)
#             H_index = get_index_adatom('H', geomPath)
#             Cu_index1 = get_index_surface_atom('O', geomPath)
#             Cu_index2 = get_index_surface_atom('H', geomPath)
#             dist_Cu_O = get_bond_dist('O', geomPath)
#             dist_Cu_H = get_bond_dist('H', geomPath)

#             adsorbed = read(geomPath)
#             bigSlab = slab * repeats
#             nbigSlab = len(bigSlab)
#             TS_candidate = adsorbed[nbigSlab:]
#         # write('maciek_emt_1.png', bigSlab + TS_candidate)

#         # List of bonds we want O-Cu; H-Cu
#             bonds = ((O_index, Cu_index1),
#                      (H_index, Cu_index2))
#         # For now, assume the optimal Cu-C and Cu-O bond distances are both 2.0
#         # Ang
#             avDist1 = get_av_dist('H', 'Cu_111')
#             avDist2 = get_av_dist('O', 'Cu_111')
#             dists = (avDist1, avDist2)
#             geom = geom[:-4]
#             # print(bonds)
#             optimizePenalty(SaveDir_final, geom)
# '''




# # def checkSymm(path, TSdir):
# #     good_adsorbate = []
# #     result_list = []
# #     gpath = os.path.join(path, TSdir)
# #     for geom in sorted(os.listdir(gpath), key=str):
# #         geomDir = os.path.join(gpath, geom)
# #         if os.path.isdir(geomDir):
# #             for traj in os.listdir(geomDir):
# #                 if traj.endswith('.traj'):
# #                     adsorbed = read(os.path.join(geomDir, traj))
# #                     adsorbed.pbc = True
# #                     comparator = SymmetryEquivalenceCheck()
# #                     result = comparator.compare(adsorbed, good_adsorbate)
# #                     result_list.append(result)
# #                     if result is False:
# #                         good_adsorbate.append(adsorbed)
# #     unique_index = []
# #     for num, res in enumerate(result_list):
# #         if res is False:
# #             unique_index.append(str(num).zfill(3))
# #     return unique_index


#     # for geom in sorted(os.listdir(path), key=str):
#     #     geomDir = os.path.join(path, geom)
#     #     for traj in os.listdir(geomDir):
#     #         if traj.endswith('.xyz'):
#     #             adsorbed = read(os.path.join(geomDir, traj))
#     #             adsorbed.pbc = True
#     #             comparator = SymmetryEquivalenceCheck()
#     #             result = comparator.compare(adsorbed, good_adsorbate)
#     #             result_list.append(result)
#     #             if result is False:
#     #                 good_adsorbate.append(adsorbed)
#     # unique_index = []
#     # for num, res in enumerate(result_list):
#     #     if res is False:
#     #         unique_index.append(str(num).zfill(3))
#     # return unique_index

# def create_unique_TS(facetpath, TSdir):
#     # checkSymmPath = os.path.join(path, 'TS_estimate')

#     # gd_ads_index = checkSymm(facetpath, TSdir)
#     # new checksym with globbing
#     gd_ads_index = checkSymm(os.path.join(facetpath, TSdir))
#     for i, index in enumerate(gd_ads_index):
#         uniqueTSdir = os.path.join(
#             facetpath, TSdir + '_unique', str(i).zfill(2))
#         if os.path.isdir(uniqueTSdir):
#             shutil.rmtree(uniqueTSdir)
#             os.makedirs(uniqueTSdir, exist_ok=True)
#         else:
#             os.makedirs(uniqueTSdir, exist_ok=True)
#         gpath = os.path.join(facetpath, TSdir)
#         for geom in os.listdir(gpath):
#             geomDir = os.path.join(gpath, geom)
#             if os.path.isdir(geomDir):
#                 for traj in sorted(os.listdir(geomDir)):
#                     if traj.startswith(gd_ads_index[i]) and traj.endswith('.traj'):
#                         srcFile = os.path.join(geomDir, traj)
#                         destFile = os.path.join(uniqueTSdir, traj[:-5] + '_ts')
#                         write(destFile + '.xyz', read(srcFile))
#                         write(destFile + '.png', read(srcFile))
#         #         shutil.copy2(os.path.join(path, geom), uniqueTSdir)
#         for ts in os.listdir(uniqueTSdir):
#             # if neb.endswith('.xyz'):
#             TS_xyz_Dir = os.path.join(uniqueTSdir, ts)
#             newTS_xyz_Dir = os.path.join(
#                 uniqueTSdir, str(i).zfill(2) + ts[3:])
#             # NEB_png_Dir = os.path.join(uniqueTSdir, str(
#             #     i).zfill(3) + neb[2:][:-4] + '.png')
#             os.rename(TS_xyz_Dir, newTS_xyz_Dir)
#             # write(NEB_png_Dir, read(newNEB_xyz_Dir))


# # def getBarrierAfterIRC(path, distCutoff):
#         # barrier is calculate in referece to the forward or reverse trajectory as
#         # after the irc calculations. Note, that some forward or reverse
#         # strucutres are not converge to the local minima - further optimization
#         # required.
# #     confList = []
# #     for conf in sorted(os.listdir(path), key=str):
# #         confPath = os.path.join(path, conf)
# #         if os.path.isdir(confPath):
# #             confList.append(conf)
# #     # print(confList)

# #     reactants = []  # forward
# #     reactmp = []
# #     products = []  # reverse
# #     prodtmp = []
# #     TS_abs_ener = []

# #     for conf in sorted(os.listdir(path), key=str):
# #         confPath = os.path.join(path, conf)
# #         if os.path.isdir(confPath):
# #             for f in os.listdir(confPath):
# #                 if f.endswith('irc_f.traj'):
# #                     ftrajPath = os.path.join(confPath, f)
# #                     ftrajFile = read(ftrajPath)[36:]
# #                     if ftrajFile.get_distance(0, 1) < distCutoff:
# #                         reactmp.append(f)
# #                     elif ftrajFile.get_distance(0, 1) > distCutoff:
# #                         prodtmp.append(f)
# #                 elif f.endswith('irc_r.traj'):
# #                     ftrajPath = os.path.join(confPath, f)
# #                     ftrajFile = read(ftrajPath)[36:]
# #                     # for el in reactmp:
# #                     #     if f_prefix in el:
# #                     #         print(f)
# #                     if ftrajFile.get_distance(0, 1) < distCutoff:
# #                         reactmp.append(f)
# #                     elif ftrajFile.get_distance(0, 1) > distCutoff:
# #                         prodtmp.append(f)
# #     # print(reactmp)
# #     # print(prodtmp)

# #     # This is a messy part. It works but far from beeing clean
# #     tmp_r = []
# #     good_index_react = []
# #     good_reactants = []
# #     for num, r in enumerate(reactmp):
# #         r_prefix = r[:3]
# #         if r_prefix not in tmp_r:
# #             tmp_r.append(r_prefix)
# #             good_index_react.append(num)
# #     for inx in good_index_react:
# #         good_reactants.append(reactmp[inx])

# #     set_diff_react = set(reactmp) - set(good_reactants)
# #     diff_react = list(set_diff_react)
# #     prodtmp = diff_react + prodtmp

# #     tmp_p = []
# #     good_index_prod = []
# #     good_products = []
# #     for num, p in enumerate(prodtmp):
# #         p_prefix = p[:3]
# #         if p_prefix not in tmp_p:
# #             tmp_p.append(p_prefix)
# #             good_index_prod.append(num)
# #     for inx in good_index_prod:
# #         good_products.append(prodtmp[inx])
# #     good_products = sorted(good_products)

# #     set_diff_prod = set(prodtmp) - set(good_products)
# #     diff_prod = list(set_diff_prod)
# #     good_reactants = sorted(diff_prod + good_reactants)
# #     # print(good_reactants)
# #     # print(good_products)
# #     #

# #     for geom in good_reactants:
# #         ffile = os.path.join(path, geom[:-5] + '.out')
# #         with open(ffile, 'r') as f:
# #             enerVal = None
# #             reactant_ener = []
# #             for line in f:
# #                 if 'IRC' in line:
# #                     line = line.split()
# #                     # get the energy of irc_f structure (last matching
# #                     # line)
# #                     enerVal = float(line[3])
# #                     reactant_ener.append(line[3])
# #             if enerVal is not None:
# #                 reactants.append(enerVal)
# #         f.close()
# #         print(reactant_ener)
# #         # TS_abs_ener.append(float(reactant_ener[0]))

# #     for geom in good_products:
# #         ffile = os.path.join(path, geom[:-5] + '.out')
# #         with open(ffile, 'r') as f:
# #             enerVal = None
# #             for line in f:
# #                 if 'IRC' in line:
# #                     line = line.split()
# #                     # get the energy of irc_f structure (last matching
# #                     # line)
# #                     enerVal = float(line[3])
# #             if enerVal is not None:
# #                 products.append(enerVal)
# #         f.close()

# #     react_ener = []
# #     for react, prod in zip(reactants, products):
# #         react_ener.append("{:8.3f}".format((prod - react) * 23.06035 * 4.184))

# #     TS_ener = []
# #     for react, ts in zip(reactants, TS_abs_ener):
# #         TS_ener.append("{:8.3f}".format((ts - react) * 23.06035 * 4.184))

# #     print('CO* + O* --> CO2* + *' + '\t' + 'CO* + O* --> TS')
# #     for rener, ts in zip(react_ener, TS_ener):
# #         print(rener + ' kJ/mol' + '\t' + '\t' + ts + ' kJ/mol')




# # def set_up_irc_choosen(path, list_struc, pytemplate_f, pytemplate_r):
# #     ts_path_list = []
# #     for struc in list_struc:
# #         struc_path = os.path.join(path, struc)
# #         ts_path_list.append(struc_path)

# #     for tsDir in ts_path_list:
# #         for traj in os.listdir(tsDir):
# #             if traj.endswith('traj'):
# #                 trajPath = os.path.join(tsDir, traj)
# #                 xyzPath = os.path.join(tsDir, traj[:-5] + '_final.xyz')
# #                 write(xyzPath, read(trajPath))
# #             if traj.endswith('final.xyz'):
# #                 fname = os.path.join(traj[:-10] + '_ts')
# #                 rxn = fname[4:][:4]
# #                 prefix = traj[:3]
# #                 facetpath = os.path.split(path)
# #                 irc_dir = os.path.join(facetpath[0], 'IRC', prefix)
# #                 irc_py_file = os.path.join(facetpath[0], 'IRC')
# #                 os.makedirs(irc_dir, exist_ok=True)
# #                 ts_src_TSxyz_path = os.path.join(tsDir, traj)
# #                 ts_dest_path = os.path.join(irc_dir, fname)
# #                 write(ts_dest_path + '.xyz', read(ts_src_TSxyz_path))
# #                 write(ts_dest_path + '.png', read(ts_src_TSxyz_path))
# #                 TS_xyz = os.path.join(prefix, fname + '.xyz')
# #                 # create run job scripts
# #                 with open(pytemplate_f, 'r') as f:
# #                     template = f.read()
# #                     job_name_f = os.path.join(irc_py_file, fname + '_irc_f.py')
# #                     with open(job_name_f, 'w') as f:
# #                         f.write(template.format(
# #                             prefix=prefix, rxn=rxn, TS_xyz=TS_xyz))
# #                         f.close()
# #                 f.close()
# #                 with open(pytemplate_r, 'r') as r:
# #                     template = r.read()
# #                     job_name_r = os.path.join(irc_py_file, fname + '_irc_r.py')
# #                     with open(job_name_r, 'w') as r:
# #                         r.write(template.format(
# #                             prefix=prefix, rxn=rxn, TS_xyz=TS_xyz))
# #                         r.close()
# #                 r.close()


# # def set_up_irc_ts_estimate(facetpath, TSdir, pytemplate_f, pytemplate_r):
# #     # get ts geoms and create png
# #     ts_path = os.path.join(facetpath, TSdir + '_unique/')
# #     for conf in os.listdir(ts_path):
# #         conf_path = os.path.join(ts_path, conf)
# #         if os.path.isdir(conf_path):
# #             for TStraj in os.listdir(conf_path):
# #                 if TStraj.endswith('.traj'):
# #                     srcTStraj = os.path.join(conf_path, TStraj)
# #                     destTStraj = os.path.join(
# #                         conf_path, TStraj[:-5] + '_final.xyz')
# #                     write(destTStraj, read(srcTStraj))
# #                 if TStraj.endswith('.xyz'):
# #                     fname = os.path.join(TStraj[:8] + '_ts')
# #                     # fname_irc = TSxyz[:6]
# #                     rxn = fname[4:][:4]
# #                     prefix = TStraj[:3]
# #                     irc_dir = os.path.join(facetpath, 'IRC_' + TSdir, prefix)
# #                     irc_py_file = os.path.join(facetpath, 'IRC_' + TSdir)
# #                     os.makedirs(irc_dir, exist_ok=True)
# #                     ts_src_TSxyz_path = os.path.join(conf_path, TStraj)
# #                     ts_dest_path = os.path.join(irc_dir, fname)
# #                     write(ts_dest_path + '.xyz', read(ts_src_TSxyz_path))
# #                     write(ts_dest_path + '.png', read(ts_src_TSxyz_path))
# #                     TS_xyz = os.path.join(prefix, fname + '.xyz')
# #                     # create run job scripts
# #                     with open(pytemplate_f, 'r') as f:
# #                         template = f.read()
# #                         job_name_f = os.path.join(irc_py_file, fname + '_irc_f.py')
# #                         with open(job_name_f, 'w') as f:
# #                             f.write(template.format(
# #                                 prefix=prefix, rxn=rxn, TS_xyz=TS_xyz))
# #                             f.close()
# #                     f.close()
# #                     with open(pytemplate_r, 'r') as r:
# #                         template = r.read()
# #                         job_name_r = os.path.join(irc_py_file, fname + '_irc_r.py')
# #                         with open(job_name_r, 'w') as r:
# #                             r.write(template.format(
# #                                 prefix=prefix, rxn=rxn, TS_xyz=TS_xyz))
# #                             r.close()
# #                     r.close()



# def gen_xyz_from_traj(avDistPath, species):
#     # if species == 'C':
#     #     species_path = os.path.join(avDistPath, species + 'O')
#     # else:
#     species_path = os.path.join(avDistPath, species)
#     for traj in sorted(os.listdir(species_path), key=str):
#         if traj.endswith('.traj'):
#             srctraj_path = os.path.join(species_path, traj)
#             destraj_path = os.path.join(species_path, traj[:-5] + '_final.xyz')
#             write(destraj_path, read(srctraj_path))


# def get_unique_minima_index_after_opt(facetpath, minima_dir):
#     good_minima = []
#     result_list = []
#     unique_minima_index = []
#     gpath = os.path.join(facetpath, minima_dir)
#     trajlist = sorted(Path(gpath).glob('*final.xyz'), key=str)
#     for traj in trajlist:
#         minima = read(traj)
#         minima.pbc = True
#         comparator = SymmetryEquivalenceCheck()
#         result = comparator.compare(minima, good_minima)
#         result_list.append(result)
#         if result is False:
#             good_minima.append(minima)
#     for num, res in enumerate(result_list):
#         if res is False:
#             unique_minima_index.append(str(num).zfill(2))
#     return unique_minima_index


# def get_av_dist(avDistPath, species, scfactor_surface, scaled=False):
#     surface_atoms_indices = []
#     adsorbate_atoms_indices = []
#     all_dists = []
#     gen_xyz_from_traj(avDistPath, species)
#     species_path = os.path.join(avDistPath, species)
#     if species in ['CH3O', 'CH2O']:
#         species = 'O'
#     if len(species) > 1:
#         species = species[:1]
#     # get unique minima indices
#     unique_minima_indices = get_unique_minima_index_after_opt(
#         'Cu_111', os.path.join('minima', species))
#     # go through all indices of *final.xyz file
#     # e.g. 00_final.xyz, 01_final.xyz
#     for index in unique_minima_indices:
#         paths_to_uniq_minima_final_xyz = Path(
#             species_path).glob('{}*final.xyz'.format(index))
#         for unique_minimum_final_xyz in paths_to_uniq_minima_final_xyz:
#             unique_minimum_atom = read(unique_minimum_final_xyz)
#             # open the *final.xyz file as xyz_file
#             with open(unique_minimum_final_xyz, 'r') as f:
#                 xyz_file = f.readlines()
#                 # find all Cu atoms and adsorbate atoms
#                 for num, line in enumerate(xyz_file):
#                     # For Cu(111) we can put ' 1 ' instead of 'Cu ' to limit
#                     # calculations to surface Cu atoms only. For Cu(211) ase is
#                     # not generating tags like this, so full calculation have
#                     # to be performed, i.e. all surface atoms'''
#                     if 'Cu ' in line:
#                         surface_atoms_indices.append(num - 2)
#                     elif species in line and not 'Cu' in line:
#                         # We need to have additional statement
#                         # 'not 'Cu' in line'
#                         # because for C it does not work without it'''
#                         adsorbate_atoms_indices.append(num - 2)
#             f.close()
#             # find the shortest distance between the adsorbate and the surface
#             dist = float(min(unique_minimum_atom.get_distances(
#                 adsorbate_atoms_indices[0], surface_atoms_indices)))
#         all_dists.append(dist)
#     # apply scaling factor if required
#     if scaled:
#         av_dist = mean(all_dists) * scfactor_surface
#     else:
#         av_dist = mean(all_dists)
#     return av_dist


# def get_bond_dist(ads_atom, geom):
#     ''' Specify adsorbate atom symbol and bond distance with the closest surface metal atom will be calculated '''
#     surface_atom = []
#     adsorbate_atom = []
#     struc = read(geom)
#     if len(ads_atom) > 1:
#         ads_atom = ads_atom[:-1]
#     # if ads_atom == 'C':
#     #     ads_atom == 'CO'
#     with open(geom, 'r') as f:
#         xyz_geom_file = f.readlines()
#         for num, line in enumerate(xyz_geom_file):
#             if ' 1 ' in line:
#                 # reading each line and adding only the one with tags = 1
#                 # (surface atoms). It is necessary to subtract 2 as indices in
#                 # atom object starts with 0 and we also have to take into
#                 # account that the first line in xyz file contains non xyz
#                 # information
#                 surface_atom.append(num - 2)
#             elif ads_atom in line:
#                 if not "Cu" in line:
#                     adsorbate_atom.append(num - 2)
#     f.close()

#     if len(adsorbate_atom) > 1:
#         dist_Cu_adsorbate = min(struc.get_distances(
#             adsorbate_atom[0], surface_atom))
#         return dist_Cu_adsorbate
#     else:
#         dist_Cu_adsorbate = min(
#             struc.get_distances(adsorbate_atom, surface_atom))
#         return dist_Cu_adsorbate


# def get_index_adatom(ads_atom, geom):
#     ''' Specify adsorbate atom symbol and its index will be returned '''
#     adsorbate_atom = []
#     if len(ads_atom) > 1:
#         ads_atom = ads_atom[:1]
#     with open(geom, 'r') as f:
#         xyz_geom_file = f.readlines()
#         for num, line in enumerate(xyz_geom_file):
#             if ads_atom in line:
#                 if not 'Cu' in line:
#                     adsorbate_atom.append(num - 2)
#     f.close()
#     return adsorbate_atom[0]


# def get_index_surface_atom(ads_atom, geom):
#     ''' Specify adsorbate atom symbol and index of the neares metal atom will be returned '''
#     surface_atom = []
#     adsorbate_atom = []
#     if len(ads_atom) > 1:
#         ads_atom = ads_atom[:1]
#     # print(ads_atom)
#     struc = read(geom)
#     with open(geom, 'r') as f:
#         xyz_geom_file = f.readlines()
#         for num, line in enumerate(xyz_geom_file):
#             if 'Cu ' in line:
#                 surface_atom.append(num - 2)
#             elif ads_atom in line:
#                 if not 'Cu' in line:
#                     adsorbate_atom.append(num - 2)
#     f.close()

#     all_dist_surface_adsorbate = struc.get_distances(
#         adsorbate_atom[0], surface_atom)
#     min_dist_surface_adsorbate = min(
#         struc.get_distances(adsorbate_atom[0], surface_atom))
#     # get index of the surface atom for which distance to adsorbate is the
#     # lowest
#     index = np.where(all_dist_surface_adsorbate == min_dist_surface_adsorbate)

#     return surface_atom[index[0][0]]


# def check_symm(path):
#     good_adsorbate = []
#     result_list = []
#     geomlist = sorted(Path(path).glob('**/*.traj'))
#     for geom in geomlist:
#         adsorbed = read(geom)
#         adsorbed.pbc = True
#         comparator = SymmetryEquivalenceCheck()
#         result = comparator.compare(adsorbed, good_adsorbate)
#         result_list.append(result)
#         if result is False:
#             good_adsorbate.append(adsorbed)
#     unique_index = []
#     for num, res in enumerate(result_list):
#         if res is False:
#             unique_index.append(str(num).zfill(3))
#     return unique_index


''' The code below is rather useless - double check it '''


# def check_symm_before_xtb(path):
#     good_adsorbate = []
#     result_list = []
#     geomlist = sorted(Path(path).glob('*.xyz'))
#     for geom in geomlist:
#         adsorbed = read(geom)
#         adsorbed.pbc = True
#         comparator = SymmetryEquivalenceCheck()
#         result = comparator.compare(adsorbed, good_adsorbate)
#         result_list.append(result)
#         if result is False:
#             good_adsorbate.append(adsorbed)
#     not_unique_index = []
#     for num, res in enumerate(result_list):
#         if res is True:
#             ''' Better to have all symmetry equivalent site here in a list. The workflow will remove them in getTSestimate function keeping all symmetry distinct sites '''
#             not_unique_index.append(str(num).zfill(3))
#     return not_unique_index


# def create_all_TS(facetpath):
#     all_TS_candidate_path = os.path.join(facetpath, 'TS_candidate')
#     os.makedirs(all_TS_candidate_path, exist_ok=True)
#     gpath = os.path.join(facetpath, 'TS_estimate')
#     for geom in os.listdir(gpath):
#         geomDir = os.path.join(gpath, geom)
#         if os.path.isdir(geomDir):
#             for traj in sorted(os.listdir(geomDir)):
#                 if traj.endswith('.traj'):
#                     # print(str(num/2).zfill(3))
#                     xyz_dir_path = os.path.join(
#                         all_TS_candidate_path, traj[:3])
#                     os.makedirs(xyz_dir_path, exist_ok=True)
#                     srcFile = os.path.join(geomDir, traj)
#                     destFile = os.path.join(xyz_dir_path, traj[:-5])
#                     write(destFile + '.xyz', read(srcFile))
#                     write(destFile + '_initial.png', read(srcFile))


# def create_all_TS_job_files(facetpath, pytemplate):
#     all_TS_candidate_path = os.path.join(facetpath, 'TS_candidate')
#     with open(pytemplate, 'r') as f:
#         pytemplate = f.read()

#     for struc in os.listdir(all_TS_candidate_path):
#         TSdir = os.path.join(all_TS_candidate_path, struc)
#         if os.path.isdir(TSdir):
#             TSpath = os.path.join(all_TS_candidate_path, struc)
#             for file in os.listdir(TSpath):
#                 if file.endswith('.xyz'):
#                     fname = os.path.join(
#                         all_TS_candidate_path, f'{struc}_{file[:-4][4:]}_relax.py')
#                     with open(fname, 'w') as f:
#                         f.write(pytemplate.format(TS=os.path.join(
#                             struc, file), rxn=file[:-4][4:], prefix=file[:3]))
#                     f.close()
#     f.close()


# def create_TS_unique_job_files(facetpath, TSdir, pytemplate):
#     unique_TS_candidate_path = os.path.join(facetpath, TSdir + '_unique')
#     with open(pytemplate, 'r') as f:
#         pytemplate = f.read()

#     for struc in os.listdir(unique_TS_candidate_path):
#         TSdir = os.path.join(unique_TS_candidate_path, struc)
#         if os.path.isdir(TSdir):
#             TSpath = os.path.join(unique_TS_candidate_path, struc)
#             for fl in os.listdir(TSpath):
#                 if fl.endswith('.xyz'):
#                     fname = os.path.join(
#                         unique_TS_candidate_path, fl[:-4] + '.py')
#                     with open(fname, 'w') as f:
#                         f.write(pytemplate.format(TS=os.path.join(
#                             struc, fl), rxn=fl[3:-7], prefix=fl[:2]))
#                     f.close()
#     f.close()


# def set_up_penalty_xtb(path, pytemplate, repeats, slab, species_list, scfactor_surface, scaled1, scaled2):
#     '''Species3 is/should be optional '''
#     # print(species3)
#     fPath = os.path.split(path)
#     path_to_minima = os.path.join(fPath[0], 'minima')
#     avDist1 = get_av_dist(path_to_minima, species_list[0], scfactor_surface, scaled1)
#     avDist2 = get_av_dist(path_to_minima, species_list[1], scfactor_surface, scaled2)
#     '''the code belowe does not work in the loop, probably two differet path_to_minima variable needed to accouut for the poscible scenaario that one species was already calculated (other set of calculations - different reactions, whereas the second in calculated here for the first time) '''
#     # checkMinimaPath = os.path.dirname(os.getcwd())
#     # spList = [species1, species2]

#     # for species in spList:
#     #     is_it_calculated = CheckIfMinimasAlreadyCalculated(
#     #         checkMinimaPath, species)
#     #     if is_it_calculated is False:
#     #         fPath = os.path.split(path)
#     #         path_to_minima = os.path.join(fPath[0], 'minima')
#     #     else:
#     #         path_to_minima = is_it_calculated[1]

#     # species_list = [species1, species2, species3]

#     with open(pytemplate, 'r') as f:
#         pytemplate = f.read()

#     for geom in sorted(os.listdir(path)):
#         if geom.endswith('.xyz'):
#             bonds = []
#             geom_path = os.path.join(path, geom)

#             for species in species_list:
#                 sp_index = get_index_adatom(species, geom_path)
#                 Cu_index = get_index_surface_atom(species, geom_path)
#                 # print(sp_index, Cu_index)
#                 bonds.append((sp_index, Cu_index))

#             # sp1_index = get_index_adatom(species1, geom_path)
#             # sp2_index = get_index_adatom(species2, geom_path)
#             # # sp3_index = get_index_adatom(species3, geom_path)
#             # Cu_index1 = get_index_surface_atom(species1, geom_path)
#             # Cu_index2 = get_index_surface_atom(species2, geom_path)
#             # # Cu_index3 = get_index_surface_atom(species3, geom_path)

#             prefix = geom.split('_')
#             calcDir = os.path.join(path, prefix[0])
#             os.makedirs(calcDir, exist_ok=True)

#             geom_name = geom[:-4]

#             # List of bonds we want O-Cu; H-Cu
#             # bonds = ((sp1_index, Cu_index1),
#             #          (sp2_index, Cu_index2))
#             # print(type(bonds))
#             avDists = (avDist1, avDist2)
#             traj_path = os.path.join(geom[:-4] + '.traj')
#             init_png = os.path.join(calcDir, geom[:-4] + '_initial.png')
#             write(init_png, read(geom_path))
#             fname = os.path.join(calcDir, geom[:-4] + '.py')
#             with open(fname, 'w') as f:
#                 f.write(pytemplate.format(geom=geom, bonds=bonds,
#                                           avDists=avDists, traj_path=traj_path,
#                                           repeats=repeats, prefix=prefix[0],
#                                           geom_name=geom_name, slabopt=slab))
#             f.close()
#             shutil.move(geom_path, calcDir)
#     f.close()

#     try:
#         rmpath = os.path.join(path, 'initial_png')
#         shutil.rmtree(rmpath)
#     except FileNotFoundError:
#         pass
#         # print('No files to delete')


# # def create_unique_TS(facetpath, TSdir):
# #     # checkSymmPath = os.path.join(path, 'TS_estimate')

# #     # gd_ads_index = check_symm(facetpath, TSdir)
# #     # new checksym with globbing
# #     gd_ads_index = check_symm(os.path.join(facetpath, TSdir))
# #     for i, index in enumerate(gd_ads_index):
# #         uniqueTSdir = os.path.join(
# #             facetpath, TSdir + '_unique', str(i).zfill(2))
# #         if os.path.isdir(uniqueTSdir):
# #             shutil.rmtree(uniqueTSdir)
# #             os.makedirs(uniqueTSdir, exist_ok=True)
# #         else:
# #             os.makedirs(uniqueTSdir, exist_ok=True)
# #         gpath = os.path.join(facetpath, TSdir)
# #         for geom in os.listdir(gpath):
# #             geomDir = os.path.join(gpath, geom)
# #             if os.path.isdir(geomDir):
# #                 for traj in sorted(os.listdir(geomDir)):
# #                     if traj.startswith(gd_ads_index[i]) and traj.endswith('.traj'):
# #                         srcFile = os.path.join(geomDir, traj)
# #                         destFile = os.path.join(uniqueTSdir, traj[:-5] + '_ts')
# #                         write(destFile + '.xyz', read(srcFile))
# #                         write(destFile + '.png', read(srcFile))
# #         #         shutil.copy2(os.path.join(path, geom), uniqueTSdir)
# #         for ts in os.listdir(uniqueTSdir):
# #             # if neb.endswith('.xyz'):
# #             TS_xyz_Dir = os.path.join(uniqueTSdir, ts)
# #             newTS_xyz_Dir = os.path.join(
# #                 uniqueTSdir, str(i).zfill(2) + ts[3:])
# #             # NEB_png_Dir = os.path.join(uniqueTSdir, str(
# #             #     i).zfill(3) + neb[2:][:-4] + '.png')
# #             os.rename(TS_xyz_Dir, newTS_xyz_Dir)
# #             # write(NEB_png_Dir, read(newNEB_xyz_Dir))
