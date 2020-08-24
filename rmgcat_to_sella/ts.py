from rmgcat_to_sella.gratoms import Gratoms
from rmgcat_to_sella.adsorption import Builder
from rmgcat_to_sella.molecules import molecule, get_3D_positions

from rmgcat_to_sella.adsorbates import Adsorbates
from rmgcat_to_sella.graph_utils import node_test
from rmgcat_to_sella.main import WorkFlow

from ase.io import read, write
from ase.utils.structure_comparator import SymmetryEquivalenceCheck

import numpy as np
import yaml
import os
import shutil
from statistics import mean
from pathlib import Path
import networkx as nx
from spglib import get_symmetry


class TS():
    def __init__(
        self, facetpath, slab, ts_dir, yamlfile, repeats, creation_dir
    ):
        ''' Initializing

        Parameters:
        ___________
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        slab : str
            a '.xyz' file name with the optimized slab
            e.g.
            'Cu_111_slab_opt.xyz'
        ts_dir : str
            a path to directory with TSs
            e.g. 'TS_estimate'
        yamlfile : str
            a name of the .yaml file with reaction list
        repeats: tuple
            specify reapeats in (x, y, z) direction,
            eg. (3, 3, 1)

        '''
        self.facetpath = facetpath
        self.slab = slab
        self.ts_dir = ts_dir
        self.yamlfile = yamlfile
        self.repeats = repeats
        self.creation_dir = creation_dir

    def prepare_ts_estimate(self, scfactor, scfactor_surface,
                            rotAngle, pytemplate_xtb, species,
                            scaled1, scaled2):
        ''' Prepare TS estimates for subsequent xTB calculations

        Parameters:
        ___________
        scfator : float
            a scaling factor to scale a bond distance between
            atoms taking part in the reaction
            e.g. 1.4
        scfactor_surface : float
            a scaling factor to scale the target bond distance, i.e.
            the average distance between adsorbed atom and the nearest
            surface atom. Helpful e.g. when H is far away form the surface
            in TS, whereas for minima it is close to the surface
            e.g. 1.0
        rotAngle : float
            an angle (deg) of rotation  of the TS guess adduct on the surface
            e.g. 60
        pytemplate_xtb : python script
            a template file for penalty function minimization job
        species : list(str)
            a list of max 2 species that take part in the reaction
            e.g. ['O', 'H'] or ['CO2', 'H']
        scaled1 : bool
            specify whether use the optional scfactor_surface
            for the species 1 (sp1)
        scaled2 : bool
            specify whether use the optional scfactor_surface
            for the species 2 (sp2)
        '''

        ts_estimate_path = os.path.join(self.facetpath, self.ts_dir)
        if os.path.exists(ts_estimate_path):
            shutil.rmtree(ts_estimate_path)
            os.makedirs(ts_estimate_path)
        else:
            os.makedirs(ts_estimate_path)

        r_name_list, p_name_list, images = TS.prepare_react_list(self)
        rxn_name = TS.get_rxn_name(self)

        TS.TS_placer(self, scfactor, rotAngle,
                     rxn_name, r_name_list, p_name_list, images)

        TS.filtered_out_equiv_ts_estimate(self, rxn_name)
        TS.set_up_penalty_xtb(self, pytemplate_xtb,
                              species, scaled1, scaled2, scfactor_surface)

    def get_max_rot_angle(self):
        ''' Get the maximum angle of rotation for a given slab that will
        generate all symmetrically distinct TS estimates.

        Returns:
        ________
        max_rot_angle : float
            a maximum angle of rotation of TS adducts on the given slab that
            would generate symmetrically dostinct configurations
            i.e. once the TS adduct is rotated by bigger angle than
            max_rot_angle, some of the generated structures will be
            symetrically equivalent to the others

        '''
        # convert slab to ASE's Atom object
        slab = read(self.slab)
        # get all symmetry operations
        symm = get_symmetry(slab)
        # count rotations operations
        nrot = len(symm['rotations'])
        # the angle of rotation is 360 devided by nrot/2 as thera are equal
        # number of symmetry operations around the z axis
        # (e.g. for Cu_111 therea are 6 z and 6 -z operations)
        max_rot_angle = 360 / (nrot / 2)
        # max_rot_angle = 360/(nrot/4) # lets try 120 angle
        return max_rot_angle

    def prepare_react_list(self):
        '''Convert yaml file to more useful format

        Returns
        _______
        r_name_list : list(str)
            a list with all reactants for the given reaction
        p_name_list : list(str)
            a list with all products for the given reaction
        images : list(Gratoms)
            a list of CatKit's Gratom object (both reactants and products)

        '''
        with open(self.yamlfile, 'r') as f:
            yamltxt = f.read()
        reactions = yaml.safe_load(yamltxt)
        speciesInd = []
        bonds = []
        unique_species = []
        unique_bonds = []
        images = []

        put_adsorbates = Adsorbates(
            self.facetpath, self.slab, self.repeats, self.yamlfile,
            self.creation_dir
        )

        for rxn in reactions:
            # transforming reactions data to gratom objects
            reactants, rbonds = put_adsorbates.rmgcat_to_gratoms(
                rxn['reactant'].split('\n'))
            products, pbonds = put_adsorbates.rmgcat_to_gratoms(
                rxn['product'].split('\n'))
            speciesInd += reactants + products
            bonds += rbonds + pbonds

        # check if any products are the same as any reactants
        for species1, bond in zip(speciesInd, bonds):
            for species2 in unique_species:
                if nx.is_isomorphic(species1.graph, species2.graph, node_test):
                    break
            else:
                images.append(get_3D_positions(species1))
                unique_species.append(species1)
                unique_bonds.append(bond)

        r_name_list = [str(species.symbols) for species in reactants]
        p_name_list = [str(species.symbols) for species in products]

        return r_name_list, p_name_list, images

    def get_rxn_name(self):
        ''' Get the reaction name

        Returns
        _______
        The name of the reaction in the following format:
        OH_H+O
        '''
        r_name_list, p_name_list, _ = TS.prepare_react_list(self)

        r_name = '+'.join([species for species in r_name_list])
        p_name = '+'.join([species for species in p_name_list])

        rxn_name = r_name + '_' + p_name
        return rxn_name

    def TS_placer(self, scfactor, rotAngle, rxn_name,
                  r_name_list, p_name_list, images):
        ''' Place adsorbates on the surface to estimate TS

        Parameters
        __________
        scfator : float
            a scaling factor to scale a bond distance between
            atoms taking part in the reaction
            e.g. 1.4
        rotAngle : float
            an angle (deg) of rotation  of the TS guess adduct on the surface
            e.g. 60
        rxn_name : str
            a reaction name
            e.g OH_O+H
        r_name_list : list(str)
            a list with all reactants for the given reaction
        p_name_list : list(str)
            a list with all products for the given reaction
        images : list(Gratoms)
            a list of CatKit's Gratom object (both reactants and products)

        '''
        ts_estimate_path = os.path.join(self.facetpath, self.ts_dir)
        slab = read(self.slab)
        slab.pbc = (True, True, False)
        # ADSORBATES
        # This part here is a bit hard coded, especially when dealing with
        # reactants with >= 3 atoms. The code works for couple of species
        # included in if else statement below. Probably will fail for other
        # species. This is becouse ase.build.molecule function numbers atoms
        # very diferently depending on the molecule specify.
        # To be resolved somehow.

        atom1 = 0  # atom1 should always have index 0, i.e. the first atom
        atom2 = 1
        # well, here is the problem becouse ase.build.molecule counts atoms
        # differently than my input in the yaml file. Probably the problem with
        # the way yamlfile is converted to species - to be resolved
        if len(r_name_list) <= len(p_name_list):
            TS_candidate = molecule(r_name_list[0])[0]
        else:
            # images[2] keeps info about product. It's a Gratom object.
            # Cannot use catkit's molecule method becouse it generates
            # gemetry based on wierd order chemical formula string. Sometimes
            # it works to use molecule method but in general it is better to
            # use get_3D_positions() as it converts yaml info into 3d
            # position. In images I keep the oryginal order of elements with
            # connectivity info. In TS_candidate.get_chemical_formula()
            # elements are sorted and grouped
            TS_candidate = images[2]
            # reactName = r_name

        ''' Different cases '''
        if TS_candidate.get_chemical_formula() == 'CHO':
            atom2 = 2
            bondedThrough = [0]
        elif TS_candidate.get_chemical_formula() == 'COH':
            atom2 = 2
            bondedThrough = [0]
        elif TS_candidate.get_chemical_formula() == 'CHO2':
            bondedThrough = [2]  # connect through oxygen

        elif TS_candidate.get_chemical_formula() == 'C2H5O2':
            atom0 = 0
            atom1 = 1
            atom2 = 5
            bondedThrough = [6]  # connect through oxygen
        else:
            bondedThrough = [0]

        bondlen = TS_candidate.get_distance(atom1, atom2)

        # Final ridgid rotations to orientate the TS_candidate on the surface
        if len(TS_candidate.get_tags()) < 3:
            TS_candidate.rotate(90, 'y')
            TS_candidate.set_distance(
                atom1, atom2, bondlen * scfactor, fix=0)
        elif len(TS_candidate.get_tags()) == 3:
            TS_candidate.rotate(90, 'z')
            TS_candidate.set_distance(
                atom1, atom2, bondlen * scfactor, fix=0)
        # else:
        elif TS_candidate.get_chemical_formula() == 'C2H5O2':
            # TS_candidate.rotate(60, 'y')
            TS_candidate.rotate(90, 'z')
            ''' Should work for H2CO*OCH3, i.e. COH3+HCOH '''
            TS_candidate.set_angle(atom2, atom1, atom0, -45,
                                   indices=[0, 1, 2, 3, 4], add=True)
            TS_candidate.set_distance(
                atom1, atom2, bondlen * scfactor, fix=1,
                indices=[0, 1, 2, 3, 4]
            )
            # indices=[0, 1, 2, 3, 4]
        elif TS_candidate.get_chemical_formula() == 'CHO2':
            TS_candidate.rotate(90, 'z')
            TS_candidate.set_distance(
                atom1, atom2, bondlen * scfactor, fix=0)
        # double check this
        put_adsorbates = Adsorbates(
            self.facetpath, self.slab, self.repeats, self.yamlfile,
            self.creation_dir
        )
        slabedges, tags = put_adsorbates.get_edges(self)
        # double check this
        grslab = Gratoms(numbers=slab.numbers,
                         positions=slab.positions,
                         cell=slab.cell,
                         pbc=slab.pbc,
                         edges=slabedges)
        grslab.arrays['surface_atoms'] = tags

        # building adsorbtion structures
        ads_builder = Builder(grslab)

        # max_angle = int(TS.get_max_rot_angle(self))
        # do a full 360 degree rotation - bridge have different symmetry than
        # hollows and top sites on Cu(111), so cannot use 60 degree for all
        # of them. So the approach is to do a full 360 degree scan with
        # 5 degree increment, check the symmetry using check_symm_before_xtb()
        # method, run the penalty function minimization for the remaining
        # structures. Once done, check symmetry again using check_symm() method
        # That would avoid bugs if other surface is applied.
        max_angle = 360
        angle = 0
        count = 0
        step_size = 5
        while angle <= max_angle:
            structs = ads_builder.add_adsorbate(
                TS_candidate, bondedThrough, -1, auto_construct=False)
            # change to True will make bondedThrough work.
            # Now it uses TS_candidate,rotate...
            # to generate adsorbed strucutres
            big_slab = slab * self.repeats
            nslab = len(slab)

            for i, struc in enumerate(structs):
                big_slab_ads = big_slab + struc[nslab:]
                write(os.path.join(ts_estimate_path, '{}'.format(
                    str(i + len(structs) * count).zfill(3)
                ) + '_' + rxn_name + '.xyz'), big_slab_ads)

            TS_candidate.rotate(step_size, 'z')
            angle += step_size
            count += 1

    def filtered_out_equiv_ts_estimate(self, rxn_name):
        '''Filtering out symmetry equivalent sites

        Parameters
        __________
        rxn_name : str
            a reaction name
            e.g. OH_O+H
        '''
        ts_estimate_path = os.path.join(self.facetpath, self.ts_dir)
        filtered_equivalen_sites = TS.check_symm_before_xtb(
            self, ts_estimate_path)
        for eqsites in filtered_equivalen_sites:
            try:
                fileToRemove = os.path.join(
                    ts_estimate_path, eqsites + '_' + rxn_name + '.xyz')
                os.remove(fileToRemove)
            except OSError:
                print("Error while deleting file : ", fileToRemove)

        for prefix, noneqsites in enumerate(
            sorted(os.listdir(ts_estimate_path))
        ):
            prefix = str(prefix).zfill(3)
            oldfname = os.path.join(ts_estimate_path, noneqsites)
            newfname = os.path.join(
                ts_estimate_path, prefix + noneqsites[3:])
            os.rename(oldfname, newfname)

    def get_average_distance_all(self, species, geom_path):
        # probably do not need it - check
        raise NotImplementedError

    def set_up_penalty_xtb(self, pytemplate, species_list,
                           scaled1, scaled2, scfactor_surface):
        ''' Prepare calculations of the penalty function

        Parameters
        __________

        pytemplate : python script
            a template for the penalty function calculations
        species_list : list(str)
            a list of max 2 species that take part in the reaction
            e.g. ['O', 'H'] or ['CO2', 'H']
            TODO: remove this limitation
        scaled1 : bool
            specify whether use the optional scfactor_surface
            for the species 1 (sp1)
        scaled2 : bool
            specify whether use the optional scfactor_surface
            for the species 2 (sp2)
        scfactor_surface : float
            a scaling factor to scale the target bond distance, i.e.
            the average distance between adsorbed atom and the nearest
            surface atom. Helpful e.g. when H is far away form the surface
            in TS, whereas for minima it is close to the surface
            e.g. 1.0

        '''
        ts_estimate_path = os.path.join(self.facetpath, self.ts_dir)
        path_to_minima = os.path.join(self.facetpath, 'minima')
        rxn_name = TS.get_rxn_name(self)

        with open(pytemplate, 'r') as f:
            pytemplate = f.read()

        # get a list with the average distances (only symmetrically distinct
        # sites considered) for all species. It can be done outside the nested
        # loop, as it is enough to calculate it only once.
        average_distance_list = []
        for species in species_list:
            av_dist = TS.get_av_dist(self, path_to_minima, species,
                                     scfactor_surface, scaled1)
            average_distance_list.append(av_dist)

        # get all ts_estimatex_xyz files in alphabetic order
        ts_estimates_xyz_files = sorted(os.listdir(ts_estimate_path))

        # loop through all .xyz files
        for prefix, xyz_file in enumerate(ts_estimates_xyz_files):
            if xyz_file.endswith('.xyz'):
                bonds = []
                xyz_file_path = os.path.join(ts_estimate_path, xyz_file)
                # loop through all species
                for species in species_list:
                    sp_index = TS.get_index_adatom(
                        self, species, xyz_file_path)
                    Cu_index = TS.get_index_surface_atom(
                        self, species, xyz_file_path)
                    bonds.append((sp_index, Cu_index))

                # set up variables
                av_dists_tuple = tuple(average_distance_list)
                prefix = str(prefix).zfill(3)
                calcDir = os.path.join(ts_estimate_path, prefix)
                os.makedirs(calcDir, exist_ok=True)
                geom_name = prefix + '_' + rxn_name
                traj_path = os.path.join(xyz_file[:-4] + '.traj')
                fname = os.path.join(calcDir, xyz_file[:-4] + '.py')

                # create job_file
                with open(fname, 'w') as f:
                    f.write(pytemplate.format(geom=xyz_file, bonds=bonds,
                                              av_dists_tuple=av_dists_tuple,
                                              traj_path=traj_path,
                                              repeats=self.repeats,
                                              prefix=prefix,
                                              geom_name=geom_name,
                                              slabopt=self.slab))
                f.close()
                # write .png files
                init_png = os.path.join(
                    calcDir, xyz_file[:-4] + '_initial.png')
                write(init_png, read(xyz_file_path))
                # remove .xyz files
                shutil.move(xyz_file_path, calcDir)
        f.close()

        # remove initial_png directory, if it exists
        try:
            rmpath = os.path.join(ts_estimate_path, 'initial_png')
            shutil.rmtree(rmpath)
        except FileNotFoundError:
            print('No files to delete')
            pass

    def get_av_dist(self, path_to_minima, species, scfactor_surface,
                    scaled=False):
        ''' Get the average bond distance between adsorbate atom and
        the nearest surface atom for all symmetrically distinct minima

        Parameters:
        ___________
        path_to_minima : str
            a path to minima
            e.g. 'Cu_111/minima'
        species : str
            a species symbol
            e.g. 'H' or 'CO'
        scfactor_surface : float
            a scaling factor to scale the target bond distance, i.e.
            the average distance between adsorbed atom and the nearest
            surface atom. Helpful e.g. when H is far away form the surface
            in TS, whereas for minima it is close to the surface
            e.g. 1.0
        scaled : bool
            specify whether to use the optional scfactor_surface
            for the given species
            default = False

        Returns:
        ________
        av_dist : float
            an average bond distance between the given species and the nearest
            surface atom for all symmetrically dictinct minima

        '''
        surface_atoms_indices = []
        adsorbate_atoms_indices = []
        all_dists = []
        TS.get_xyz_from_traj(self, path_to_minima, species)
        species_path = os.path.join(path_to_minima, species)
        # treat the special cases
        if species in ['CH3O', 'CH2O']:
            species = 'O'
        # deal with the multiatomic molecules and look only for the surface
        # bonded atom
        if len(species) > 1:
            sp_bonded = species[:1]
        else:
            # for one atomic species it is trivial
            sp_bonded = species
        # get unique minima indices
        unique_minima_indices = TS.get_unique_minima_indicies_after_opt(
            self, path_to_minima, species
        )
        # go through all indices of *final.xyz file
        # e.g. 00_final.xyz, 01_final.xyz
        for index in unique_minima_indices:
            paths_to_uniq_minima_final_xyz = Path(
                species_path).glob('{}*final.xyz'.format(index))
            for unique_minimum_final_xyz in paths_to_uniq_minima_final_xyz:
                unique_minimum_atom = read(unique_minimum_final_xyz)
                # open the *final.xyz file as xyz_file
                with open(unique_minimum_final_xyz, 'r') as f:
                    xyz_file = f.readlines()
                    # find all Cu atoms and adsorbate atoms
                    for num, line in enumerate(xyz_file):
                        # For Cu(111) we can put ' 1 ' instead of 'Cu ' to
                        # limit calculations to surface Cu atoms only.
                        # For Cu(211) ASE is not generating tags like this,
                        # so full calculation have to be performed, i.e. all
                        # surface atoms
                        if 'Cu ' in line:
                            surface_atoms_indices.append(num - 2)
                        elif sp_bonded in line and 'Cu ' not in line:
                            # We need to have additional statement
                            # 'not 'Cu' in line'
                            # because for C it does not work without it'''
                            adsorbate_atoms_indices.append(num - 2)
                f.close()
                # find the shortest distance between the adsorbate
                # and the surface
                dist = float(min(unique_minimum_atom.get_distances(
                    adsorbate_atoms_indices[0], surface_atoms_indices)))
            all_dists.append(dist)
        # apply scaling factor if required
        if scaled:
            av_dist = mean(all_dists) * scfactor_surface
        else:
            av_dist = mean(all_dists)
        return av_dist

    def get_xyz_from_traj(self, path_to_minima, species):
        ''' Convert all ASE's traj files to .xyz files for a given species

        Parameters:
        ___________
        path_to_minima : str
            a path to minima
            e.g. 'Cu_111/minima'
        species : str
            a species symbol
            e.g. 'H' or 'CO'

        '''
        species_path = os.path.join(path_to_minima, species)
        for traj in sorted(os.listdir(species_path), key=str):
            if traj.endswith('.traj'):
                src_traj_path = os.path.join(species_path, traj)
                des_traj_path = os.path.join(
                    species_path, traj[:-5] + '_final.xyz')
                write(des_traj_path, read(src_traj_path))

    def get_unique_minima_indicies_after_opt(self, path_to_minima, species):
        ''' Get the indicies of the symmetrically distinct minima
        for a given species

        Parameters:
        ___________
        path_to_minima : str
            a path to minima
            e.g. 'Cu_111/minima'
        species : str
            a species symbol
            e.g. 'H' or 'CO'

        Returns:
        ________

        unique_minima_indices : list(str)
            a list with indecies of all unique minima for a given species
            e.g. ['01', '02', '04']

        '''
        good_minima = []
        result_list = []
        unique_minima_indices = []
        path_to_species = os.path.join(path_to_minima, species)
        trajlist = sorted(Path(path_to_species).glob('*final.xyz'), key=str)
        for traj in trajlist:
            minima = read(traj)
            minima.pbc = True
            comparator = SymmetryEquivalenceCheck()
            result = comparator.compare(minima, good_minima)
            result_list.append(result)
            if result is False:
                good_minima.append(minima)
        for prefix, result in enumerate(result_list):
            if result is False:
                unique_minima_indices.append(str(prefix).zfill(2))
        return unique_minima_indices

    def copy_minimas_prev_calculated(self, current_dir, species_list,
                                     minima_dir):
        ''' If minimas have been already calculated in different set of
         reactions, they are copied to the current workflow and used instead
         of calculating it again

        Parameters
        __________

        current_dir : str
            a path to the current directory
            e.g './Cu_111'
        species_list : list(str)
            a list of max 2 species that take part in the reaction
            e.g. ['O', 'H'] or ['CO2', 'H']
        minima_dir : str
            a path to minima directory
            e.g. Cu_111/minima

        '''
        for species in species_list:
            is_it_calculated = WorkFlow().check_if_minima_already_calculated(
                current_dir, species, self.facetpath)
            if is_it_calculated[0] is False:
                pass
            else:
                try:
                    _, copy_path_dft, copy_path_outfiles = is_it_calculated
                    dst_dir = os.path.join(minima_dir, species)
                    # copy DFT files
                    shutil.copytree(copy_path_dft, dst_dir)
                    # copy Sella's .out files
                    for outfile in copy_path_outfiles:
                        shutil.copy2(outfile, minima_dir)
                except FileExistsError:
                    raise FileExistsError('Files already copied')

    def create_unique_TS(self):
        ''' Create unique TS files for calculations '''
        gd_ads_index = TS.check_symm(
            self, os.path.join(self.facetpath, self.ts_dir))
        for i, index in enumerate(gd_ads_index):
            uniqueTSdir = os.path.join(
                self.facetpath, self.ts_dir + '_unique', str(i).zfill(2))
            if os.path.isdir(uniqueTSdir):
                shutil.rmtree(uniqueTSdir)
                os.makedirs(uniqueTSdir, exist_ok=True)
            else:
                os.makedirs(uniqueTSdir, exist_ok=True)
            gpath = os.path.join(self.facetpath, self.ts_dir)
            for geom in os.listdir(gpath):
                geomDir = os.path.join(gpath, geom)
                if os.path.isdir(geomDir):
                    for traj in sorted(os.listdir(geomDir)):
                        if (
                            traj.startswith(gd_ads_index[i])
                            and traj.endswith('.traj')
                        ):
                            srcFile = os.path.join(geomDir, traj)
                            destFile = os.path.join(
                                uniqueTSdir, traj[:-5] + '_ts')
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

    def create_TS_unique_job_files(
        self, pytemplate,
        pseudopotentials, pseudo_dir,
        balsam_exe_settings,
        calc_keywords
    ):
        ''' Create job submission files'''
        unique_TS_candidate_path = os.path.join(
            self.facetpath, self.ts_dir + '_unique')
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
                                struc, fl), rxn=fl[3:-7], prefix=fl[:2],
                                pseudopotentials=pseudopotentials,
                                pseudo_dir=pseudo_dir,
                                balsam_exe_settings=balsam_exe_settings,
                                calc_keywords=calc_keywords,
                                creation_dir=self.creation_dir
                            ))
                        f.close()
        f.close()

    def get_bond_dist(self, ads_atom, geom):
        ''' Specify adsorbate atom symbol and bond distance with the closest
            surface metal atom will be calculated.

        Parameters:
        ___________
        ads_atom : str
            an atom of the species bonded to the surface, e.g. 'O' for OH
        geom : str
            a .xyz of .traj file name with geometry of the structure
            to be analysed

        Returns:
        ________
        dist_Cu_adsorbate : float
            distance between ads_atom and the closest surface atom

        '''
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
                    # TODO: possibly a bug when other surface is used,
                    # e.g. Cu_211
                    #
                    # reading each line and adding only the one with tags = 1
                    # (surface atoms). It is necessary to subtract 2 as indices
                    # in atom object starts with 0 and we also have to take
                    # into account that the first line in xyz file contains
                    # non xyz information
                    surface_atom.append(num - 2)
                elif ads_atom in line:
                    if "Cu" not in line:
                        adsorbate_atom.append(num - 2)
        f.close()

        if len(adsorbate_atom) > 1:
            dist_Cu_adsorbate = min(struc.get_distances(
                adsorbate_atom[0], surface_atom))
            return dist_Cu_adsorbate
        else:
            dist_Cu_adsorbate = min(
                struc.get_distances(adsorbate_atom, surface_atom))
            return dist_Cu_adsorbate

    def get_index_adatom(self, ads_atom, geom):
        ''' Specify adsorbate atom symbol and its index will be returned.

        Parameters:
        ___________
        ads_atom : str
            an atom of the species bonded to the surface, e.g. 'O' for OH
        geom : str
            a .xyz of .traj file name with geometry of the structure
            to be analysed

        Returns:
        ________
        adsorbate_atom[0] : float
            index of the adsorbed atom

        '''
        # TODO: currently it works only for one ads_atom, so it will return
        # only index at [0]. In the future I plan to add support for many
        # indices
        adsorbate_atom = []
        if len(ads_atom) > 1:
            ads_atom = ads_atom[:1]
        with open(geom, 'r') as f:
            xyz_geom_file = f.readlines()
            for num, line in enumerate(xyz_geom_file):
                if ads_atom in line:
                    if 'Cu' not in line:
                        adsorbate_atom.append(num - 2)
        f.close()
        return adsorbate_atom[0]

    def get_index_surface_atom(self, ads_atom, geom):
        ''' Specify adsorbate atom symbol and index of the nearest metal atom
            will be returned.

        Parameters:
        ___________
        ads_atom : str
            an atom of the species bonded to the surface, e.g. 'O' for OH
        geom : str
            a .xyz of .traj file name with geometry of the structure
            to be analysed

        Returns:
        ________
        surface_atom[index[0][0]] : float
            index of the metal atom to which ads_atom is bonded

        '''
        surface_atom = []
        adsorbate_atom = []
        if len(ads_atom) > 1:
            ads_atom = ads_atom[:1]
        struc = read(geom)
        with open(geom, 'r') as f:
            xyz_geom_file = f.readlines()
            for num, line in enumerate(xyz_geom_file):
                if 'Cu ' in line:
                    surface_atom.append(num - 2)
                elif ads_atom in line:
                    if 'Cu' not in line:
                        adsorbate_atom.append(num - 2)
        f.close()

        all_dist_surface_adsorbate = struc.get_distances(
            adsorbate_atom[0], surface_atom)
        min_dist_surface_adsorbate = min(
            struc.get_distances(adsorbate_atom[0], surface_atom))
        # get index of the surface atom for which distance to adsorbate is the
        # lowest
        index = np.where(all_dist_surface_adsorbate
                         == min_dist_surface_adsorbate)

        return surface_atom[index[0][0]]

    def check_symm(self, path):
        ''' Check for the symmetry equivalent structures in the given path

        Parameters:
        ___________
        path : str
            a path to a directory where are files to be checked,
            e.g. Cu_111/TS_estimate_unique

        Returns:
        ________
        unique_index : list(str)
            a list with prefixes of all symmetrically distinct sites

        '''
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

    def check_symm_before_xtb(self, path):
        ''' Check for the symmetry equivalent structures in the given path
            before executing penalty function minimization

        Parameters:
        ___________
        path : str
            a path to a directory where are files to be checked,
            e.g. Cu_111/TS_estimate

        Returns:
        ________
        not_unique_index : list(str)
            a list with prefixes of all symmetry equivalent structures

        '''
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
        not_unique_index = []
        for num, res in enumerate(result_list):
            if res is True:
                # Better to have all symmetry equivalent site here in a list.
                # The workflow will remove them in getTSestimate function
                # keeping all symmetry distinct sites
                not_unique_index.append(str(num).zfill(3))
        return not_unique_index
