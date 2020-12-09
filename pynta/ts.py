from typing import List, Tuple, Dict
from pynta.excatkit.gratoms import Gratoms
from pynta.excatkit.adsorption import Builder
from pynta.general_ts_guesses import GeneralTSGuessesGenerator
from pynta.adsorbates import Adsorbates
from pynta.io import IO

from ase.io import read, write
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from ase.collections import g2
from ase import Atoms

from collections import Counter

import numpy as np
import os
import shutil
from statistics import mean
from pathlib import Path, PosixPath
from spglib import get_symmetry


class TS():
    def __init__(
            self,
            facetpath: str,
            slab: str,
            ts_estimate_dir: str,
            yamlfile: str,
            repeats: Tuple[int, int, int],
            creation_dir: PosixPath) -> None:
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
        ts_estimate_dir : str
            a path to directory with TSs
            e.g. 'TS_estimate'
        yamlfile : str
            a name of the .yaml file with a reactions list
        repeats : tuple(int, int, int)
            how to replicate unit cell in (x, y, z) direction,
            e.g. (3, 3, 1)
        creation_dir : posix
            a posix path to the main working directory

        '''
        self.facetpath = facetpath
        self.slab = slab
        self.ts_estimate_dir = ts_estimate_dir
        self.yamlfile = yamlfile
        self.repeats = repeats
        self.creation_dir = creation_dir
        self.n_kpts = IO().get_kpoints(self.repeats)
        self.path_to_slab = os.path.join(self.creation_dir, self.slab)
        big_slab = read(self.path_to_slab) * self.repeats
        self.nslab = len(big_slab)

    def prepare_ts_estimate(
            self,
            rxn: Dict[str, str],
            scfactor: float,
            scfactor_surface: float,
            pytemplate_xtb: str,
            species_list: List[str],
            reacting_atoms: List[str],
            metal_atom: str,
            scaled1: bool,
            scaled2: bool) -> None:
        ''' Prepare TS estimates for subsequent xTB calculations

        Parameters:
        ___________

        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
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
        pytemplate_xtb : python script
            a template file for penalty function minimization job
        species_list : list(str)
            a list of species which atoms take part in the reaction,
            i.e. for ['CO2'] ['C'] is taking part in reaction
            e.g. ['O', 'H'] or ['CO2', 'H']
        easier_to_build : list(str)
            a list of species that are considerd as the reactiong one
        metal_atom : str
            a checmical symbol for the surface atoms (only metallic surfaces
            are allowed)
        scaled1 : bool
            specify whether use the optional scfactor_surface
            for the species 1 (sp1)
        scaled2 : bool
            specify whether use the optional scfactor_surface
            for the species 2 (sp2)

        '''
        r_name_list, p_name_list = IO.get_reactants_and_products(rxn)
        rxn_name = IO.get_rxn_name(rxn)

        ts_estimate_path = os.path.join(
            self.creation_dir,
            self.facetpath,
            rxn_name,
            self.ts_estimate_dir)

        self.TS_placer(
            ts_estimate_path,
            scfactor,
            rxn,
            rxn_name,
            r_name_list,
            p_name_list,
            reacting_atoms)

        self.filtered_out_equiv_ts_estimate(
            ts_estimate_path,
            rxn_name)

        self.set_up_penalty_xtb(
            ts_estimate_path,
            pytemplate_xtb,
            species_list,
            reacting_atoms,
            metal_atom,
            scaled1,
            scfactor_surface)

    def get_max_rot_angle(self) -> None:
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

    def TS_placer(
            self,
            ts_estimate_path: str,
            scfactor: float,
            rxn: Dict[str, str],
            rxn_name: str,
            r_name_list: List[str],
            p_name_list: List[str],
            reacting_atoms: Dict[str, int]) -> None:
        ''' Place adsorbates on the surface to estimate TS

        Parameters:
        ___________

        ts_estimate_path : str
            a path to TS_estimate directory,
            e.g. {creation_dir}/Cu_111/TS_estimate
        scfator : float
            a scaling factor to scale a bond distance between
            atoms taking part in the reaction
            e.g. 1.4
        rxn_name : str
            a reaction name
            e.g OH_O+H
        r_name_list : list(str)
            a list with all reactants for the given reaction
        p_name_list : list(str)
            a list with all products for the given reaction

        '''
        # create TS_estimate directory
        os.makedirs(ts_estimate_path, exist_ok=True)

        slab_atom = read(self.path_to_slab)
        slab_atom.pbc = (True, True, False)

        # decide whether reactant or product is easier to use
        # to build a ts_guess
        if len(r_name_list) <= len(p_name_list):
            easier_to_build = 'reactant'
            ts_estimators = r_name_list
        else:
            easier_to_build = 'product'
            ts_estimators = p_name_list

        # Developing assuming there is only one element in the list, e.g 'OH'
        # So the reaction is A + B -> C
        # Potentially, one can imagine a reaction where A + B -> C + D with
        # all species adsorbed to the surface - there will be 2 elements in
        # ts_estimators
        for ts_est in ts_estimators:
            ts_guess_generator = GeneralTSGuessesGenerator(
                ts_est, rxn, rxn_name, easier_to_build, scfactor)
            ts_guess, bonded_idx = ts_guess_generator.decide(
                reacting_atoms.values())

            # convert slab (Atom) to grslab(Gratom)
            ads_builder = self.prepare_slab_for_ts_guess(slab_atom)

            # put ts_guess on the surface in all possible spots and rotate
            self.ts_guess_place_and_rotate(ts_guess, ads_builder,
                                           bonded_idx, slab_atom,
                                           ts_estimate_path, rxn_name)

    def prepare_slab_for_ts_guess(
            self,
            slab_atom: Atoms) -> Builder:
        ''' Convert slab_atom to Gratom object and finally to Builder object
            to have a suitable format for ts_guess placement

        Parameters
        ----------
        slab_atom : Atoms
            an initial slab converted to Atoms object, e.g.
            ase.io.read('Cu_111_slab_opt.xyz')

        Returns
        -------
        ads_builder : Builder
            a slab converted to a Builder object, suitable for ts_guess
            placement

        '''
        put_adsorbates = Adsorbates(self.facetpath, self.slab, self.repeats,
                                    self.yamlfile, self.creation_dir)
        slabedges, tags = put_adsorbates.get_edges(self)
        grslab = Gratoms(numbers=slab_atom.numbers,
                         positions=slab_atom.positions,
                         cell=slab_atom.cell,
                         pbc=slab_atom.pbc,
                         edges=slabedges)
        grslab.arrays['surface_atoms'] = tags

        ads_builder = Builder(grslab)
        return ads_builder

    def ts_guess_place_and_rotate(
            self,
            ts_guess: Gratoms,
            ads_builder: Builder,
            bonded_idx: int,
            slab_atom: Atoms,
            ts_estimate_path: str,
            rxn_name: str,
            max_angle: int = 360,
            angle: int = 0,
            angle_increment: int = 5) -> None:
        ''' Place ts_est on the surface (ads_builder) in all possible
            locations, rotate and save all obtained structures as a
            seperate .xyz files

        Parameters
        ----------
        ts_guess : Gratoms
            an estimate TS structure that will be placed on the surface
        ads_builder : Builder
            a slab converted to a Builder object, suitable for ts_guess
            placement
        bonded_idx : int
            an index of ts_guess atom(s) connected to the surface
        slab_atom : Atoms
            an initial slab converted to Atoms object, e.g.
            ase.io.read('Cu_111_slab_opt.xyz')
        ts_estimate_path : str
            a path where save all generated .xyz files
        rxn_name : str
            a reaction name
            e.g. OH_O+H
        max_angle : int, optional
            maximum angle of rotation around z axis, by default 360
        angle : int, optional
            initial angle of rotation around z axis, by default 0
        angle_increment : int, optional
            an increment of angle of rotation around z axis, by default 5

        '''
        count = 0
        while angle <= max_angle:
            structs = ads_builder.add_adsorbate(
                ts_guess, [bonded_idx], -1, auto_construct=False)
            # change to True will make bonded_through work.
            # Now it uses ts_guess,rotate...
            # to generate adsorbed strucutres
            big_slab = slab_atom * self.repeats
            nslab = len(slab_atom)

            for i, struc in enumerate(structs):
                big_slab_ads = big_slab + struc[nslab:]
                prefix = str(i + len(structs) * count).zfill(3)
                xyz_fname = os.path.join(
                    prefix + '_' + self.facetpath + '_' + rxn_name + '.xyz')
                path_to_xyz = os.path.join(ts_estimate_path, xyz_fname)
                write(path_to_xyz, big_slab_ads)

            ts_guess.rotate(angle_increment, 'z')
            angle += angle_increment
            count += 1

    def filtered_out_equiv_ts_estimate(
            self,
            ts_estimate_path: str,
            rxn_name: str) -> None:
        '''Filtered out symmetry equivalent sites and remove them keeping
            only symmetry disctinct structures.

        Parameters
        __________

        ts_estimate_path : str
            a path to TS_estimate directory,
            e.g. {creation_dir}/Cu_111/TS_estimate
        rxn_name : str
            a reaction name
            e.g. OH_O+H

        '''
        # check the symmetry
        filtered_equivalent_sites = TS.check_symm_before_xtb(
            ts_estimate_path)

        # remove all symmetry equivalent structures
        for eqsites in filtered_equivalent_sites:
            file_to_remove = os.path.join(
                ts_estimate_path, eqsites + '_' + self.facetpath + '_' +
                rxn_name + '.xyz')
            try:
                os.remove(file_to_remove)
            except OSError:
                print("Error while deleting file : ", file_to_remove)

        # rename and organize symmetry distinct structures
        for prefix, noneqsites in enumerate(
            sorted(os.listdir(ts_estimate_path))
        ):
            prefix = str(prefix).zfill(3)
            old_fname = os.path.join(ts_estimate_path, noneqsites)
            new_fname = os.path.join(
                ts_estimate_path, prefix + noneqsites[3:])
            os.rename(old_fname, new_fname)

    def set_up_penalty_xtb(
            self,
            ts_estimate_path: str,
            pytemplate: str,
            species_list: List[str],
            reacting_atoms: List[str],
            metal_atom: str,
            scaled1: bool,
            scfactor_surface,) -> None:
        ''' Prepare calculations of the penalty function

        Parameters:
        __________

        ts_estimate_path : str
            a path to TS_estimate directory,
            e.g. {creation_dir}/Cu_111/TS_estimate
        rxn_name : str
            a reaction name
            e.g. OH_O+H
        pytemplate : python script
            a template for the penalty function calculations
        species_list : list(str)
            a list of species which atoms take part in the reaction,
            i.e. for ['CO2'] ['C'] is taking part in reaction
            e.g. ['O', 'H'] or ['CO2', 'H']
        easier_to_build : list(str)
            a list of species that are considerd as the reactiong one
        metal_atom : str
            a checmical symbol for the surface atoms (only metallic surfaces
            are allowed)
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
        # load pytemplate
        with open(pytemplate, 'r') as f:
            pytemplate = f.read()

        # set up path to optimized minima files
        path_to_minima = os.path.join(
            self.creation_dir, self.facetpath, 'minima')

        # get a dictionary with average distances for all species in
        # species_list, e.g. {'CO2': 4.14.., 'H': 1.665..., 'O': 1.847...}
        sp_surf_av_dists = TS.get_av_dists_dict(
            species_list, metal_atom, path_to_minima, scfactor_surface,
            scaled1)

        # convert sp_surf_av_dists dict to a tuple and take into accout that
        # the same type of species can be included into penalty function
        # calculations many times, e.g. ['C', 'H', 'O', 'O'], so the av_dist
        # for the 'O' have to be specified twice (order not important)
        av_dists_tuple = TS.get_av_dists_tuple(
            reacting_atoms.values(), sp_surf_av_dists)

        # get all .xyz files with TS estimates
        ts_estimates_xyz_files = []
        ts_est = Path(ts_estimate_path).glob('*.xyz')
        for ts in ts_est:
            ts_estimates_xyz_files.append(str(ts))

        # sort it in increasing order
        ts_estimates_xyz_files = sorted(ts_estimates_xyz_files)

        # take a first file and use it as a template to get info about
        # surface atom and adsorbate atoms indices
        tmp_ts_atom = read(ts_estimates_xyz_files[0])

        # create surface_atoms_idx dict with all surface atoms and theirs
        # corresponding indicies
        surface_atoms_idxs = {
            atom.symbol + '_' + str(atom.index): atom.index
            for atom in tmp_ts_atom if atom.symbol == metal_atom}

        # Main loop #
        # Loop through all .xyz files
        for prefix, xyz_file in enumerate(ts_estimates_xyz_files):
            bonds = []

            # Loop through all relevant_species
            for idx in reacting_atoms.values():
                # sp_index = TS.get_sp_index(
                #     species, visited_species, adsorbate_atoms_idxs)
                adsorbed_atom_idx = idx + self.nslab

                metal_idx = TS.get_index_surface_atom(
                    adsorbed_atom_idx, surface_atoms_idxs, xyz_file)

                bonds.append((adsorbed_atom_idx, metal_idx))

            # set up some variables
            prefix = str(prefix).zfill(3)
            calc_dir = os.path.join(ts_estimate_path, prefix)
            os.makedirs(calc_dir, exist_ok=True)
            f_name_xyz = os.path.basename(xyz_file)[:-4]
            traj_path = os.path.join(f_name_xyz + '.traj')
            fname = os.path.join(calc_dir, f_name_xyz + '.py')

            # create job_file
            with open(fname, 'w') as f:
                f.write(pytemplate.format(geom=os.path.basename(xyz_file),
                                          bonds=bonds,
                                          av_dists_tuple=av_dists_tuple,
                                          creation_dir=self.creation_dir,
                                          traj_path=traj_path,
                                          repeats=self.repeats,
                                          prefix=prefix,
                                          geom_name=f_name_xyz,
                                          slabopt=self.slab))
            # move .xyz file
            shutil.move(xyz_file, calc_dir)

    @staticmethod
    def is_valid_sp_index(
            species: str,
            sp_index: int,
            adsorbate_atoms_idxs: Dict[str, int]) -> int:
        ''' Return True if given sp_index is possible for a given species

        Parameters
        ----------
        species : str
            a species symbol
            e.g. 'H', 'C' or 'O', etc.
        sp_index : int
            an index for a given species
        adsorbate_atoms_idxs : dict(str:int)
            a dictionary with all adsorbate atoms and theirs corresponding
            indicies

        Returns
        -------
        bool
            True if given sp_index is possible for a given species.
            False otherwise

        '''
        key = species + '_' + str(sp_index)
        try:
            adsorbate_atoms_idxs[key]
            return True
        except KeyError:
            return False

    @staticmethod
    def get_av_dists_tuple(
            easier_to_build: List[str],
            sp_surf_av_dists: Dict[str, float]):
        ''' Create a av_dists_tuple with all relevant average bond distances.

            This method loops through sp_surf_av_dists dictionary. If
            particular key exists n > 1 times in easier_to_build, this
            entry is added n times to a new dictionary. Otherwise, a new
            dictionary is updated with the keys and values of sp_surf_av_dists.
            At the end, the valuses of the new dict are transformed into tuple

        Parameters
        ----------
        easier_to_build : list(str)
            a list of species that are considerd as the reactiong one
        sp_surf_av_dists : dict(str:float)
            a dictionary with keys being species name and average distances
            as values

        Returns
        -------
        av_dist_tuple : tuple
            a tuple with all relevant average bond distances

        '''
        count = Counter(easier_to_build)
        av_dists_dict = {}
        for key, value in sp_surf_av_dists.items():
            n = count[key]
            if n > 1:
                for i in range(1, n):
                    av_dists_dict[key + '_' + str(i)] = value
            av_dists_dict[key] = value

        av_dist_tuple = tuple(av_dists_dict.values())
        return av_dist_tuple

    @staticmethod
    def get_av_dists_dict(
            species_list: str,
            metal_atom: str,
            path_to_minima: str,
            scfactor_surface: float,
            scaled1: bool) -> Dict[str, float]:
        ''' Get a dictionary with average distances for all species in
            species_list

        Parameters
        ----------
        species_list : list(str)
            a list of species which atoms take part in the reaction,
            i.e. for ['CO2'] ['C'] is taking part in reaction
            e.g. ['O', 'H'] or ['CO2', 'H']
        minima_dir : str
            a path to minima directory
            e.g. Cu_111/minima
        path_to_minima : str
            a path to minima
            e.g. 'Cu_111/minima'
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

        Returns
        -------
        sp_surf_av_dists : dict(str:float)
            a dictionary with keys being species name and average distances
            as values

        '''
        sp_surf_av_dists = {}
        for species in species_list:
            av_dist = TS.get_av_dist(
                path_to_minima, species, metal_atom, scfactor_surface, scaled1)
            sp_surf_av_dists[species] = av_dist
        return sp_surf_av_dists

    @staticmethod
    def get_av_dist(
            path_to_minima: str,
            species: str,
            metal_atom: str,
            scfactor_surface: float,
            scaled: bool = False) -> float:
        ''' Get the average bond distance between a given adsorbate atom and
        the nearest surface atom for all symmetrically distinct minima

        Parameters:
        ___________
        path_to_minima : str
            a path to minima
            e.g. 'Cu_111/minima'
        species : str
            a species symbol
            e.g. 'H' or 'CO'
        metal_atom : str
            a checmical symbol for the surface atoms (only metallic surfaces
            are allowed)
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
        path_to_species = os.path.join(path_to_minima, species)

        # get unique minima prefixes
        unique_minima_prefixes = IO.get_unique_prefixes(path_to_species)

        # choose a representative traj file based on which surface_atom_idxs
        # and adsorbate_atom_idxs will be created
        path_to_tmp_traj = os.path.join(path_to_species, '00.traj')
        tmp_traj = read(path_to_tmp_traj)

        # create a list with all surface atom indicies
        surface_atom_idxs = [
            atom.index for atom in tmp_traj if atom.symbol == metal_atom]

        # create a list with all adosorbate atom indicies
        adsorbate_atom_idxs = {
            atom.symbol + '_' + str(atom.index): atom.index
            for atom in tmp_traj if atom.symbol != metal_atom}

        # loop through all unique traj files, e.g. 00.traj, 01.traj ...
        all_dists_bonded = []
        for index in unique_minima_prefixes:
            path_to_unique_minima_traj = os.path.join(
                path_to_species, '{}.traj'.format(index))
            uq_species_atom = read(path_to_unique_minima_traj)
            ads_atom_surf_dist = {}
            for key, ads_atom_idx in adsorbate_atom_idxs.items():
                # create a dict storing the shortest distance between given
                # adsorbate index atom and surface atoms
                ads_atom_surf_dist[key] = min(uq_species_atom.get_distances(
                    ads_atom_idx, surface_atom_idxs))
            # the shortest distance in bonded_ads_atom_surf_dist.values()
            # is considered as a distance between atom bonded to the surface
            # and the surface
            bonded_ads_atom_surf_dist = min(ads_atom_surf_dist.values())
            all_dists_bonded.append(bonded_ads_atom_surf_dist)
        # apply the scalling factor
        if scaled:
            av_dist = mean(all_dists_bonded) * scfactor_surface
        else:
            av_dist = mean(all_dists_bonded)
        return av_dist

    def create_unique_ts_all(
            self,
            ts_estimate_path: str,
            rxn_name: str,
            pytemplate: str,
            pseudopotentials: Dict[str, str],
            pseudo_dir: str,
            balsam_exe_settings: Dict[str, str],
            calc_keywords: Dict[str, str]) -> None:
        ''' Create all TS_estimate_unique files

        Parameters:
        ___________

        ts_estimate_path : str
            a path to TS_estimate directory,
            e.g. {creation_dir}/Cu_111/TS_estimate
        rxn_name : str
            a reaction name
            e.g. OH_O+H
        pytemplate : python script
            a template file for saddle point minimization with Sella
        pseudopotentials : dict(str: str)
            a dictionary with QE pseudopotentials for all species.
            e.g.
            dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                H='H.pbe-kjpaw_psl.1.0.0.UPF',
                O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                C='C.pbe-n-kjpaw_psl.1.0.0.UPF',
                )
        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            '/home/mgierad/espresso/pseudo'
        balsam_exe_settings : dict(str:str)
            a dictionary with balsam settings
        calc_keywords : dict(str:str)
            a dictionary with keywords to Quantum Espresso calculations

        '''
        # create .xyz and .png files
        TS.create_unique_ts_xyz_and_png(ts_estimate_path)
        # create job files (.py scripts)
        self.create_ts_unique_py_file(pytemplate,
                                      rxn_name,
                                      pseudopotentials,
                                      pseudo_dir,
                                      ts_estimate_path,
                                      balsam_exe_settings,
                                      calc_keywords
                                      )

    @staticmethod
    def create_unique_ts_xyz_and_png(
            ts_estimate_path: str) -> None:
        ''' Create unique TS files for saddle point calculations
            for a given scfactor

        Parameters:
        ___________

        ts_estimate_path : str
            a path to TS_estimate directory,
            e.g. {creation_dir}/Cu_111/TS_estimate

        '''
        # check symmetry of all TS estimates in ts_estimate_path
        gd_ads_index = TS.check_symm(ts_estimate_path)

        for i, index in enumerate(gd_ads_index):
            # name of the directory with unique TS_estimates for which saddle
            # point calculations are to be perfomed
            ts_estimate_unique_path = os.path.join(ts_estimate_path +
                                                   '_unique', str(i).zfill(2))
            # create TS_estimate_unique directory
            if os.path.isdir(ts_estimate_unique_path):
                shutil.rmtree(ts_estimate_unique_path)
                os.makedirs(ts_estimate_unique_path, exist_ok=True)
            else:
                os.makedirs(ts_estimate_unique_path, exist_ok=True)

            # search for trajectories of symmetry distinct structures
            uq_traj_search = '**/{}*traj'.format(gd_ads_index[i])
            trajs = Path(ts_estimate_path).glob(uq_traj_search)
            # loop through all unique trajectory and create .xyz and .png
            for traj in trajs:
                traj = str(traj)
                # split the path, get file name, remove last 5 characters and
                # add sufix '_ts'
                fname = os.path.split(traj)[1][:-5] + '_ts'
                uq_ts_file = os.path.join(ts_estimate_unique_path, fname)
                write(uq_ts_file + '.xyz', read(traj))
                write(uq_ts_file + '.png', read(traj))
            # rename TS to have prefixes in order with no gaps
            # e.g. was 027_OH_O+H_ts.xyz; is 03_OH_O+H_ts.png
            for ts in os.listdir(ts_estimate_unique_path):
                old_ts_name = os.path.join(ts_estimate_unique_path, ts)
                new_ts_name = os.path.join(
                    ts_estimate_unique_path, str(i).zfill(2) + ts[3:])
                os.rename(old_ts_name, new_ts_name)

    def create_ts_unique_py_file(
            self,
            pytemplate: str,
            rxn_name: str,
            pseudopotentials: Dict[str, str],
            pseudo_dir: str,
            ts_estimate_path: str,
            balsam_exe_settings: Dict[str, str],
            calc_keywords: Dict[str, str]) -> None:
        ''' Create job submission files for TS calculation with Sella

        Parameters:
        ___________

        pytemplate : python script
            a template file for saddle point minimization with Sella
        rxn_name : str
            a reaction name
            e.g. OH_O+H
        pseudopotentials : dict(str: str)
            a dictionary with QE pseudopotentials for all species.
            e.g.
            dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                H='H.pbe-kjpaw_psl.1.0.0.UPF',
                O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                C='C.pbe-n-kjpaw_psl.1.0.0.UPF',
                )
        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            '/home/mgierad/espresso/pseudo'
        ts_estimate_path : str
            a path to TS_estimate directory,
            e.g. {creation_dir}/Cu_111/TS_estimate
        balsam_exe_settings : dict(str:str)
            a dictionary with balsam settings
        calc_keywords : dict(str:str)
            a dictionary with keywords to Quantum Espresso calculations

        '''
        ts_estimate_unique_path = ts_estimate_path + '_unique'

        with open(pytemplate, 'r') as f:
            pytemplate = f.read()

        tss = Path(ts_estimate_unique_path).glob('**/*ts.xyz')
        for ts in tss:
            ts = str(ts)
            ts_dir, ts_fname = os.path.split(ts)
            py_dir, _ = os.path.split(ts_dir)
            py_fname = ts_fname[:-4] + '.py'
            py_file = os.path.join(py_dir, py_fname)

            with open(py_file, 'w') as f:
                f.write(pytemplate.format(
                    ts_fname=ts_fname,
                    rxn_name=rxn_name,
                    prefix=ts_fname[:2],
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=self.creation_dir,
                    n_kpts=self.n_kpts,
                    repeats=self.repeats
                ))

    @staticmethod
    def get_index_adatom(
            ads_atom: str,
            adsorbate_atoms_idxs: Dict[str, float]) -> int:
        ''' Specify adsorbate atom symbol and its index will be returned.

        Parameters:
        ___________
        ads_atom : str
            an atom of the species bonded to the surface, e.g. 'O' for OH
        adsorbate_atoms_idxs : dict(str:int)
            a dictionary with all adsorbate atoms and theirs corresponding
            indicies

        Returns:
        ________
        ads_index : int
            index of the adsorbed atom

        '''
        vals = []
        for key, val in adsorbate_atoms_idxs.items():
            if key.startswith(ads_atom):
                vals.append(val)
        if len(vals) > 1:
            ads_index = vals[0] + len(vals) - 1
        else:
            ads_index = vals[0]
        return ads_index

    @staticmethod
    def get_index_surface_atom(
            sp_index: int,
            surface_atoms_idxs: Dict[str, int],
            geom: str) -> int:
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
        surface_atoms[index[0][0]] : float
            index of the metal atom to which ads_atom is bonded

        '''
        surface_atoms_idxs_list = list(surface_atoms_idxs.values())

        ts_est_atom = read(geom)

        all_dist_surface_adsorbate = ts_est_atom.get_distances(
            sp_index, surface_atoms_idxs_list)

        min_dist_surface_adsorbate = min(all_dist_surface_adsorbate)
        # get index of the surface atom for which distance to adsorbate is the
        # lowest
        index = np.where(all_dist_surface_adsorbate
                         == min_dist_surface_adsorbate)
        return surface_atoms_idxs_list[index[0][0]]

    @staticmethod
    def check_symm(
            path: str) -> List[int]:
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

    @staticmethod
    def check_symm_before_xtb(
            path: str) -> List[int]:
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
            a list with prefixes of all symmetry equivalent structures, i.e.,
            files with these prefixes should be deleted

        '''
        comparator = SymmetryEquivalenceCheck()
        geomlist = sorted(Path(path).glob('*.xyz'))

        good_adsorbates = []
        result_list = []

        for geom in geomlist:
            adsorbed = read(geom)
            adsorbed.pbc = True
            result = comparator.compare(adsorbed, good_adsorbates)
            result_list.append(result)
            if result is False:
                # if compared structures are different, add the current one to
                # good_adsorbates
                good_adsorbates.append(adsorbed)
        not_unique_index = []
        for num, res in enumerate(result_list):
            if res is True:
                # Better to have all symmetry equivalent site here in a list.
                # The workflow will remove them in getTSestimate function
                # keeping all symmetry distinct sites
                not_unique_index.append(str(num).zfill(3))
        return not_unique_index
