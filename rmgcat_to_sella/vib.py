from sys import path
from rmgcat_to_sella.io import IO
from ase.io import read, write
from pathlib import Path, PosixPath
from typing import Tuple, Dict, List
from numpy import floor
import os
import shutil


class AfterTS():
    def __init__(
            self,
            facetpath: str,
            yamlfile: str,
            slab: str,
            repeats: Tuple[int, int, int],
            creation_dir: PosixPath) -> None:

        self.facetpath = facetpath
        self.yamlfile = yamlfile
        self.repeats = repeats
        self.creation_dir = creation_dir
        self.slab_atom_path = os.path.join(self.creation_dir, slab)
        self.nslab = len(read(self.slab_atom_path) * self.repeats)
        self.n_kpts = IO().get_kpoints(self.repeats)

    def set_up_ts_vib(
            self,
            rxn: Dict[str, str],
            pytemplate: str,
            balsam_exe_settings: Dict[str, int],
            calc_keywords: Dict[str, str],
            pseudopotentials: Dict[str, str],
            pseudo_dir: str) -> None:
        '''Set up files for TSs vibration calculations

        Parameters
        ----------
        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        pytemplate : python script
            a template for TS_vib calculations
        balsam_exe_settings : dict{str:int}
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.
            balsam_exe_settings = {'num_nodes': 1,
                                   'ranks_per_node': 48,
                                   'threads_per_rank': 1}
        calc_keywords : dict{str:str}
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}
        pseudopotentials : dict{str:str}
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

        '''
        rxn_name = IO().get_rxn_name(rxn)

        ts_estimate_unique_dir = os.path.join(
            self.creation_dir,
            self.facetpath,
            rxn_name,
            'TS_estimate_unique')
        ts_vib_dir = os.path.join(
            self.creation_dir,
            self.facetpath,
            rxn_name,
            'TS_estimate_unique_vib')

        ts_final_geoms = Path(ts_estimate_unique_dir).glob('**/*final.xyz')

        for ts_final_geom in ts_final_geoms:
            ts_final_geom = str(ts_final_geom)
            prefix = ts_final_geom.split('/')[-2]
            ts_vib_dir_prefix = os.path.join(ts_vib_dir, prefix)
            os.makedirs(ts_vib_dir_prefix, exist_ok=True)

            # copy *ts_final.xyz files to ts_vib_dir_prefix - for debug
            shutil.copy2(ts_final_geom, ts_vib_dir_prefix)
            _, geom = os.path.split(ts_final_geom)

            py_fname = ts_vib_dir_prefix + '_' + \
                self.facetpath + '_' + rxn_name + '_ts_vib.py'

            self.create_ts_vib_py_files(
                pytemplate,
                geom,
                py_fname,
                balsam_exe_settings,
                calc_keywords,
                pseudopotentials,
                pseudo_dir)

    def create_ts_vib_py_files(
            self,
            pytemplate: str,
            geom: str,
            py_fname: str,
            balsam_exe_settings: Dict[str, int],
            calc_keywords: Dict[str, str],
            pseudopotentials: Dict[str, str],
            pseudo_dir: str,
            nimages: int = 30,
            n: int = 0) -> None:
        ''' Create job submission .py files for frequency calculation of TSs

        Parameters
        ----------
        pytemplate : python script
            a template for TS_vib calculations
        py_fname : str
            path to the the .py (including the .py file name) file that is
            about to be created
        balsam_exe_settings : dict(str:int)
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.
            balsam_exe_settings = {'num_nodes': 1,
                                'ranks_per_node': 48,
                                'threads_per_rank': 1}
        calc_keywords : dict(str:str)
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}
        pseudopotentials : dict(str:str)
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
        nimages : int, optional
            how many strucutres to use to construct a trajectory visualizing
            oscilations, by default 30
        n : int, optional
            mode of oscilation, i.e. which vibration to analyze.
            0 is the first vibration, should be imaginary, by default 0

        '''
        with open(pytemplate, 'r') as f:
            pytemplate = f.read()
        with open(py_fname, 'w') as f:
            f.write(pytemplate.format(
                geom=geom,
                balsam_exe_settings=balsam_exe_settings,
                calc_keywords=calc_keywords,
                creation_dir=self.creation_dir,
                pseudopotentials=pseudopotentials,
                pseudo_dir=pseudo_dir,
                nimages=nimages,
                n=n,
                repeats=self.repeats,
                n_kpts=self.n_kpts
            ))

    def prepare_opt_after_ts(
            self,
            rxn: Dict[str, str],
            pytemplate: str,
            balsam_exe_settings: Dict[str, int],
            calc_keywords: Dict[str, str],
            pseudopotentials: Dict[str, str],
            pseudo_dir: str) -> None:
        ''' Create files for after_TS calculations - to verify TS structures
            and get corresponding reactant and product minima[summary]

        Parameters
        ----------
        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        pytemplate : python script
            a template for after_TS calculations
        balsam_exe_settings : dict{str:int}
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.
            balsam_exe_settings = {'num_nodes': 1,
                                   'ranks_per_node': 48,
                                   'threads_per_rank': 1}
        calc_keywords : dict{str:str}
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}
        pseudopotentials : dict{str:str}
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

        '''
        rxn_name = IO().get_rxn_name(rxn)
        ts_vib_dir = os.path.join(self.creation_dir,
                                  self.facetpath,
                                  rxn_name,
                                  'TS_estimate_unique_vib')

        vib_traj_files = Path(ts_vib_dir).glob('**/*traj')
        for vib_traj in vib_traj_files:
            vib_traj = str(vib_traj)
            prefix = vib_traj.split('/')[-2]

            after_ts_dir = os.path.join(self.creation_dir,
                                        self.facetpath,
                                        rxn_name,
                                        'after_TS',
                                        prefix)
            os.makedirs(after_ts_dir, exist_ok=True)

            fname = os.path.join(
                prefix + '_' + self.facetpath + '_' + rxn_name + '_after_ts')

            fname_forward = os.path.join(after_ts_dir, fname + '_f')
            fname_reverse = os.path.join(after_ts_dir, fname + '_r')

            AfterTS.get_forward_and_reverse(
                vib_traj,
                fname_forward,
                fname_reverse)

            self.create_after_ts_py_files(
                pytemplate,
                fname_forward,
                fname_reverse,
                balsam_exe_settings,
                calc_keywords,
                pseudopotentials,
                pseudo_dir)

    @staticmethod
    def get_forward_and_reverse(
            vib_traj: str,
            fname_forward: str,
            fname_reverse: str,
            nimages: int = 30) -> None:
        ''' Get forward and reverse .xyz file by nudging TS towards imaginary
        mode of oscilations summary

        Parameters
        ----------
        vib_traj : str
            a path to trajectory file with optimized TS
        fname_forward : str
            a path for forward calculations
        fname_reverse : str
            a path for reverse calculations
            [description]
        nimages : int, optional
            how many strucutres to use to construct a trajectory visualizing
            oscilations, by default 30

        Raises
        ------
        ValueError
            raised if there are more than one *traj file visualizing vibrations

        '''

        # vib_.traj file contains nimages structures visualizing imaginary
        # frequency mode. 0 and nimages/2 are the same structures - no
        # displacement. Max and min displacement is defined as
        # floor(nimages/4) and nimages - floor(nimages/4), respectively.
        # for inmages = 16 it would be like this
        # start...max...start...min...
        index_forward = int(floor(nimages/4))
        index_reverse = int(nimages - index_forward)

        fname = [fname_forward, fname_reverse]
        indices = [index_forward, index_reverse]
        extensions = ['.xyz', '.png']

        # get forward and reverse displacement
        for fn, ind in zip(fname, indices):
            for ext in extensions:
                write(fn + ext, read(vib_traj, index=ind))

    def create_after_ts_py_files(
            self,
            pytemplate: str,
            fname_forward: str,
            fname_reverse: str,
            balsam_exe_settings: Dict[str, int],
            calc_keywords: Dict[str, str],
            pseudopotentials: Dict[str, str],
            pseudo_dir: str) -> None:
        ''' Create job submission files for minimization displaced structures
            after TS

        Parameters
        ----------
        pytemplate : python script
            a template for after_TS calculations
        fname_forward : str
            a path for forward calculations
        fname_reverse : str
            a path for reverse calculations
        balsam_exe_settings : dict(str:int)
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.
            balsam_exe_settings = {'num_nodes': 1,
                                'ranks_per_node': 48,
                                'threads_per_rank': 1}
        calc_keywords : dict(str:str)
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}
        pseudopotentials : dict(str:str)
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

        '''
        with open(pytemplate, 'r') as f:
            pytemplate = f.read()

        for fn in [fname_forward, fname_reverse]:
            tmp, geom = os.path.split(fn)
            fname = os.path.join(os.path.split(tmp)[0], geom + '.py')
            with open(fname, 'w') as f:
                f.write(pytemplate.format(
                    geom=geom,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=self.creation_dir,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    repeats=self.repeats,
                    n_kpts=self.n_kpts
                ))

    def get_all_distances(self) -> None:
        ''' Get distances between reacting species for ts, forward and
        reverse structure

        '''
        all_rxn_names = IO().get_list_all_rxns_names(self.yamlfile)
        for rxn_name in all_rxn_names:
            ts_dist_dict = self.get_ts_dist(rxn_name)
            f_dist_dict, r_dist_dict = AfterTS.get_forward_and_reverse_dist(
                rxn_name)
            AfterTS.print_table(ts_dist_dict, f_dist_dict, r_dist_dict)

    @staticmethod
    def print_table(
            ts_dist_dict: Dict[str, float],
            forward_dist_dict: Dict[str, float],
            reverse_dist_dict: Dict[str, float]) -> None:
        ''' Print information about bond distances (reacting atoms) for TS,
        forward and reverse .xyz file

        Parameters
        ----------
        ts_dist_dict : dict(str:float)
            a dictionary with keys being TS file names while values are
            bond distance beteween of reacting species
        f_dist_dict : dict(str:float)
            a dictionary with keys being paths to forward files while values
            are bond distance beteween of reacting species
        r_dist_dict : dict(str:float)
            a dictionary with keys being paths to reverse files while values
            are bond distance beteween of reacting species

        '''
        keys = ts_dist_dict.keys()
        ts_val = ts_dist_dict.values()
        f_val = forward_dist_dict.values()
        r_val = reverse_dist_dict.values()

        print('TS \t \t TS_dist \t F_dist \t R_dist')
        for key, ts, f, r in zip(keys, ts_val, f_val, r_val):
            print('{} \t {:2f} \t {:2f} \t {:2f}'.format(key, ts, f, r))

    def get_ts_dist(
            self,
            rxn_name: str) -> Dict[str, float]:
        ''' For given rxn_name get distances between reacting species
            in TS structure

        Parameters
        ----------
        rxn_name : str
            a name of the reaction in the following format:
            'OH_H+O'

        Returns
        -------
        ts_dist_dict : dict(str:float)
            a dictionary with keys being TS file names while values are
            bond distance beteween of reacting species

        '''
        ts_dist_dict = {}

        ts_estimate_unique_dir = os.path.join(
            self.facetpath, rxn_name, 'TS_estimate_unique')

        traj_files = Path(ts_estimate_unique_dir).glob('**/*traj')

        for traj in sorted(traj_files):
            traj = str(traj)
            traj_atoms = read(traj)[self.nslab:]
            dist = traj_atoms.get_distance(0, 1)
            traj_fname = traj.split('/')[-1]
            ts_dist_dict[traj_fname] = dist
        return ts_dist_dict

    def get_forward_and_reverse_dist(
            self,
            rxn_name: str) -> Tuple[Dict[str, float], Dict[str, float]]:
        ''' For given rxn_name get distances between reacting species in
            forward and reverse .xyz file

        Parameters
        ----------
        rxn_name : str
            a name of the reaction in the following format:
            'OH_H+O'

        Returns
        -------
        f_dist_dict : dict(str:float)
            a dictionary with keys being paths to forward files while values
            are bond distance beteween of reacting species
        r_dist_dict : dict(str:float)
            a dictionary with keys being paths to reverse files while values
            are bond distance beteween of reacting species

        '''
        f_dist_dict = {}
        r_dist_dict = {}

        after_ts_dir = os.path.join(
            self.facetpath, rxn_name, 'after_TS')

        xyz_files = Path(after_ts_dir).glob('**/*xyz')

        for xyz in sorted(xyz_files):
            xyz = str(xyz)
            if 'forward' in xyz:
                xyz_atom = read(xyz)[self.nslab:]
                dist = xyz_atom.get_distance(0, 1)
                f_dist_dict[xyz] = dist
            elif 'reverse' in xyz:
                xyz_atom = read(xyz)[self.nslab:]
                dist = xyz_atom.get_distance(0, 1)
                r_dist_dict[xyz] = dist
        return f_dist_dict, r_dist_dict


class minimaVib():
    def __init__(
            self,
            facetpath: str,
            creation_dir: PosixPath) -> None:

        self.facetpath = facetpath
        self.creation_dir = creation_dir
        self.minima_path = os.path.join(
            self.creation_dir, self.facetpath, 'minima')

    def create_minima_vib_all(
            self,
            species: str,
            pytemplate: str,
            balsam_exe_settings: Dict[str, int],
            pseudo_dir: str,
            pseudopotentials: Dict[str, str],
            calc_keywords: Dict[str, str],
            creation_dir: PosixPath) -> None:
        ''' Create all files for frequency calculations for the most stable
            conformer for a given species.

        Parameters
        ----------
        species : str
            a chemical symbol of the adsorbate
        species_list : List[str]
            a list with all species taking part in all reactions,
            e.g. ['H', 'C', 'CH', 'O', 'OH']
        pytemplate : str
            a pytemplate for frequency calculations of minima
        balsam_exe_settings : Dict[str, int]
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.
            balsam_exe_settings = {'num_nodes': 1,
                                   'ranks_per_node': 48,
                                   'threads_per_rank': 1}
        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            '/home/mgierad/espresso/pseudo'
        pseudopotentials : Dict[str, str]
            a dictionary with QE pseudopotentials for all species.
            e.g.
            dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                H='H.pbe-kjpaw_psl.1.0.0.UPF',
                O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                C='C.pbe-n-kjpaw_psl.1.0.0.UPF',
                )
        calc_keywords : Dict[str, str]
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}
        creation_dir : PosixPath
            a posix path to the working directory

        '''
        minima_vib_path = os.path.join(
            self.creation_dir, self.facetpath, 'minima_vib')
        os.makedirs(minima_vib_path, exist_ok=True)
        path_to_minima_species = os.path.join(self.minima_path, species)
        path_to_vib_species = os.path.join(minima_vib_path, species)
        os.makedirs(path_to_vib_species, exist_ok=True)

        traj_fname = os.path.join(path_to_vib_species, species + '.traj')

        self.create_minima_vib_xyz(
            path_to_minima_species, traj_fname)
        self.create_minima_vib_py_files(
            species, traj_fname, minima_vib_path, pytemplate,
            balsam_exe_settings, pseudo_dir, pseudopotentials,
            calc_keywords, creation_dir)

    def create_minima_vib_xyz(
            self,
            path_to_minima_species: str,
            traj_fname: str) -> None:
        ''' Create traj file with optimized minima for a vibrational frequency
            calculations

        Parameters
        ----------
        path_to_minima_species : str
            a path where all minima calculations for a given species are
        traj_fname : str
            a name of a new traj file

        '''
        path_to_min_ener_species = self.get_min_conformer(
            path_to_minima_species)
        shutil.copy2(path_to_min_ener_species, traj_fname)

    def get_min_conformer(
            self,
            path_to_minima_species: str) -> str:
        ''' Get path to the most stable conformer of a given species

        Parameters
        ----------
        path_to_minima_species : str
            a path where all minima calculations for a given species are

        Returns
        -------
        path_to_min_ener_species: str
            a path to the most stable conformer of a given species

        '''
        minima_energies = minimaVib.get_minima_energies(path_to_minima_species)
        path_to_min_ener_species = min(
            minima_energies, key=minima_energies.get)
        return path_to_min_ener_species

    @staticmethod
    def get_minima_energies(
            path_to_minima_species: str) -> Dict[str, float]:
        '''[summary]

        Parameters
        ----------
        path_to_minima_species : str
            a path where all minima calculations for a given species are

        Returns
        -------
        energies: Dict[str, float]
            a dict with path to a given minima and potential energy
            associated with it
        '''
        energies = {}
        final_xyz_files = Path(path_to_minima_species).glob('*traj')
        for xyz in final_xyz_files:
            species_atom = read(xyz)
            xyz = str(xyz)
            energies[xyz] = species_atom.get_potential_energy()
        return energies

    def create_minima_vib_py_files(
            self,
            species: str,
            traj_fname: str,
            minima_vib_path: str,
            pytemplate: str,
            balsam_exe_settings: Dict[str, int],
            pseudo_dir: str,
            pseudopotentials: Dict[str, str],
            calc_keywords: Dict[str, str],
            creation_dir: PosixPath,
            nimages: int = 30,
            n: int = 0) -> None:
        ''' Create a .py files for a vibrational frequiency calculations
            for a given species

        Parameters
        ----------
        species : str
            a chemical symbol of a reacting species
        traj_fname : str
            a name of a new traj file
        minima_vib_path : str
            a path to vibrational frequency calculations of minima,
            e.g. 'Cu_111/minima_vib'
        pytemplate : str
            a pytemplate for frequency calculations of minima
        balsam_exe_settings : Dict[str, int]
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.
            balsam_exe_settings = {'num_nodes': 1,
                                   'ranks_per_node': 48,
                                   'threads_per_rank': 1}
        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            '/home/mgierad/espresso/pseudo'
        pseudopotentials : Dict[str, str]
            a dictionary with QE pseudopotentials for all species.
            e.g.
            dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                H='H.pbe-kjpaw_psl.1.0.0.UPF',
                O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                C='C.pbe-n-kjpaw_psl.1.0.0.UPF',
                )
        calc_keywords : Dict[str, str]
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}
        creation_dir : PosixPath
            a posix path to the working directory
        nimages : int, optional
            how many strucutres to use to construct a trajectory visualizing
            oscilations, by default 30
        n : int, optional
            mode of oscilation, i.e. which vibration to analyze.
            0 is the first vibration, should be imaginary, by default 0

        '''
        with open(pytemplate, 'r') as f:
            pytemplate_txt = f.read()
            py_file_name = os.path.join(
                self.facetpath + '_' + species + '_vib.py')
            py_file = os.path.join(minima_vib_path, py_file_name)
            with open(py_file, 'w') as c:
                c.write(pytemplate_txt.format(
                    geom=traj_fname,
                    balsam_exe_settings=balsam_exe_settings,
                    creation_dir=creation_dir,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    calc_keywords=calc_keywords,
                    nimages=nimages,
                    n=n,
                ))

    def dependency_minima_vib(
            self,
            species: str) -> List[str]:
        ''' Get a list with all dependant species minima calculations that have
            to be finished before vibrational frequency for that species
            can start

        Parameters
        ----------
        species : str
            a chemical symbol of the adsorbate

        Returns
        -------
        dependancy_minima_list : List[str]
            a list with all minima calculations for that particular species,
            that have to be finished before starting vibrational frequency
            calculations

        '''
        dependancy_minima_list = []
        dependancy_minima = Path(self.minima_path).glob(
            '*_{}_*.py'.format(species))
        for pyscript in dependancy_minima:
            dependancy_minima_list.append(os.path.basename(pyscript))
        return dependancy_minima_list
