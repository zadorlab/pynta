from pynta.io import IO
from ase.io import read, write
from pathlib import Path, PosixPath
from typing import Tuple, Dict
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
            socket_calculator: str,
            rxn: Dict[str, str],
            pytemplate: str,
            balsam_exe_settings: Dict[str, int],
            calc_keywords: Dict[str, str],
            pseudopotentials: Dict[str, str],
            pseudo_dir: str) -> None:
        ''' Set up files for TSs vibration calculations

        Parameters
        ----------
        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            :literal:`*.yaml` file
        pytemplate : python script
            a template for TS_vib calculations
        balsam_exe_settings : Dict[str, int],
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.

            >>> balsam_exe_settings = {'num_nodes': 1,
                                    'ranks_per_node': 48,
                                    'threads_per_rank': 1}

        calc_keywords : Dict[str, str]
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            >>> calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                                'smearing':  'marzari-vanderbilt',
                                'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                                'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}

        pseudopotentials : Dict[str, str]
            a dictionary with QE pseudopotentials for all species.
            e.g.

            >>> dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                    H='H.pbe-kjpaw_psl.1.0.0.UPF',
                    O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                    C='C.pbe-n-kjpaw_psl.1.0.0.UPF')

        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            ``'/home/mgierad/espresso/pseudo'``

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
                socket_calculator,
                pytemplate,
                geom,
                py_fname,
                balsam_exe_settings,
                calc_keywords,
                pseudopotentials,
                pseudo_dir)

    def create_ts_vib_py_files(
            self,
            socket_calculator: str,
            pytemplate: str,
            geom: str,
            py_fname: str,
            balsam_exe_settings: Dict[str, int],
            calc_keywords: Dict[str, str],
            pseudopotentials: Dict[str, str],
            pseudo_dir: str,
            nimages: int = 30,
            n: int = 0) -> None:
        ''' Create job submission .py files for frequency calculation of TSs.

        Parameters
        ----------
        pytemplate : python script
            a template for TS_vib calculations
        py_fname : str
            path to the the .py (including the .py file name) file that is
            about to be created
        balsam_exe_settings : Dict[str, int]
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.

            >>> balsam_exe_settings = {'num_nodes': 1,
                                    'ranks_per_node': 48,
                                    'threads_per_rank': 1}

        calc_keywords : Dict[str, str]
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            >>> calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                                'smearing':  'marzari-vanderbilt',
                                'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                                'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}

        pseudopotentials : Dict[str, str]
            a dictionary with QE pseudopotentials for all species.
            e.g.

            >>> dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                    H='H.pbe-kjpaw_psl.1.0.0.UPF',
                    O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                    C='C.pbe-n-kjpaw_psl.1.0.0.UPF')

        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            ``'/home/mgierad/espresso/pseudo'``
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
                socket_calculator=socket_calculator,
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
            socket_calculator: str,
            rxn: Dict[str, str],
            pytemplate: str,
            balsam_exe_settings: Dict[str, int],
            calc_keywords: Dict[str, str],
            pseudopotentials: Dict[str, str],
            pseudo_dir: str) -> None:
        ''' Create files for after_TS calculations - to verify TS structures
        and get corresponding reactant and product minima

        Parameters
        ----------
        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            :literal:`*.yaml` file
        pytemplate : python script
            a template for after_TS calculations
        balsam_exe_settings : Dict[str, int]
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.

            >>> balsam_exe_settings = {'num_nodes': 1,
                                    'ranks_per_node': 48,
                                    'threads_per_rank': 1}

        calc_keywords : Dict[str, str]
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            >>> calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                                'smearing':  'marzari-vanderbilt',
                                'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                                'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}

        pseudopotentials : Dict[str, str]
            a dictionary with QE pseudopotentials for all species.
            e.g.

            >>> dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                    H='H.pbe-kjpaw_psl.1.0.0.UPF',
                    O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                    C='C.pbe-n-kjpaw_psl.1.0.0.UPF')

        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            ``'/home/mgierad/espresso/pseudo'``

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
                socket_calculator,
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
        ''' Get forward and reverse :literal:`*.xyz` file by nudging TS
        towards imaginary mode of oscilations summary.

        Parameters
        ----------
        vib_traj : str
            a path to trajectory file with optimized TS
        fname_forward : str
            a path for forward calculations
        fname_reverse : str
            a path for reverse calculations
        nimages : int, optional
            how many strucutres to use to construct a trajectory visualizing
            oscilations, by default 30

        Raises
        ------
        ValueError
            raised if there are more than one :literal:`*.traj` files
            collecting vibration trajectory

        '''

        # vib_.traj file contains nimages structures visualizing imaginary
        # frequency mode. 0 and nimages/2 are the same structures - no
        # displacement. Max and min displacement is defined as
        # floor(nimages/4) and nimages - floor(nimages/4), respectively.
        # for inmages = 16 it would be like this
        # start...max...start...min...
        index_forward = int(floor(nimages / 4))
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
            socket_calculator: str,
            pytemplate: str,
            fname_forward: str,
            fname_reverse: str,
            balsam_exe_settings: Dict[str, int],
            calc_keywords: Dict[str, str],
            pseudopotentials: Dict[str, str],
            pseudo_dir: str) -> None:
        ''' Create job submission files for minimization displaced structures
        after TS.

        Parameters
        ----------
        pytemplate : python script
            a template for after_TS calculations
        fname_forward : str
            a path for forward calculations
        fname_reverse : str
            a path for reverse calculations
        balsam_exe_settings : Dict[str, int]
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.

            >>> balsam_exe_settings = {'num_nodes': 1,
                                    'ranks_per_node': 48,
                                    'threads_per_rank': 1}

        calc_keywords : Dict[str, str]
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            >>> calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                                'smearing':  'marzari-vanderbilt',
                                'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                                'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}

        pseudopotentials : Dict[str, str]
            a dictionary with QE pseudopotentials for all species.
            e.g.

            >>> dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                    H='H.pbe-kjpaw_psl.1.0.0.UPF',
                    O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                    C='C.pbe-n-kjpaw_psl.1.0.0.UPF')

        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            ``'/home/mgierad/espresso/pseudo'``

        '''
        with open(pytemplate, 'r') as f:
            pytemplate = f.read()

        for fn in [fname_forward, fname_reverse]:
            tmp, geom = os.path.split(fn)
            fname = os.path.join(os.path.split(tmp)[0], geom + '.py')
            with open(fname, 'w') as f:
                f.write(pytemplate.format(
                    socket_calculator=socket_calculator,
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
        reverse structure.

        '''
        all_rxn_names = IO().get_list_all_rxns_names(self.yamlfile)
        for rxn_name in all_rxn_names:
            ts_dist_dict = self.get_ts_dist(rxn_name)
            f_dist_dict, r_dist_dict = self.get_forward_and_reverse_dist(
                rxn_name)
            AfterTS.print_table(ts_dist_dict, f_dist_dict, r_dist_dict)

    @staticmethod
    def print_table(
            ts_dist_dict: Dict[str, float],
            forward_dist_dict: Dict[str, float],
            reverse_dist_dict: Dict[str, float]) -> None:
        ''' Print information about bond distances (reacting atoms only) for TS,
        forward and reverse :literal:`*.xyz` files

        Parameters
        ----------
        ts_dist_dict : Dict[str:float]
            a dictionary with keys being TS file names while values are
            bond distance beteween of reacting species
        f_dist_dict : Dict[str:float]
            a dictionary with keys being paths to forward files while values
            are bond distance beteween of reacting species
        r_dist_dict : Dict[str:float]
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
        in TS structure.

        Parameters
        ----------
        rxn_name : str
            a name of the reaction in the following format,
            e.g.:
            ``'OH_H+O'``

        Returns
        -------
        ts_dist_dict : Dict[str:float]
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
            ```'OH_H+O'```

        Returns
        -------
        f_dist_dict : Dict[str:float]
            a dictionary with keys being paths to forward files while values
            are bond distance beteween of reacting species
        r_dist_dict : Dict[str:float]
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
            socket_calculator: str,
            facetpath: str,
            yamlfile_path: PosixPath,
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
        facetpath : str
            a path to the workflow's main dir
            e.g. ``'Cu_111'``
        yamlfile_path : posix
            a path to .yaml file with all reactions e.g. 'reactions.yaml'
        pytemplate : str
            a pytemplate for frequency calculations of minima
        balsam_exe_settings : Dict[str, int]
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.

            >>> balsam_exe_settings = {'num_nodes': 1,
                                    'ranks_per_node': 48,
                                    'threads_per_rank': 1}

        calc_keywords : Dict[str, str]
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            >>> calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                                'smearing':  'marzari-vanderbilt',
                                'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                                'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}

        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            ``'/home/mgierad/espresso/pseudo'``

        pseudopotentials : Dict[str, str]
            a dictionary with QE pseudopotentials for all species.
            e.g.

            >>> dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                    H='H.pbe-kjpaw_psl.1.0.0.UPF',
                    O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                    C='C.pbe-n-kjpaw_psl.1.0.0.UPF')

        calc_keywords : Dict[str, str]
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            >>> calc_keywords = {'kpts': (3, 3, 1),
                    'occupations': 'smearing',
                    'smearing': 'marzari-vanderbilt',
                    'degauss': 0.01,
                    'ecutwfc': 40,
                    'nosym': True,
                    'conv_thr': 1e-11,
                    'mixing_mode':'local-TF'}

        creation_dir : PosixPath
            a posix path to the working directory

        '''
        minima_vib_path = os.path.join(
            self.creation_dir, self.facetpath, 'minima_vib')
        os.makedirs(minima_vib_path, exist_ok=True)

        unique_adsorbates_prefixes = IO().get_unique_adsorbates_prefixes(
            facetpath, yamlfile_path, creation_dir)
        for adsorbate, unique_prefixes in unique_adsorbates_prefixes.items():
            for prefix in unique_prefixes:
                path_to_minimum_traj = os.path.join(
                    self.minima_path, adsorbate, prefix + '.traj')

                path_to_vib_species = os.path.join(
                    minima_vib_path, adsorbate, prefix)
                os.makedirs(path_to_vib_species, exist_ok=True)

                traj_to_start_vib = os.path.join(
                    path_to_vib_species,
                    '{}_{}.traj'.format(prefix, adsorbate))
                shutil.copy2(path_to_minimum_traj, traj_to_start_vib)

                self.create_minima_vib_py_files(
                    socket_calculator, adsorbate, prefix, traj_to_start_vib,
                    minima_vib_path, pytemplate, balsam_exe_settings,
                    pseudo_dir, pseudopotentials, calc_keywords, creation_dir)

    def create_minima_vib_py_files(
            self,
            socket_calculator: str,
            adsorbate: str,
            prefix: str,
            traj_to_start_vib: str,
            minima_vib_path: str,
            pytemplate: str,
            balsam_exe_settings: Dict[str, int],
            pseudo_dir: str,
            pseudopotentials: Dict[str, str],
            calc_keywords: Dict[str, str],
            creation_dir: PosixPath,
            nimages: int = 30,
            n: int = 0) -> None:
        ''' Create a :literal:`*.py` files for a vibrational frequiency
        calculations for a given species

        Parameters
        ----------
        adsorbate : str
            a chemical symbol of a reacting species
        prefix : str
            a prefix for unique minimum for which vibrational frequencies
            will be calculated
        traj_fname : str
            a name of a new traj file
        minima_vib_path : str
            a path to vibrational frequency calculations of minima,
            e.g. ``'Cu_111/minima_vib'``
        pytemplate : str
            a pytemplate for frequency calculations of minima
        balsam_exe_settings : Dict[str, int]
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.

            >>> balsam_exe_settings = {'num_nodes': 1,
                                    'ranks_per_node': 48,
                                    'threads_per_rank': 1}

        calc_keywords : Dict[str, str]
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            >>> calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                                'smearing':  'marzari-vanderbilt',
                                'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                                'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}

        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            ``'/home/mgierad/espresso/pseudo'``

        pseudopotentials : Dict[str, str]
            a dictionary with QE pseudopotentials for all species.
            e.g.

            >>> dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                    H='H.pbe-kjpaw_psl.1.0.0.UPF',
                    O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                    C='C.pbe-n-kjpaw_psl.1.0.0.UPF')

        calc_keywords : Dict[str, str]
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            >>> calc_keywords = {'kpts': (3, 3, 1),
                    'occupations': 'smearing',
                    'smearing': 'marzari-vanderbilt',
                    'degauss': 0.01,
                    'ecutwfc': 40,
                    'nosym': True,
                    'conv_thr': 1e-11,
                    'mixing_mode':'local-TF'}

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
            py_file_name = os.path.join('{}_{}_{}_vib.py'.format(
                self.facetpath, prefix, adsorbate))
            py_file = os.path.join(minima_vib_path, py_file_name)
            with open(py_file, 'w') as c:
                c.write(pytemplate_txt.format(
                    socket_calculator=socket_calculator,
                    geom=traj_to_start_vib,
                    balsam_exe_settings=balsam_exe_settings,
                    creation_dir=creation_dir,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    calc_keywords=calc_keywords,
                    nimages=nimages,
                    n=n,
                ))
