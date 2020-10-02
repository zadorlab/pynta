from rmgcat_to_sella.io import IO

from ase.calculators.emt import EMT
from ase.vibrations import Vibrations
from ase.io import read, write

from pathlib import Path

from numpy import floor

import os


class AfterTS():
    def __init__(
            self,
            facetpath,
            yamlfile,
            slab,
            repeats):

        self.facetpath = facetpath
        self.yamlfile = yamlfile
        self.slab = slab
        self.repeats = repeats
        self.nslab = len(read(self.slab)*self.repeats)
        self.io = IO()

    def set_up_ts_vib(
            self,
            rxn,
            pytemplate):

        rxn_name = self.io.get_rxn_name(rxn)
        ts_estimate_unique_dir = os.path.join(
            self.facetpath, rxn_name, 'TS_estimate_unique')
        traj_files = Path(ts_estimate_unique_dir).glob('**/*traj')
        for traj in traj_files:
            traj = str(traj)
            prefix = traj.split('/')[-2]

            ts_vib_dir = os.path.join(
                self.facetpath, rxn_name, 'TS_estimate_unique_vib', prefix)
            os.makedirs(ts_vib_dir, exist_ok=True)

            fname_vib = os.path.join(
                ts_vib_dir, prefix + '_' + rxn_name + '_vib')

    def prepare_opt_after_ts(
            self,
            rxn,
            pytemplate,
            balsam_exe_settings,
            calc_keywords,
            creation_dir,
            pseudopotentials,
            pseudo_dir):
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
        rxn_name = self.io.get_rxn_name(rxn)
        ts_estimate_unique_dir = os.path.join(
            self.facetpath, rxn_name, 'TS_estimate_unique')
        traj_files = Path(ts_estimate_unique_dir).glob('**/*traj')
        for traj in traj_files:
            traj = str(traj)
            prefix = traj.split('/')[-2]

            after_ts_dir = os.path.join(
                self.facetpath, rxn_name, 'after_TS', prefix)
            os.makedirs(after_ts_dir, exist_ok=True)

            fname_forward = os.path.join(
                after_ts_dir, prefix + '_' + rxn_name + '_after_ts_f')
            fname_reverse = os.path.join(
                after_ts_dir, prefix + '_' + rxn_name + '_after_ts_r')

            self.get_forward_and_reverse(traj, fname_forward, fname_reverse)
            self.create_after_ts_py_files(
                pytemplate, fname_forward, fname_reverse, balsam_exe_settings,
                calc_keywords, creation_dir, pseudopotentials, pseudo_dir)

    def create_after_ts_py_files(
            self,
            pytemplate,
            fname_forward,
            fname_reverse,
            balsam_exe_settings,
            calc_keywords,
            creation_dir,
            pseudopotentials,
            pseudo_dir):
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
                    creation_dir=creation_dir,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir
                ))

    def get_forward_and_reverse(
            self,
            traj,
            fname_forward,
            fname_reverse,
            n=0,
            nimages=30):
        ''' Get forward and reverse .xyz file by nudging TS towards imaginary
        mode of oscilations summary

        Parameters
        ----------
        traj : str
            a path to trajectory file with optimized TS
        fname_forward : str
            a path for forward calculations
        fname_reverse : str
            a path for reverse calculations
            [description]
        n : int, optional
            mode of oscilation 0 is the first (imaginary), by default 0
        nimages : int, optional
            how many strucutres to use to construct a trajectory visualizing
            oscilations, by default 30

        Raises
        ------
        ValueError
            raised if there are more than one *traj file visualizing vibrations

        '''
        index_forward = int(floor(nimages/4))
        index_reverse = int(nimages - index_forward)

        vib_traj_path, _ = os.path.split(fname_forward)

        traj_atoms = read(traj)
        traj_atoms.calc = EMT()
        # vibrate only adsorbates, not surface
        indices = [atom.index for atom in traj_atoms if atom.symbol != 'Cu']
        name_vib_files = os.path.join(vib_traj_path, 'vib')
        vib = Vibrations(traj_atoms, indices=indices, name=name_vib_files)
        vib.run()
        vib.summary()
        vib.clean()
        vib.write_mode(n=n, nimages=nimages)
        traj_list = []
        for vtraj in os.listdir(vib_traj_path):
            if vtraj.endswith('traj'):
                vib_traj = os.path.join(vib_traj_path, vtraj)
                traj_list.append(vib_traj)
        if len(traj_list) > 1:
            print('!!!')
            print('Only one vibrational trajectory is allowed.')
            print('!!!')
            raise ValueError

        write(fname_forward + '.xyz', read(vib_traj, index=index_forward))
        write(fname_forward + '.png', read(vib_traj, index=index_forward))
        write(fname_reverse + '.xyz', read(vib_traj, index=index_reverse))
        write(fname_reverse + '.png', read(vib_traj, index=index_reverse))

    def get_all_distances(self):
        ''' Get distances between reacting species for ts, forward and
        reverse structure

        '''
        all_rxn_names = IO().get_list_all_rxns_names(self.yamlfile)
        for rxn_name in all_rxn_names:
            ts_dist_dict = self.get_ts_dist(rxn_name)
            f_dist_dict, r_dist_dict = self.get_forward_and_reverse_dist(
                rxn_name)
            self.print_table(ts_dist_dict, f_dist_dict, r_dist_dict)

    def print_table(self, ts_dist_dict, forward_dist_dict, reverse_dist_dict):
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

    def get_ts_dist(self, rxn_name):
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

    def get_forward_and_reverse_dist(self, rxn_name):
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
