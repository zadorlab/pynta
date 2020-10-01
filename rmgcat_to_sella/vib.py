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

    def prepare_opt_after_ts(
            self,
            rxn,
            pytemplate,
            balsam_exe_settings,
            calc_keywords,
            creation_dir,
            pseudopotentials,
            pseudo_dir):

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

            # self.get_forward_and_reverse(traj, fname_forward, fname_reverse)
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
            after TS '''

        with open(pytemplate, 'r') as f:
            pytemplate = f.read()
        for fn in [fname_forward, fname_reverse]:
            fname = fn + '.py'
            with open(fname, 'w') as f:
                f.write(pytemplate.format(
                    geom=fn+'xyz',
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
        mode of oscilations '''

        index_forward = int(floor(nimages/4))
        index_reverse = int(nimages - index_forward)

        vib_traj_path, _ = os.path.split(fname_forward)

        traj_atom = read(traj)
        traj_atom.calc = EMT()
        name_vib_files = os.path.join(vib_traj_path, 'vib')
        vib = Vibrations(traj_atom, name=name_vib_files)
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
        reverse structure '''
        all_rxn_names = IO().get_list_all_rxns_names(self.yamlfile)
        for rxn_name in all_rxn_names:
            ts_dist_dict = self.get_ts_dist(rxn_name)
            forward_dist_dict, reverse_dist_dict = self.get_forward_and_reverse_dist(
                rxn_name)
            self.print_table(ts_dist_dict, forward_dist_dict,
                             reverse_dist_dict)

    def print_table(self, ts_dist_dict, forward_dist_dict, reverse_dist_dict):
        keys = ts_dist_dict.keys()
        ts_val = ts_dist_dict.values()
        f_val = forward_dist_dict.values()
        r_val = reverse_dist_dict.values()

        print('TS \t \t TS_dist \t F_dist \t R_dist')
        for key, ts, f, r in zip(keys, ts_val, f_val, r_val):

            print('{} \t {:2f} \t {:2f} \t {:2f}'.format(key, ts, f, r))

    def get_ts_dist(self, rxn_name):
        ts_dist_dict = {}

        ts_estimate_unique_dir = os.path.join(
            self.facetpath, rxn_name, 'TS_estimate_unique')

        traj_files = Path(ts_estimate_unique_dir).glob('**/*traj')

        for traj in sorted(traj_files):
            traj = str(traj)
            traj_atom = read(traj)[self.nslab:]
            dist = traj_atom.get_distance(0, 1)
            traj_fname = traj.split('/')[-1]
            ts_dist_dict[traj_fname] = dist
        return ts_dist_dict

    def get_forward_and_reverse_dist(self, rxn_name):
        forward_dist_dict = {}
        reverse_dist_dict = {}

        after_ts_dir = os.path.join(
            self.facetpath, rxn_name, 'after_TS')

        xyz_files = Path(after_ts_dir).glob('**/*xyz')

        for xyz in sorted(xyz_files):
            xyz = str(xyz)
            if 'forward' in xyz:
                xyz_atom = read(xyz)[self.nslab:]
                dist = xyz_atom.get_distance(0, 1)
                forward_dist_dict[xyz] = dist
            elif 'reverse' in xyz:
                xyz_atom = read(xyz)[self.nslab:]
                dist = xyz_atom.get_distance(0, 1)
                reverse_dist_dict[xyz] = dist
        return forward_dist_dict, reverse_dist_dict
