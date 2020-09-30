from rmgcat_to_sella.io import IO

from ase.calculators.emt import EMT
from ase.vibrations import Vibrations
from ase.io import read, write

from pathlib import Path

from numpy import floor

import os


class AfterTS():
    def __init__(self, facetpath, yamlfile):
        self.facetpath = facetpath
        self.yamlfile = yamlfile

    def prepare_all(self):
        all_rxn_names = IO().get_list_all_rxns_names(self.yamlfile)
        for rxn_name in all_rxn_names:
            self.prepare_after_ts(rxn_name)

    def prepare_after_ts(self, rxn_name):
        ts_estimate_unique_dir = os.path.join(
            self.facetpath, rxn_name, 'TS_estimate_unique')
        traj_files = Path(ts_estimate_unique_dir).glob('**/*traj')
        for traj in traj_files:
            traj = str(traj)
            prefix = traj.split('/')[-2]

            after_ts_dir = os.path.join(
                self.facetpath, rxn_name, 'after_TS', prefix)
            os.makedirs(after_ts_dir, exist_ok=True)

            fname_forward = os.path.join(after_ts_dir, 'forward_in')
            fname_reverse = os.path.join(after_ts_dir, 'reverse_in')

            self.get_forward_and_reverse(traj, fname_forward, fname_reverse)

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

        traj_atom = read(traj)
        traj_atom.calc = EMT()
        vib = Vibrations(traj_atom)
        vib.run()
        vib.summary()
        vib.clean()
        vib.write_mode(n=n, nimages=nimages)

        write(fname_forward + '.xyz', read(traj, index=index_forward))
        write(fname_forward + '.png', read(traj, index=index_forward))
        write(fname_reverse + '.xyz', read(traj, index=index_reverse))
        write(fname_reverse + '.png', read(traj, index=index_reverse))

    def get_all_distances(self):
        ''' Get distances between reacting species for ts, forward and
        reverse structure '''
        all_rxn_names = IO().get_list_all_rxns_names(self.yamlfile)
        for rxn_name in all_rxn_names:
            ts_dist = self.get_ts_dist(rxn_name)
            print(ts_dist)
            # forward_dist, reverse_dist = self.get_forward_and_reverse_dist()
            # self.get_forward_and_reverse_dist(rxn_name)

    def get_ts_dist(self, rxn_name):
        ts_dist_dict = {}

        ts_estimate_unique_dir = os.path.join(
            self.facetpath, rxn_name, 'TS_estimate_unique')

        traj_files = Path(ts_estimate_unique_dir).glob('**/*traj')

        for traj in sorted(traj_files):
            traj = str(traj)
            traj_atom = read(traj)[27:]
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
            if 'reverse' in xyz:
                xyz_atom = read(xyz)[27:]
                dist = xyz_atom.get_distance(0, 1)
                print(dist)
