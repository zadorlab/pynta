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

            fname_forward = os.path.join(after_ts_dir, 'forward.xyz')
            fname_reverse = os.path.join(after_ts_dir, 'reverse.xyz')

            self.get_forward_and_reverse(traj, fname_forward, index_forward)
            self.get_forward_and_reverse(traj, fname_reverse, index_reverse)

        # for ts in os.listdir(ts_estimate_unique_dir):
        #     ts_path = os.path.join(ts_estimate_unique_dir, ts)
        #     if os.path.isdir(ts_path):
        #         prefix = ts_path

    def get_forward_and_reverse(self, traj, fname, index, nimages=30):
        ''' Get forward and reverse .xyz file by nudging TS towards imaginary
        mode of oscilations '''

        index_forward = int(floor(nimages/4))
        index_reverse = int(nimages - index_forward)

        traj_atom = read(traj)
        traj_atom.calc(EMT())
        vib = Vibrations(traj_atom)
        vib.run()
        vib.summary()
        vib.clean()
        vib.write_mode(0, nimages=nimages)
        write(fname, read(traj, index=index))
