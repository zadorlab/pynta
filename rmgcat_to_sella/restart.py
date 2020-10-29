import os
from sys import path
from ase.io import read, write
from ase.io.formats import UnknownFileTypeError
from numpy.core.fromnumeric import prod
# from balsam.launcher.dag import BalsamJob

all_unfinished = {'05_Cu_100_set_up_TS_vib_CH_C+H.py': 'AWAITING_PARENTS', '04_Cu_111_set_up_TS_vib_OH_O+H.py': 'AWAITING_PARENTS', '03_Cu_100_run_TS_CH_C+H.py': 'AWAITING_PARENTS', '03_Cu_100_run_TS_OH_O+H.py': 'AWAITING_PARENTS', 'Cu_111_HO_03_relax.py': 'RUNNING', 'Cu_111_O_00_relax.py': 'RUNNING', '02_Cu_100_set_up_TS_with_xtb_OH_O+H.py': 'AWAITING_PARENTS', '02_Cu_111_set_up_TS_with_xtb_CH_C+H.py': 'AWAITING_PARENTS', '03_Cu_111_run_TS_CH_C+H.py': 'AWAITING_PARENTS', '04_Cu_111_set_up_TS_vib_CH_C+H.py': 'AWAITING_PARENTS', '03_Cu_111_run_TS_OH_O+H.py': 'AWAITING_PARENTS', '04_Cu_100_set_up_TS_vib_CH_C+H.py': 'AWAITING_PARENTS', '05_Cu_100_set_up_TS_vib_OH_O+H.py': 'AWAITING_PARENTS', 'Cu_100_CH_01_relax.py': 'RUNNING',
                  'Cu_100_HO_02_relax.py': 'RUNNING', 'Cu_100_O_01_relax.py': 'RUNNING', '02_Cu_111_set_up_TS_with_xtb_OH_O+H.py': 'AWAITING_PARENTS', '05_Cu_111_set_up_TS_vib_OH_O+H.py': 'AWAITING_PARENTS', 'Cu_100_HO_01_relax.py': 'RUNNING', 'Cu_100_C_01_relax.py': 'RUNNING', '02_Cu_100_set_up_TS_with_xtb_CH_C+H.py': 'AWAITING_PARENTS', '04_Cu_100_set_up_TS_vib_OH_O+H.py': 'AWAITING_PARENTS', 'Cu_111_C_03_relax.py': 'RUNNING', 'Cu_111_O_01_relax.py': 'RUNNING', 'Cu_111_O_02_relax.py': 'RUNNING', 'Cu_111_O_03_relax.py': 'RUNNING', 'Cu_111_H_00_relax.py': 'RUNNING', 'Cu_111_H_01_relax.py': 'RUNNING', 'Cu_111_H_02_relax.py': 'RUNNING', 'Cu_111_H_03_relax.py': 'RUNNING', '05_Cu_111_set_up_TS_vib_CH_C+H.py': 'AWAITING_PARENTS', '00_Cu_111_OH_O+H_ts.py': 'RUNNING'}


class LowLevelRestart():
    def __init__(self):
        self.current_dir = os.getcwd()
        # get all python (ASE) jobs
        # self.ase_jobs = BalsamJob.objects.filter(
        #     application__contains='python')

    def get_jobs_to_restart(self):
        ''' Get a dictionary with all ase_jobs that did not finish '''

        all_unfinished = {}
        for job in self.ase_jobs:
            if job.state != 'JOB_FINISHED':
                all_unfinished[job.name] = job.state
        return all_unfinished

    def get_minima_to_restart(self):
        ''' Go through all_unfinished dictionary and create a list with all
            minima ase jobs that started but are not yet finished and require
            to be restarted.

        Returns
        -------
        unfinished_minima : list(str)
            a list with all *py files for minima optimization which
            did not finished

        '''
        # all_unfinished = self.get_jobs_to_restart()
        unfinished_minima = []
        for key in all_unfinished.keys():
            if 'relax' in key:
                unfinished_minima.append(key)
        return unfinished_minima

    def get_tss_to_restart(self):
        ''' Go through all_unfinished dictionary and create a list with all
            TS ase jobs that started but are not yet finished and require to be
            restarted.

        Returns
        -------
        unfinished_tss : list(str)
            a list with all *py files for started but yet finished saddle
            point optimizations

        '''
        # all_unfinished = self.get_jobs_to_restart()
        unfinished_tss = []
        for key, value in all_unfinished.items():
            if 'ts' in key and 'AWAITING_PARENTS' not in value:
                unfinished_tss.append(key)
        return unfinished_tss

    def prepare_minima_to_restart(self):
        ''' If there is at least one optimization step in a .traj file
            for a given minima, this method will create a new *.xyz file for
            that job from the last converged opt step.
            So, when HighLevelRestart is executed, the optimization can resume
            from the last converged step, taking advantade of previous calcs.

        '''
        unfinished_minima = self.get_minima_to_restart()
        for minimum in unfinished_minima:
            metal_symbol, facet, species, prefix, _ = minimum.split('_')
            facetpath = metal_symbol + '_' + facet
            path_to_species = os.path.join(
                facetpath, 'minima', species, prefix)

            try:
                # try convert last step in *traj file to a new .xyz file
                write(path_to_species + '.xyz', path_to_species + '.traj')
            except UnknownFileTypeError:
                # continue if *traj file is empty
                # hard HighLevelRestart is required
                continue

    def prepare_ts_to_restart(self):
        unfinished_tss = self.get_tss_to_restart()
        for ts in unfinished_tss:
            prefix, metal_symbol, facet, reactant, product, _ = ts.split(
                '_')
            facetpath = metal_symbol + '_' + facet
            rxn_name = reactant + '_' + product
            ts_name = os.path.join(
                prefix + '_' + facetpath + '_' + rxn_name + '_ts')
            path_to_ts = os.path.join(
                facetpath, rxn_name, 'TS_estimate_unique', prefix, ts_name)

            try:
                # try convert last step in *traj file to a new .xyz file
                write(path_to_ts + '.xyz', read(path_to_ts + '.traj'))
            except UnknownFileTypeError:
                # continue if *traj file is empty
                # hard HighLevelRestart is required
                continue


class HighLevelRestart():
    def __init__(self):
        # get all python (ASE) jobs
        self.ase_jobs = BalsamJob.objects.filter(
            application__contains='python')

    def restart(self):
        ''' Prepare all unfinished jobs to restart

        '''
        # remove all balsam calculator objects
        BalsamJob.objects.filter(name__contains='balsam',
                                 workflow='QE_Socket').delete()

        # update state of every not finished jobs to 'READY'
        for job in self.ase_jobs:
            if job.state not in ['JOB_FINISHED', 'AWAITING_PARENTS']:
                job.state = 'READY'
                job.save()
