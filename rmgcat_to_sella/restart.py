import os
from ase.io import read, write
from ase.io.formats import UnknownFileTypeError
# from balsam.launcher.dag import BalsamJob

all_unfinished = {'05_Cu_100_set_up_TS_vib_CH_C+H.py': 'AWAITING_PARENTS', '04_Cu_111_set_up_TS_vib_OH_O+H.py': 'AWAITING_PARENTS', '03_Cu_100_run_TS_CH_C+H.py': 'AWAITING_PARENTS', '03_Cu_100_run_TS_OH_O+H.py': 'AWAITING_PARENTS', 'Cu_111_HO_03_relax.py': 'RUNNING', 'Cu_111_O_00_relax.py': 'RUNNING', '02_Cu_100_set_up_TS_with_xtb_OH_O+H.py': 'AWAITING_PARENTS', '02_Cu_111_set_up_TS_with_xtb_CH_C+H.py': 'AWAITING_PARENTS', '03_Cu_111_run_TS_CH_C+H.py': 'AWAITING_PARENTS', '04_Cu_111_set_up_TS_vib_CH_C+H.py': 'AWAITING_PARENTS', '03_Cu_111_run_TS_OH_O+H.py': 'AWAITING_PARENTS', '04_Cu_100_set_up_TS_vib_CH_C+H.py': 'AWAITING_PARENTS', '05_Cu_100_set_up_TS_vib_OH_O+H.py': 'AWAITING_PARENTS', 'Cu_100_CH_01_relax.py': 'RUNNING',
                  'Cu_100_HO_02_relax.py': 'RUNNING', 'Cu_100_O_01_relax.py': 'RUNNING', '02_Cu_111_set_up_TS_with_xtb_OH_O+H.py': 'AWAITING_PARENTS', '05_Cu_111_set_up_TS_vib_OH_O+H.py': 'AWAITING_PARENTS', 'Cu_100_HO_01_relax.py': 'RUNNING', 'Cu_100_C_01_relax.py': 'RUNNING', '02_Cu_100_set_up_TS_with_xtb_CH_C+H.py': 'AWAITING_PARENTS', '04_Cu_100_set_up_TS_vib_OH_O+H.py': 'AWAITING_PARENTS', 'Cu_111_C_03_relax.py': 'RUNNING', 'Cu_111_O_01_relax.py': 'RUNNING', 'Cu_111_O_02_relax.py': 'RUNNING', 'Cu_111_O_03_relax.py': 'RUNNING', 'Cu_111_H_00_relax.py': 'RUNNING', 'Cu_111_H_01_relax.py': 'RUNNING', 'Cu_111_H_02_relax.py': 'RUNNING', 'Cu_111_H_03_relax.py': 'RUNNING', '05_Cu_111_set_up_TS_vib_CH_C+H.py': 'AWAITING_PARENTS'}


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
        # all_unfinished = self.get_jobs_to_restart()
        unfinished_minima = []
        for key in all_unfinished.keys():
            if 'relax' in key:
                unfinished_minima.append(key)
        return unfinished_minima

    def prepare_minima_to_restart(self):
        unfinished_minima = self.get_minima_to_restart()
        for minimum in unfinished_minima:
            metal_symbol, facet, species, prefix, _ = minimum.split('_')
            facetpath = metal_symbol + '_' + facet
            path_to_species = os.path.join(
                facetpath, 'minima', species, prefix)
            new_xyz_fname = os.path.join(
                path_to_species + '.xyz')

            try:
                atom_traj_file = read(path_to_species + '.traj')
                write(new_xyz_fname, atom_traj_file)
            except UnknownFileTypeError:
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
