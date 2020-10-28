import os
from balsam.launcher.dag import BalsamJob
from pathlib import Path


class LowLevelRestart():
    def __init__(self):
        self.current_dir = os.getcwd()
        # get all python (ASE) jobs
        self.ase_jobs = BalsamJob.objects.filter(
            application__contains='python')

    def jobs_to_restart(self):
        ''' Get a dictionary with all ase_jobs that did not finish '''
        unfinished = {}
        for job in self.ase_jobs:
            if job.state != 'JOB_FINISHED':
                unfinished[job.name] = job.state
        return unfinished

    # def get_minima_to_restart(self):
    #     minima_out_files = Path(self.current_dir).glob('**/*.out')
    #     for mininma_out_file in minima_out_files:
    #         print(mininma_out_file)


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
