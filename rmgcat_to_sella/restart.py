#!/usr/bin/env python3
from balsam.launcher.dag import BalsamJob
from collections import Counter


class Restart():
    def __init__(self):
        # get all python (ASE) jobs
        self.ase_jobs = BalsamJob.objects.filter(
            application__contains='python')

    def restart(self):
        # remove all balsam calculator objects
        BalsamJob.objects.filter(name__contains='balsam',
                                 workflow='QE_Socket').delete()

        # update state of every not finished jobs to 'READY'
        for job in self.ase_jobs:
            if job.state not in ['JOB_FINISHED', 'AWAITING_PARENTS']:
                job.state = 'READY'
                job.save()

    def how_many_still_running(self):
        running_jobs = []
        for job in self.ase_jobs:
            if job.state != 'JOB_FINISHED':
                running_jobs.append(job)
        return len(running_jobs)

    def describe(self):
        ''' Print info about the current status of the Balsam DB, i.e.
            How many jobs are running? How many already finished? etc...

        '''
        job_states = []
        for job in self.ase_jobs:
            job_states.append(job.state)
        current_state = Counter(job_states)
        for key, val in current_state.items():
            print('{:>15} : {:>4}'.format(key, val))

    def not_finished(self):
        not_finished = {}
        for job in self.ase_jobs:
            if job.state != 'JOB_FINISHED':
                not_finished[job.state] = 'Workflow : {}. Jobname : {}'.format(
                    job.workflow, job.name)
        # return not_finished
        for key, val in not_finished.items():
            print('{:>15} : {:>4}'.format(key, val))
