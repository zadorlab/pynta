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
        job_states = []
        for job in self.ase.jobs:
            job_states.append(job.state)
        current_state = Counter(job_states)
        for key, val in current_state.items():
            print('{:>15} : {:>4}'.format(key, val))
