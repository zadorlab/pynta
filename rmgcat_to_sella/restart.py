#!/usr/bin/env python3
from balsam.launcher.dag import BalsamJob


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
