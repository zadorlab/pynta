#!/usr/bin/env python3
from balsam.launcher.dag import BalsamJob


def restart():
    # remove all balsam calculator objects
    BalsamJob.objects.filter(name__contains='balsam',
                             workflow='QE_Socket').delete()

    # get all python (ASE) jobs
    ase_jobs = BalsamJob.objects.filter(application__contains='python')

    # update state of every not finished jobs to 'READY'
    for job in ase_jobs:
        if job.state not in ['JOB_FINISHED', 'AWAITING_PARENTS']:
            job.state = 'READY'
            job.save()
