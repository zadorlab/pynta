#!/usr/bin/env python3
from collections import Counter


class Show():
    def __init__(self) -> None:
        from balsam.launcher.dag import BalsamJob
        # get all python (ASE) jobs
        self.ase_jobs = BalsamJob.objects.filter(
            application__contains='python')

    def how_many_still_running(self):
        ''' Show how many jobs required restart

        Returns
        -------
        len(running_jobs) : int
            a number of jobs that has to be restarted

        '''
        running_jobs = []
        for job in self.ase_jobs:
            if job.state != 'JOB_FINISHED':
                running_jobs.append(job)
        return len(running_jobs)

    def status(self) -> None:
        ''' Show info about the current status of the Balsam DB, i.e.
        How many jobs are running? How many already finished? etc...

        '''
        job_states = []
        for job in self.ase_jobs:
            job_states.append(job.state)
        current_state = Counter(job_states)
        for key, val in current_state.items():
            print('{:>17} : {:>4}'.format(key, val))

    def not_finished(self) -> None:
        ''' Show info about all jobs that did not finish:
        State, workflow name jobname

        '''
        not_finished = {}
        for job in self.ase_jobs:
            if job.state != 'JOB_FINISHED':
                # TODO improve because keys are overwritten
                key = 'Workflow : {:>30} Jobname : {:>50}'.format(
                    job.workflow, job.name)
                not_finished[key] = job.state
        # return not_finished
        for key, val in not_finished.items():
            print('{:>4} : {:>18}'.format(key, val))
