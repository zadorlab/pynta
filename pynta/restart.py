import os
from ase.io import read, write
from ase.io.formats import UnknownFileTypeError
from typing import Dict, List


class LowLevelRestart():
    ''' Low level restart means that each quantum chemistry job that can be
        restarted will be restarted. These are:
            * adsorbates minimization
            * TSs optimization
            * optimization of reactants and products as from TS
            (IRC alternative)

        For those jobs, once restarted, calculations will start from the last
        converged optimization step.

        For other types of jobs, including:
            * vibrational frequencies calculations of adsorbates and TSs
            * xTB + penalty function minimization to get TS guesses

        low level restart is not possible, so they will start from scratch,
        during high level restart.

    '''

    def __init__(self) -> None:
        from balsam.launcher.dag import BalsamJob
        self.current_dir = os.getcwd()
        # get all python(ASE) jobs
        self.ase_jobs = BalsamJob.objects.filter(
            application__contains='python')

    def restart(self) -> None:
        self.prepare_minima_to_restart()
        self.prepare_ts_to_restart()
        self.prepare_after_ts_to_restart()

    def get_jobs_to_restart(self) -> Dict[str, str]:
        ''' Get a dictionary with all ase_jobs that did not finish

        Returns:
        ________
        all_unfinished: Dict[str, str]
            a dictionary with all unfinished jobs, e.g.
            all_unfinished = {Cu_111_C_01_relax.py: 'RUNNING'}

        '''
        all_unfinished = {}
        for job in self.ase_jobs:
            if job.state != 'JOB_FINISHED':
                all_unfinished[job.name] = job.state
        return all_unfinished

    def get_minima_to_restart(self) -> List[str]:
        ''' Go through all_unfinished dictionary and create a list with all
            minima ase jobs that started but are not yet finished and require
            to be restarted.

        Returns
        -------
        unfinished_minima : list(str)
            a list with all *py files for minima optimization which
            did not finished

        '''
        all_unfinished = self.get_jobs_to_restart()
        unfinished_minima = []
        for key in all_unfinished.keys():
            if 'relax' in key:
                unfinished_minima.append(key)
        return unfinished_minima

    def get_tss_to_restart(self) -> List[str]:
        ''' Go through all_unfinished dictionary and create a list with all
            TS ase jobs that started but are not yet finished and require to be
            restarted.

        Returns
        -------
        unfinished_tss : list(str)
            a list with all *py files for started but yet finished saddle
            point optimizations

        '''
        all_unfinished = self.get_jobs_to_restart()
        unfinished_tss = []
        for key, value in all_unfinished.items():
            only_ts = 'ts' in key and all(
                [k not in key for k in ['ts_vib', 'after_ts']])
            # comparator = 'ts' in key and 'after_ts' not in key and 'ts_vib' not in key
            if only_ts and 'AWAITING_PARENTS' not in value:
                unfinished_tss.append(key)
        return unfinished_tss

    def get_after_ts_to_restart(self) -> List[str]:
        ''' Go through all_unfinished dictionary and create a list with all
            after_ts ase jobs that started but are not yet finished and
            require to berestarted.

        Returns
        -------
        unfinished_after_tss : list(str)
            a list with unique *py file names for started but yet finished
            minimization of reactants and products after TS calculations

        '''
        all_unfinished = self.get_jobs_to_restart()
        unfinished_after_tss = []
        for key, value in all_unfinished.items():
            if 'after_ts' in key and 'AWAITING_PARENTS' not in value:
                # remove '_f.py' or '_r.py' from the key name
                # unfinished_after_tss.append(key[:-5])
                unfinished_after_tss.append(key)
        # return list(set(unfinished_after_tss))
        return unfinished_after_tss

    def prepare_minima_to_restart(self) -> None:
        ''' If there is at least one optimization step in a .traj file
            for a given minima, this method will create a new *.xyz file for
            that job from the last converged opt step.
            So, when HighLevelRestart is executed, the optimization can resume
            from the last converged step, taking advantade of previous calcs.

        '''
        unfinished_minima = self.get_minima_to_restart()
        print('Restarting minima:')
        for minimum in unfinished_minima:
            print('    {}'.format(minimum))
            metal_symbol, facet, species, prefix, _ = minimum.split('_')
            facetpath = metal_symbol + '_' + facet
            path_to_species = os.path.join(
                facetpath, 'minima', species, prefix)

            try:
                # try to convert last step in *traj file to a new .xyz file
                write(path_to_species + '.xyz',
                      read(path_to_species + '.traj'))
            except (FileNotFoundError, UnknownFileTypeError):
                # continue if *traj file is missing or it is empty
                # hard HighLevelRestart is required
                continue

    def prepare_ts_to_restart(self) -> None:
        ''' If there is at least one optimization step in a .traj file
            for a given ts, this method will create a new *.xyz file for
            that job from the last converged saddle point opt step.
            So, when HighLevelRestart is executed, the optimization can resume
            from the last converged step, taking advantade of previous calcs.

        '''
        unfinished_tss = self.get_tss_to_restart()
        print('Restarting TSs:')
        for ts in unfinished_tss:
            print('    {}'.format(ts))
            prefix, metal_symbol, facet, react, prod, _ = ts.split(
                '_')
            facetpath = metal_symbol + '_' + facet
            rxn_name = react + '_' + prod

            # BUG to be fixed
            ts_name_write = os.path.join(
                prefix + '_' + facetpath + '_' + rxn_name + '_ts')
            ts_name = os.path.join(prefix + '_' + rxn_name)
            path_to_ts = os.path.join(
                facetpath, rxn_name, 'TS_estimate_unique', prefix, ts_name)
            path_to_ts_write = os.path.join(
                facetpath, rxn_name, 'TS_estimate_unique', prefix,
                ts_name_write)

            try:
                # try to convert last step in *traj file to a new .xyz file
                write(path_to_ts_write + '.xyz', read(path_to_ts + '.traj'))
            except (FileNotFoundError, UnknownFileTypeError):
                # continue if *traj file is missing or it is empty
                # hard HighLevelRestart is required
                continue

    def prepare_after_ts_to_restart(self) -> None:
        ''' If there is at least one optimization step in a .traj file
            for a given after_ts, this method will create a new *.xyz file for
            that job from the last converge opt step.
            So, when HighLevelRestart is executed, the optimization can resume
            from the last converged step, taking advantade of previous calcs.

        '''
        unfinished_after_tss = self.get_after_ts_to_restart()
        print('Restarting After TSs:')
        for after_ts in unfinished_after_tss:
            print('    {}'.format(after_ts))
            prefix, metal_atom, facet, react, prod, * \
                _, letter = after_ts.split('_')
            facetpath = metal_atom + '_' + facet
            rxn_name = react + '_' + prod
            after_ts_name = os.path.join(
                prefix + '_' + facetpath + '_' + rxn_name
                + '_after_ts_' + letter[0])
            path_to_after_ts = os.path.join(
                facetpath,
                rxn_name,
                'after_TS',
                prefix,
                after_ts_name)

            try:
                # try to convert last step in *traj file to a new .xyz file
                write(path_to_after_ts + '.xyz',
                      read(path_to_after_ts + '.traj'))
            except (FileNotFoundError, UnknownFileTypeError):
                # continue if *traj file is missing or it is empty
                # hard HighLevelRestart is required
                continue


class HighLevelRestart():
    ''' High level restart means that each ASE job that has an unfinished
        status will be restarted. All balsam jobs will be removed. They will be
        regenerated once ASE jobs start again.

    '''

    def __init__(self):
        # get all python (ASE) jobs
        self.balsamjob = __import__(
            'balsam.launcher.dag', fromlist=['BalsamJob'])
        self.ase_jobs = self.balsamjob.BalsamJob.objects.filter(
            application__contains='python')

    def restart(self) -> None:
        ''' Prepare all unfinished jobs to restart

        '''
        # remove all balsam calculator objects
        self.balsamjob.BalsamJob.objects.filter(name__contains='balsam',
                                                workflow='QE_Socket').delete()

        # update state of every not finished jobs to 'READY'
        for job in self.ase_jobs:
            if job.state not in ['JOB_FINISHED', 'AWAITING_PARENTS']:
                job.state = 'READY'
                job.save()
        self.set_awaiting_status()

    def set_awaiting_status(self):
        for job in self.ase_jobs:
            if len(job.parents) > 2 and job.state != 'JOB_FINISHED':
                job.state = 'AWAITING_PARENTS'
                job.save()
