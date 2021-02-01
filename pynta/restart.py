import os
from pathlib import Path
import re
import shutil
import sys
from ase.io import read, write
from ase.io.formats import UnknownFileTypeError
from typing import Dict, List


class LowLevelRestart():
    ''' Low level restart means that each quantum chemistry job that can be
    restarted will be restarted.

    These include:

        * adsorbates minimization
        * TSs optimization
        * optimization of reactants and products as from TS (IRC alternative)

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
        ''' Run low level restart - restart jobs from the last converged point,
        if possible.

        '''
        self.prepare_minima_to_restart()
        self.remove_partialy_generated_xtb_run()
        self.prepare_ts_to_restart()
        self.prepare_after_ts_to_restart()

    def get_jobs_to_restart(self) -> Dict[str, str]:
        ''' Get a dictionary with all ase_jobs that did not finish

        Returns
        -------
        all_unfinished: Dict[str, str]
            a dictionary with all unfinished jobs,
            e.g.

            >>> all_unfinished = {Cu_111_C_01_relax.py: 'RUNNING'}

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
            a list with all :literal:`*.py` files for minima optimization which
            did not finished

        '''
        all_unfinished = self.get_jobs_to_restart()
        unfinished_minima = []
        for key in all_unfinished.keys():
            if 'relax' in key:
                unfinished_minima.append(key)
        return unfinished_minima

    def get_ts_xtb_to_remove(self) -> List[str]:
        ''' Get a list of all ts_xtb jobs for which initial :literal:`*.xyz`
        files were generated but :literal:`*.py` not yet.

        This method will take care of edge cases where the workflow execution
        stopped while symmetry of initially generated :literal:`*.xyz` were
        checking

        Returns
        -------
        unfinished_ts_xtb : List[str]
            a list of all ts_xtb jobs for which :literal:`*.xyz` files have to
            be removed

        '''
        all_unfinished = self.get_jobs_to_restart()
        unfinished_ts_xtb = []
        for key, value in all_unfinished.items():
            if 'xtb' in key and 'AWAITING_PARENTS' not in value:
                unfinished_ts_xtb.append(key)
        return unfinished_ts_xtb

    def get_tss_to_restart(self) -> List[str]:
        ''' Go through all_unfinished dictionary and create a list with all
        TS ase jobs that started but are not yet finished and require to be
        restarted.

        Returns
        -------
        unfinished_tss : list(str)
            a list with all :literal:`*.py` files for started but yet
            unfinished 1st order saddle point optimizations

        '''
        all_unfinished = self.get_jobs_to_restart()
        unfinished_tss = []
        for key, value in all_unfinished.items():
            only_ts = 'ts' in key and all(
                [k not in key for k in ['ts_vib', 'after_ts']])
            # comparator = 'ts' in key and 'after_ts' not in key
            # and 'ts_vib' not in key
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
            a list with unique :literal:`*.py` file names for started but
            yet finishedminimization of reactants and products
            after TS calculations

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
        for a given minima, this method will create a new :literal:`*.xyz`
        file for that job from the last converged opt step.

        So, when :class:`pynta.restart.HighLevelRestart` is executed,
        the optimization can resume from the last converged step, taking
        advantade of previous calculations.

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
        if not unfinished_minima:
            print('    Nothing to restart!')

    def remove_partialy_generated_xtb_run(self):
        ''' Remove all partially generated :literal:`*.xyz` files as described
        in :meth:`pynta.restart.LowLevelRestart.get_ts_xtb_to_remove`

        '''
        unfinished_ts_xtb = self.get_ts_xtb_to_remove()

        print('Removing partially generated TS_xtb jobs:')
        for ts_xtb in unfinished_ts_xtb:
            rxn_name = re.search('xtb_(.*).py', ts_xtb).group(1)
            facetpath = re.search('02_(.*)_set', ts_xtb).group(1)
            path_to_ts_xtb = os.path.join(
                self.current_dir, facetpath, rxn_name)

            print('    * trying {}'.format(ts_xtb))
            try:
                shutil.rmtree(path_to_ts_xtb)
                print('        Done')
            except FileNotFoundError:
                print('        \'.xyz\' files not found. Nothing to remove')
                print('        ## check {} ##'.format(path_to_ts_xtb))
        if not unfinished_ts_xtb:
            print('    Nothing to restart!')

    def prepare_ts_to_restart(self) -> None:
        ''' If there is at least one optimization step in a .traj file
        for a given ts, this method will create a new :literal:`*.xyz` file for
        that job from the last converged saddle point opt step.

        So, when :class:`pynta.restart.HighLevelRestart` is executed,
        the optimization can resume from the last converged step, taking
        advantade of previous calculations.

        '''
        unfinished_tss = self.get_tss_to_restart()

        print('Restarting TSs:')

        for ts in unfinished_tss:
            print('    {}'.format(ts))
            prefix, metal_symbol, facet, react, prod, _ = ts.split(
                '_')
            facetpath = metal_symbol + '_' + facet
            rxn_name = react + '_' + prod

            ts_name = os.path.join(
                prefix + '_' + facetpath + '_' + rxn_name)
            path_to_ts = os.path.join(
                facetpath, rxn_name, 'TS_estimate_unique', prefix, ts_name)

            try:
                # try to convert last step in *traj file to a new .xyz file
                write(path_to_ts + '_ts.xyz',
                      read(path_to_ts + '.traj'))
            except (FileNotFoundError, UnknownFileTypeError):
                # continue if *traj file is missing or it is empty
                # hard :meth:HighLevelRestart is required
                continue
        if not unfinished_tss:
            print('    Nothing to restart!')

    def prepare_after_ts_to_restart(self) -> None:
        ''' If there is at least one optimization step in a .traj file
        for a given after_ts, this method will create a new :literal:`*.xyz`
        file for that job from the last converge opt step.

        So, when :class:`pynta.restart.HighLevelRestart` is executed,
        the optimization can resume from the last converged step, taking
        advantade of previous calculations.

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
        if not unfinished_after_tss:
            print('    Nothing to restart!')


class HighLevelRestart():
    ''' High level restart means that each ASE job that has an unfinished
    status will be restarted. All balsam jobs will be removed. They will be
    regenerated once ASE jobs start again.

    '''
    error_files = ['*/pwscf*', '*/core']

    def __init__(self):
        self.current_dir = os.getcwd()
        # get all python (ASE) jobs
        self.balsamjob = __import__(
            'balsam.launcher.dag', fromlist=['BalsamJob'])
        self.ase_jobs = self.balsamjob.BalsamJob.objects.filter(
            application__contains='python')

    def restart(self) -> None:
        ''' Prepare all unfinished jobs to restart.

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
        self.remove_error_files()
        self.remove_empty_pickle_files()

    def remove_empty_pickle_files(self) -> None:
        ''' Remove all empty pckl files from unfinished vib calculations

        '''
        print('Removing unfinished/empty *.pckl files:')
        all_pickle_files = Path(self.current_dir).glob('**/*pckl')
        removed_pckls = []
        for pckl in all_pickle_files:
            if os.stat(pckl).st_size == 0:
                rxn_name = os.path.dirname(pckl).split('/')[-3]
                f_name = os.path.basename(pckl)
                print('    rxn: {} file: {}'.format(rxn_name, f_name))
                os.remove(pckl)
                removed_pckls.append(pckl)
        if not removed_pckls:
            print('    Nothing to remove!')

    def set_awaiting_status(self) -> None:
        ''' Make sure that jobs which not yet started and depends on other jobs
        have an ``'AWAITING_PARENTS'`` status

        '''
        for job in self.ase_jobs:
            if len(job.parents) > 2 and job.state != 'JOB_FINISHED':
                job.state = 'AWAITING_PARENTS'
                job.save()

    @staticmethod
    def get_workflow_path() -> str:
        ''' Read run_me.sh` submission script and extract path where stdout
        balsam files are, e.g.

        >>> workflow_path = '/home/user_name/myWorkflow/data/QE_Socket'

        Returns
        -------
        workflow_path : str
            an absolute path to stdout balsam files

        '''
        try:
            with open('run_me.sh') as infile:
                lines = infile.readlines()
                for line in lines:
                    if line.startswith('source balsamactivate'):
                        workflow_path = line.split()[2]
                        if '~' in workflow_path:
                            workflow_path = os.path.expanduser(workflow_path)
                        return workflow_path
        except FileNotFoundError:
            print('run_me.sh not found in \n'
                  '    {}'.format(os.getcwd()))
            sys.exit()

    def remove_error_files(self):
        ''' Remove all error/unfinished files from previous calculations
        including :literal:`pwscf*`, `core`

        '''
        workflow_path = HighLevelRestart.get_workflow_path()
        path_to_balsam_out = os.path.join(workflow_path, 'data', 'QE_Socket')

        print('Removing unfinished temp files:')
        for err_files_type in self.error_files:
            err_files = Path(path_to_balsam_out).glob(err_files_type)

            for err in err_files:
                try:
                    print('Removing file: {}'.format(err))
                    os.remove(err)
                except IsADirectoryError:
                    print('Removing directory: {}'.format(err))
                    shutil.rmtree(err)
