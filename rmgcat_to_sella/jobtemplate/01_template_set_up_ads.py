#!/usr/bin/env python3
import os
from glob import glob
from pathlib import Path

from rmgcat_to_sella.adsorbates import Adsorbates
from rmgcat_to_sella.ts import TS

# from balsam.launcher.dag import BalsamJob, add_dependency


facetpath = 'Cu_111'
slab = 'Cu_111_slab_opt.xyz'
repeats = (3, 3, 1)
yamlfile = 'reactions.yaml'
pytemplate = '/Users/mgierad/00_SANDIA_WORK/05_rmgcat_to_stella/rmgcat_to_sella/pytemplate/pytemplate_relax_Cu_111_ads.py'
pseudopotentials = dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF', H='H.pbe-kjpaw_psl.1.0.0.UPF',
                        O='O.pbe-n-kjpaw_psl.1.0.0.UPF', C='C.pbe-n-kjpaw_psl.1.0.0.UPF')
pseudo_dir = '/home/mgierad/espresso/pseudo'
balsam_exe_settings = {'num_nodes': 1,
                       'ranks_per_node': 48, 'threads_per_rank': 1}
calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing', 'smearing': 'marzari-vanderbilt',
                 'degauss': 0.01, 'ecutwfc': 40, 'nosym': True, 'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}
creation_dir = '/Users/mgierad/00_SANDIA_WORK/05_rmgcat_to_stella/test/rmgcat_to_sella/00_code_test/code_test_yaml_parser/'
ts_dir = 'TS_estimate'

# put_adsorbates = Adsorbates(
#     facetpath, slab, repeats, yamlfile, creation_dir)
# put_adsorbates.adjacency_to_3d()
# put_adsorbates.create_relax_jobs(
#     pytemplate, pseudopotentials, pseudo_dir,
#     balsam_exe_settings, calc_keywords
# )

ts = TS(facetpath, slab, ts_dir, yamlfile, repeats, creation_dir)
dependancy_dict = ts.depends_on()


def jobs_to_be_finished(dependancy_dict, rxn_name):
    print('---')
    print('Reaction : {}'.format(rxn_name))
    print('Jobs to be finished: ')
    jobs_to_be_finished = dependancy_dict[rxn_name]
    print(jobs_to_be_finished)
    return jobs_to_be_finished


def submit_all_jobs(dependancy_dict):
    all_jobs = []
    for rxn_name in dependancy_dict.keys():
        for job in jobs_to_be_finished(dependancy_dict, rxn_name):
            if job not in all_jobs:
                all_jobs.append(job)

    print('===')
    print('All jobs to be submitted:')
    print(all_jobs)
    return all_jobs


submit_all_jobs(dependancy_dict)


def run_01(dependancy_dict, rxn_name):
    all_jobs = submit_all_jobs(dependancy_dict)
    cwd = Path.cwd().as_posix()
    workflow_name = yamlfile + facetpath + '01'
    # specify dependant 02 for given reactions
    dependent_workflow_name = yamlfile+facetpath+'02'+rxn_name

    # have to find a way to specify dependancy which species have to be
    # calculated for a given reaction
    pending_simulations_dep = BalsamJob.objects.filter(
        workflow__contains=dependent_workflow_name
    ).exclude(state="JOB_FINISHED")

    # pending_simulations_dep = BalsamJob.objects.filter(
    #     workflow__contains=dependent_workflow_name).exclude(state="JOB_FINISHED")

    # lookup_phrase = os.path.join(facetpath, 'minima', rxn_name, '*py')
    # for py_script in glob(lookup_phrase):

    # make a list with all jobs that have to be finished for a given reaction
    # make sure be submit jobs only once
    for py_script in jobs_to_be_finished(dependancy_dict, rxn_name):
        job_dir = cwd + '/' + '/'.join(py_script.strip().split('/')[:-1])
        script_name = py_script.strip().split('/')[-1]
        job_to_add = BalsamJob(
            name=script_name,
            workflow=workflow_name,
            application='python',
            args=cwd + '/' + py_script,
            input_files='',
            user_workdir=job_dir,
            node_packing_count=64,
            ranks_per_node=1,
        )
        job_to_add.save()
        for job in pending_simulations_dep:
            add_dependency(job_to_add, job)  # parent, child

        # job_to_add should be for O, H, OH for 02 job OH_O+H


'''
Or better, a pseudocode here

def check_for_status(args):
    ts = TS(facetpath, slab, ts_dir, yamlfile, repeats, creation_dir)
    dependancy_dict = ts.depends_on()

    finished_jobs = []

    for jobs in dependancy_dict.values()
        for job in jobs:
            if jobs.status('FINISHED'):
                finished_jobs.append(job)
        if all(finished_jobs):
            run_01(...)

def run_01(...):
    arguments = ...

    put_adsorbates = Adsorbates(
        facetpath, slab, repeats, yamlfile, creation_dir)
    put_adsorbates.adjacency_to_3d()
    put_adsorbates.create_relax_jobs(
        pytemplate, pseudopotentials, pseudo_dir,
        balsam_exe_settings, calc_keywords
    )

    cwd = Path.cwd().as_posix()
    workflow_name = yamlfile + facetpath + '01'

    dependent_workflow_name = yamlfile+facetpath+'02'
    pending_simulations_dep = BalsamJob.objects.filter(
        workflow__contains=dependent_workflow_name
    ).exclude(state="JOB_FINISHED")

    lookup_phrase = os.path.join(facetpath, 'minima', rxn_name, '*py')
    # for py_script in glob('Cu_111/minima/*.py'):
    for py_script in glob(lookup_phrase):
        job_dir = cwd + '/' + '/'.join(py_script.strip().split('/')[:-1])
        script_name = py_script.strip().split('/')[-1]
        job_to_add = BalsamJob(
            name=script_name,
            workflow=workflow_name,
            application='python',
            args=cwd + '/' + py_script,
            input_files='',
            user_workdir=job_dir,
            node_packing_count=64,
            ranks_per_node=1,
        )
        job_to_add.save()
        for job in pending_simulations_dep:
            add_dependency(job_to_add, job)  # parent, child



'''
