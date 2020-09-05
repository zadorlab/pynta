#!/usr/bin/env python3

from glob import glob
from pathlib import Path

from rmgcat_to_sella.adsorbates import Adsorbates
from rmgcat_to_sella.io import IO

from balsam.launcher.dag import BalsamJob, add_dependency

facetpath = '{facetpath}'
slab = '{slabopt}'
repeats = {repeats}
yamlfile = '{yamlfile}'
pytemplate = '{pytemplate}'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'

put_adsorbates = Adsorbates(facetpath, slab, repeats, yamlfile, creation_dir)
put_adsorbates.adjacency_to_3d()
put_adsorbates.create_relax_jobs(
    pytemplate, pseudopotentials, pseudo_dir,
    balsam_exe_settings, calc_keywords
)

io = IO()
dependancy_dict = io.depends_on(facetpath, yamlfile)

cwd = Path.cwd().as_posix()
workflow_name = yamlfile + facetpath + '01'


def jobs_to_be_finished(dependancy_dict, rxn_name):
    jobs_to_be_finished = dependancy_dict[rxn_name]
    return jobs_to_be_finished


def run_01(dependancy_dict):
    # all_jobs = all_unique_01_jobs(dependancy_dict)
    cwd = Path.cwd().as_posix()
    workflow_name = yamlfile + facetpath + '01'

    # keep track of all submitted jobs (all unique)
    all_submitted_jobs = []

    # specify dependant 02 for given reactions
    for rxn_name in dependancy_dict.keys():
        dependent_workflow_name = yamlfile+facetpath+'02'+rxn_name

        # have to find a way to specify dependancy which species have to be
        # calculated for a given reaction
        pending_simulations_dep = BalsamJob.objects.filter(
            workflow__contains=dependent_workflow_name
        ).exclude(state="JOB_FINISHED")

        # for each reaction keep track of its dependencies
        # e.g. for OH --> O + H those have to be finished OH, O and H
        # Make sure each species is calculated only once
        # e.g. H in CH --> C + H and OH --> O + H
        new_unique_submission = []
        for py_script in jobs_to_be_finished(dependancy_dict, rxn_name):
            job_dir = cwd + '/' + '/'.join(py_script.strip().split('/')[:-1])
            # get all unique submission
            if py_script not in all_submitted_jobs:
                new_unique_submission.append(py_script)
                all_submitted_jobs.append(py_script)
            if py_script in new_unique_submission:
                job_to_add = BalsamJob(
                    name=py_script,
                    workflow=workflow_name,
                    application='python',
                    args=cwd + '/' + py_script,
                    input_files='',
                    user_workdir=job_dir,
                    node_packing_count=48,
                    ranks_per_node=1,
                )
                job_to_add.save()
                for job in pending_simulations_dep:
                    add_dependency(job_to_add, job)  # parent, child
            else:
                job_to_add = BalsamJob(
                    name=py_script,
                    workflow=workflow_name,
                    application='python',
                    args=cwd + '/' + py_script,
                    input_files='',
                    user_workdir=job_dir,
                    node_packing_count=48,
                    ranks_per_node=1,
                )
                # job_to_add.save()
                for job in pending_simulations_dep:
                    add_dependency(job_to_add, job)  # parent, child
