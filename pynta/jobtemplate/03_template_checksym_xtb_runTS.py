#!/usr/bin/env python3

import os

from pathlib import Path

from pynta.ts import TS

from balsam.launcher.dag import BalsamJob, add_dependency

socket_calculator = '{socket_calculator}'
slab = '{slab}'
repeats = {repeats}
yamlfile = '{yamlfile}'
facetpath = '{facetpath}'
pytemplate = '{pytemplate}'
ts_estimate_dir = 'TS_estimate'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
rxn = {rxn}
rxn_name = '{rxn_name}'

ts_estimate_path = os.path.join(
    creation_dir, facetpath, rxn_name, ts_estimate_dir)
ts_estimate_path_uq = os.path.join(
    creation_dir, facetpath, rxn_name, ts_estimate_dir + '_unique')

ts = TS(
    facetpath,
    slab,
    ts_estimate_dir,
    yamlfile,
    repeats,
    creation_dir)

ts.create_unique_ts_all(
    socket_calculator,
    ts_estimate_path,
    rxn_name,
    pytemplate,
    pseudopotentials,
    pseudo_dir,
    balsam_exe_settings,
    calc_keywords)

workflow_name = facetpath + '_03_' + rxn_name
dependency_workflow_name = facetpath + '_02_' + rxn_name
dependent_workflow_name = facetpath + '_04_' + rxn_name

pending_simulations = BalsamJob.objects.filter(
    workflow__contains=dependency_workflow_name
).exclude(state="JOB_FINISHED")

pending_simulations_dep = BalsamJob.objects.filter(
    workflow__contains=dependent_workflow_name
).exclude(state="JOB_FINISHED")


for py_script in Path(ts_estimate_path_uq).glob('*.py'):
    job_dir, script_name = os.path.split(str(py_script))
    job_to_add = BalsamJob(
        name=script_name,
        workflow=workflow_name,
        application='python',
        args=str(py_script),
        input_files='',
        user_workdir=job_dir,
        node_packing_count={node_packing_count},
        ranks_per_node=1,
    )
    job_to_add.save()
    # all job_to_add_ are childs of 02 job for a given reaction
    for job in pending_simulations:
        add_dependency(job, job_to_add)  # parent, child
    # do not run 04 until all 03 for a given reaction are done
    for job in pending_simulations_dep:
        add_dependency(job_to_add, job)  # parent, child
