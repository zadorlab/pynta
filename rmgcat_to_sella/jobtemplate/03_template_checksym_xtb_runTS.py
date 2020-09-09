#!/usr/bin/env python3

import os

from pathlib import Path

from rmgcat_to_sella.ts import TS
from rmgcat_to_sella.io import IO

from balsam.launcher.dag import BalsamJob, add_dependency

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
cwd = Path.cwd().as_posix()
path_to_ts_estimate_uq = os.path.join(
    facetpath, rxn_name, 'TS_estimate_unique')

ts = TS(
    facetpath,
    slab,
    ts_estimate_dir,
    yamlfile,
    repeats,
    creation_dir)

ts.create_unique_ts_all(
    rxn,
    pytemplate,
    pseudopotentials,
    pseudo_dir,
    balsam_exe_settings,
    calc_keywords)

workflow_name = yamlfile + facetpath + '03' + rxn_name
dependency_workflow_name = yamlfile + facetpath + '02' + rxn_name
dependent_workflow_name = yamlfile + facetpath + '04' + rxn_name

pending_simulations = BalsamJob.objects.filter(
    workflow__contains=dependency_workflow_name
).exclude(state="JOB_FINISHED")

pending_simulations_dep = BalsamJob.objects.filter(
    workflow__contains=dependent_workflow_name
).exclude(state="JOB_FINISHED")


for py_script in Path(path_to_ts_estimate_uq).glob('*.py'):
    job_dir, script_name = os.path.split(str(py_script))
    job_to_add = BalsamJob(
        name=script_name,
        workflow=workflow_name,
        application='python',
        args=cwd + '/' + str(py_script),
        input_files='',
        user_workdir=job_dir,
        node_packing_count=48,
        ranks_per_node=1,
    )
    job_to_add.save()
    # all job_to_add_ are childs of 02 job for a given reaction
    for job in pending_simulations:
        add_dependency(job, job_to_add)  # parent, child
    # do not run 04 until all 03 for a given reaction are done
    for job in pending_simulations_dep:
        add_dependency(job_to_add, job)  # parent, child
