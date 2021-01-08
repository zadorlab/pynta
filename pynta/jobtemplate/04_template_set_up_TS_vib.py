#!/usr/bin/env python3
import os

from pathlib import Path

from pynta.vib import AfterTS

from balsam.launcher.dag import BalsamJob, add_dependency

socket_calculator = '{socket_calculator}'
facetpath = '{facetpath}'
slab = '{slab}'
repeats = {repeats}
yamlfile = '{yamlfile}'
pytemplate = '{pytemplate}'
pseudo_dir = '{pseudo_dir}'
pseudopotentials = {pseudopotentials}
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
rxn = {rxn}
rxn_name = '{rxn_name}'
cwd = Path.cwd().as_posix()
path_to_ts_vib = os.path.join(creation_dir, facetpath,
                              rxn_name, 'TS_estimate_unique_vib')

after_ts = AfterTS(
    facetpath,
    yamlfile,
    slab,
    repeats,
    creation_dir)

after_ts.set_up_ts_vib(
    socket_calculator,
    rxn,
    pytemplate,
    balsam_exe_settings,
    calc_keywords,
    pseudopotentials,
    pseudo_dir)

workflow_name = facetpath + '_04_' + rxn_name
dependency_workflow_name = facetpath + '_03_' + rxn_name
dependent_workflow_name = facetpath + '_05_' + rxn_name

pending_simulations = BalsamJob.objects.filter(
    workflow__contains=dependency_workflow_name
).exclude(state="JOB_FINISHED")

pending_simulations_dep = BalsamJob.objects.filter(
    workflow__contains=dependent_workflow_name
).exclude(state="JOB_FINISHED")


for py_script in Path(path_to_ts_vib).glob('*.py'):
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
    # all job_to_add_ are childs of 03 job for a given reaction
    for job in pending_simulations:
        add_dependency(job, job_to_add)  # parent, child
    # do not run 05 until all 04 for a given reaction are done
    for job in pending_simulations_dep:
        add_dependency(job_to_add, job)  # parent, child
