#!/usr/bin/env python3
import os

from pathlib import Path

from pynta.vib import AfterTS

from balsam.launcher.dag import BalsamJob, add_dependency

socket_calculator = '{socket_calculator}'
slab = '{slab}'
repeats = {repeats}
yamlfile = '{yamlfile}'
facetpath = '{facetpath}'
pytemplate = '{pytemplate}'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
rxn = {rxn}
rxn_name = '{rxn_name}'
path_to_after_ts = os.path.join(creation_dir, facetpath,
                                rxn_name, 'after_TS')

after_ts = AfterTS(
    facetpath,
    yamlfile,
    slab,
    repeats,
    creation_dir)

after_ts.prepare_opt_after_ts(
    socket_calculator,
    rxn,
    pytemplate,
    balsam_exe_settings,
    calc_keywords,
    pseudopotentials,
    pseudo_dir)

workflow_name = facetpath + '_05_' + rxn_name
dependency_workflow_name = facetpath + '_04_' + rxn_name

pending_simulations = BalsamJob.objects.filter(
    workflow__contains=dependency_workflow_name
).exclude(state="JOB_FINISHED")

for py_script in Path(path_to_after_ts).glob('*.py'):
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
