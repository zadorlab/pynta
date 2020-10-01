#!/usr/bin/env python3
import os

from pathlib import Path

from rmgcat_to_sella.vib import AfterTS

from balsam.launcher.dag import BalsamJob, add_dependency

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
cwd = Path.cwd().as_posix()
path_to_after_TS = os.path.join(
    facetpath, rxn_name, 'after_TS')

after_ts = AfterTS(facetpath, yamlfile, slab, repeats)
after_ts.prepare_all(rxn)

workflow_name = facetpath + '_04_' + rxn_name
dependency_workflow_name = facetpath + '_03_' + rxn_name

pending_simulations = BalsamJob.objects.filter(
    workflow__contains=dependency_workflow_name
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
