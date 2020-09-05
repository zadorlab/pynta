#!/usr/bin/env python3

import os
from pathlib import Path

from rmgcat_to_sella.ts import TS

from balsam.launcher.dag import BalsamJob, add_dependency

slab = '{slab}'
repeats = {repeats}
yamlfile = '{yamlfile}'
facetpath = '{facetpath}'
rotAngle = {rotAngle}
scfactor = {scfactor}
scfactor_surface = {scfactor_surface}
pytemplate_xtb = '{pytemplate_xtb}'
species_list = {species_list}
current_dir = os.path.dirname(os.getcwd())
minima_dir = os.path.join(facetpath, 'minima')
scaled1 = {scaled1}
scaled2 = {scaled2}
ts_dir = 'TS_estimate'
workflow_name = yamlfile+facetpath+'02'
dependency_workflow_name = yamlfile+facetpath+'01'
creation_dir = '{creation_dir}'
rxn = {rxn}
rxn_name = '{rxn_name}'

ts = TS(
    facetpath,
    slab, ts_dir,
    yamlfile,
    repeats,
    creation_dir)
# ts.copy_minimas_prev_calculated(current_dir, species, minima_dir)
ts.prepare_ts_estimate(
    rxn,
    scfactor,
    scfactor_surface,
    rotAngle,
    pytemplate_xtb,
    species_list,
    scaled1,
    scaled2)

dependent_workflow_name = yamlfile+facetpath+'03'
pending_simulations_dep = BalsamJob.objects.filter(
    workflow__contains=dependent_workflow_name
).exclude(state="JOB_FINISHED")

pending_simulations = BalsamJob.objects.filter(
    workflow__contains=dependency_workflow_name
).exclude(state="JOB_FINISHED")
cwd = Path.cwd().as_posix()
path_to_ts_estimate = os.path.join(facetpath, rxn_name, 'TS_estimate')
for py_script in Path(path_to_ts_estimate).glob('**/*.py'):
    print(py_script)
    job_dir, script_name = os.path.split(str(py_script))
    job_to_add = BalsamJob(
        name=script_name,
        workflow=workflow_name,
        application='python',
        args=cwd + '/' + str(py_script),
        input_files='',
        ranks_per_node=1,
        node_packing_count=48,
        user_workdir=job_dir,
    )
    job_to_add.save()
    for job in pending_simulations:
        add_dependency(job, job_to_add)  # parent, child
    for job in pending_simulations_dep:
        add_dependency(job_to_add, job)  # parent, child
