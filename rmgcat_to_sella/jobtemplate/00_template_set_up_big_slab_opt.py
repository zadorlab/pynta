#!/usr/bin/env python3
from pathlib import Path
import os


from balsam.launcher.dag import BalsamJob, add_dependency

facetpath = '{facetpath}'
slab_name = '{slab_name}'
repeats = {repeats}
pytemplate = '{pytemplate}'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
py_script = '{py_script}'


workflow_name = facetpath + '_00_big_slab_opt'
dependency_workflow_name = facetpath + '_00_'

# create a BalsamJob object, i.e. submit all unique jobs
job_to_add = BalsamJob(
    name=py_script,
    workflow=workflow_name,
    application='python',
    args=py_script_dir,
    input_files='',
    user_workdir=job_dir,
    node_packing_count={node_packing_count},
    ranks_per_node=1,
)
job_to_add.save()

# for a given rxn_name, get all BalsamJob objects that it depends on
dependancy.append(BalsamJob.objects.filter(
    name=py_script).exclude(state="JOB_FINISHED"))

# add dependencies
for job in pending_simulations_dep:
    for adding_job in dependancy:
        add_dependency(adding_job, job)  # parent, child
