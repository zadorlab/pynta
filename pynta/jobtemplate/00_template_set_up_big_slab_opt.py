#!/usr/bin/env python3
from pynta.main import WorkFlow
import os

from balsam.launcher.dag import BalsamJob, add_dependency

socket_calculator = '{socket_calculator}'
facetpath = '{facetpath}'
slab_name = '{slab_name}'
repeats = {repeats}
pytemplate = '{pytemplate}'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'

WorkFlow().create_big_slab_pyjob(
    socket_calculator,
    pytemplate,
    facetpath,
    slab_name,
    repeats,
    creation_dir)

workflow_name = facetpath + '_big_slab_opt'
dependency_workflow_name = facetpath + '_00_'

py_script_fname = os.path.join(facetpath + '_big_slab_opt_job.py')

pending_simulations = BalsamJob.objects.filter(
    workflow__contains=dependency_workflow_name
).exclude(state="JOB_FINISHED")

# create a BalsamJob object, i.e. submit all unique jobs
job_to_add = BalsamJob(
    name=py_script_fname,
    workflow=workflow_name,
    application='python',
    args=str(py_script_fname),
    input_files='',
    user_workdir='{job_file_dir_name}',
    node_packing_count={node_packing_count},
    ranks_per_node=1,
)
job_to_add.save()

# add dependencies
for job in pending_simulations:
    add_dependency(job, job_to_add)  # parent, child
