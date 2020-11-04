#!/usr/bin/env python3
import os
from rmgcat_to_sella.vib import minimaVib
from balsam.launcher.dag import BalsamJob, add_dependency

facetpath = '{facetpath}'
repeats = {repeats}
adsorbate = '{adsorbate}'
prefix = '{prefix}'
py_script_prev_opt = '{py_script_prev_opt}'
pytemplate = '{pytemplate}'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
node_packing_count = {node_packing_count}
submit_py = '{{}}_{{}}_{{}}_vib.py'.format(facetpath, prefix, adsorbate)
submit_py_path = os.path.join(creation_dir, facetpath, 'minima_vib', submit_py)

mv = minimaVib(facetpath, creation_dir)
mv.create_minima_vib_all(adsorbate, prefix, pytemplate, balsam_exe_settings,
                         pseudo_dir, pseudopotentials, calc_keywords,
                         creation_dir)

workflow_name = '{{}}_01_{{}}_{{}}_vib'.format(facetpath, prefix, adsorbate)

job_to_add = BalsamJob(
    name=submit_py,
    workflow=workflow_name,
    application='python',
    args=str(submit_py),
    input_files='',
    user_workdir=submit_py_path,
    node_packing_count={node_packing_count},
    ranks_per_node=1,
)
job_to_add.save()

# pending_simulation = BalsamJob.objects.filter(
#     name=py_script_prev_opt).exclude(state="JOB_FINISHED"))

# # add dependencies
# for job in pending_simulations:
#     for sub_job in job:
#         add_dependency(sub_job, job_to_add)  # parent, child
