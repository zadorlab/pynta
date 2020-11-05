#!/usr/bin/env python3
import os
from rmgcat_to_sella.vib import minimaVib
from rmgcat_to_sella.io import IO
from balsam.launcher.dag import BalsamJob, add_dependency

facetpath = '{facetpath}'
yamlfile = '{yamlfile}'
repeats = {repeats}
pytemplate = '{pytemplate}'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
node_packing_count = {node_packing_count}
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
yamlfile_path = os.path.join(creation_dir, yamlfile)

mv = minimaVib(facetpath, creation_dir)
mv.create_minima_vib_all(facetpath, yamlfile_path, pytemplate,
                         balsam_exe_settings, pseudo_dir,
                         pseudopotentials, calc_keywords,
                         creation_dir)

# job_to_add = BalsamJob(
#     name=submit_py,
#     workflow=workflow_name,
#     application='python',
#     args=str(submit_py),
#     input_files='',
#     user_workdir=submit_py_path,
#     node_packing_count=48,
#     ranks_per_node=1,
# )
# job_to_add.save()

# pending_simulations = BalsamJob.objects.filter(
#     name=py_script_prev_opt).exclude(state="JOB_FINISHED"))

# # add dependencies
# for job in pending_simulations:
#     for sub_job in job:
#         add_dependency(sub_job, job_to_add)  # parent, child
