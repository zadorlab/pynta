#!/usr/bin/env python3
from operator import sub
from pathlib import Path
import os

from rmgcat_to_sella.vib import minimaVib
from rmgcat_to_sella.io import IO

# from balsam.launcher.dag import BalsamJob, add_dependency

facetpath = '{facetpath}'
repeats = {repeats}
adsorbate = '{adsorbate}'
pytemplate = '{pytemplate}'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
node_packing_count = {node_packing_count}
submit_py = os.path.join(facetpath + '_' + adsorbate + '_vib.py')
submit_py_path = os.path.join(creation_dir, facetpath, 'minima_vib', submit_py)

mv = minimaVib(facetpath, creation_dir)
mv.create_minima_vib_all(adsorbate, pytemplate, balsam_exe_settings,
                         pseudo_dir, pseudopotentials, calc_keywords,
                         creation_dir)

workflow_name = facetpath + '_01_' + adsorbate + '_vib'

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

py_scripts = mv.dependency_minima_vib(adsorbate)
pending_simulations = []
for py_script in py_scripts:
    pending_simulations.append(BalsamJob.objects.filter(
        name=py_script).exclude(state="JOB_FINISHED"))

# add dependencies
for job in pending_simulations:
    for sub_job in job:
        add_dependency(sub_job, job_to_add)  # parent, child
