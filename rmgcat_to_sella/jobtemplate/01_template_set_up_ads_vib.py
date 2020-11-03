#!/usr/bin/env python3
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

mv = minimaVib(facetpath, creation_dir)
mv.create_minima_vib_all(adsorbate, pytemplate, balsam_exe_settings,
                         pseudo_dir, pseudopotentials, calc_keywords,
                         creation_dir)

py_scripts = mv.dependency_minima_vib(adsorbate)
print(py_scripts)

# dependancy_dict = IO().depends_on(facetpath, yamlfile, creation_dir)

# # keep track of all submitted jobs (all unique)
# all_submitted_jobs = []

# # specify dependant 02 for given reactions
# for rxn_name in dependancy_dict.keys():
#     workflow_name = facetpath + '_01_' + rxn_name
#     dependent_workflow_name = facetpath + '_02_' + rxn_name

#     # get all dependant workflow BalsamJobs objects (should be one)
#     pending_simulations_dep = BalsamJob.objects.filter(
#         workflow__contains=dependent_workflow_name
#     ).exclude(state="JOB_FINISHED")

#     # for each reaction keep track of its dependencies
#     # e.g. for OH --> O + H those have to be finished OH, O and H
#     # Make sure each species is calculated only once
#     # e.g. H in CH --> C + H and OH --> O + H
#     new_unique_submission = []
#     dependancy = []

#     jobs_to_be_finished = dependancy_dict[rxn_name]

#     for py_script in jobs_to_be_finished:
#         py_script_dir = os.path.join(
#             creation_dir, facetpath, 'minima', py_script)
#         job_dir, _ = os.path.split(py_script_dir)

#         # get all unique submission
#         if py_script not in all_submitted_jobs:
#             new_unique_submission.append(py_script)
#             all_submitted_jobs.append(py_script)

#             # create a BalsamJob object, i.e. submit all unique jobs
#             job_to_add = BalsamJob(
#                 name=py_script,
#                 workflow=workflow_name,
#                 application='python',
#                 args=py_script_dir,
#                 input_files='',
#                 user_workdir=job_dir,
#                 node_packing_count={node_packing_count},
#                 ranks_per_node=1,
#             )
#             job_to_add.save()

#         # for a given rxn_name, get all BalsamJob objects that it depends on
#         dependancy.append(BalsamJob.objects.filter(
#             name=py_script).exclude(state="JOB_FINISHED"))

#     # add dependencies
#     for job in pending_simulations_dep:
#         for adding_job in dependancy:
#             add_dependency(adding_job, job)  # parent, child
