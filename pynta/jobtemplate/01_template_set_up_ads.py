#!/usr/bin/env python3
import os

from pynta.adsorbates import Adsorbates
from pynta.io import IO

from balsam.launcher.dag import BalsamJob, add_dependency

socket_calculator = '{socket_calculator}'
facetpath = '{facetpath}'
slab = '{slab}'
repeats = {repeats}
yamlfile = '{yamlfile}'
pytemplate = '{pytemplate}'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'

put_adsorbates = Adsorbates(facetpath, slab, repeats, yamlfile, creation_dir)
put_adsorbates.adjacency_to_3d()
put_adsorbates.create_relax_jobs(
    socket_calculator, pytemplate, pseudopotentials, pseudo_dir,
    balsam_exe_settings, calc_keywords
)

dependancy_dict = IO().depends_on(facetpath, yamlfile, creation_dir)

# keep track of all submitted jobs (all unique)
all_submitted_jobs = []

# specify dependant 02 for given reactions
for rxn_name in dependancy_dict.keys():
    workflow_name = facetpath + '_01_' + rxn_name
    # ts estimate jobs
    dependent_workflow_name_1 = facetpath + '_02_' + rxn_name
    # get all dependant workflow BalsamJobs objects (should be one)
    pending_simulations_dep_1 = BalsamJob.objects.filter(
        workflow__contains=dependent_workflow_name_1
    ).exclude(state="JOB_FINISHED")

    # for each reaction keep track of its dependencies
    # e.g. for OH --> O + H those have to be finished OH, O and H
    # Make sure each species is calculated only once
    # e.g. H in CH --> C + H and OH --> O + H
    new_unique_submission = []
    dependancy = []

    jobs_to_be_finished = dependancy_dict[rxn_name]

    for py_script in jobs_to_be_finished:
        py_script_dir = os.path.join(
            creation_dir, facetpath, 'minima', py_script)
        job_dir, _ = os.path.split(py_script_dir)

        # get all unique submission
        if py_script not in all_submitted_jobs:
            new_unique_submission.append(py_script)
            all_submitted_jobs.append(py_script)

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
    for job in pending_simulations_dep_1:
        for adding_job in dependancy:
            # handle double species like O2_O+O
            try:
                add_dependency(adding_job, job)  # parent, child
            except RuntimeError:
                pass
# ads_vib jobs dependancies
dependent_workflow_name_2 = facetpath + '_vib'
pending_simulations_dep_2 = BalsamJob.objects.filter(
    workflow__contains=dependent_workflow_name_2
).exclude(state="JOB_FINISHED")

for pending_job in pending_simulations_dep_2:
    for submitted_job in all_submitted_jobs:
        balsam_submitted_job = BalsamJob.objects.filter(
            name=submitted_job
        ).exclude(state="JOB_FINISHED")
        add_dependency(balsam_submitted_job, pending_job)  # parent, child
