#!/usr/bin/env python3
import os
from pathlib import Path

from pynta.ts import TS
from pynta.io import IO

from balsam.launcher.dag import BalsamJob, add_dependency


slab = '{slab}'
repeats = {repeats}
yamlfile = '{yamlfile}'
facetpath = '{facetpath}'
scfactor = {scfactor}
scfactor_surface = {scfactor_surface}
pytemplate_xtb = '{pytemplate_xtb}'
species_list = {species_list}
reacting_atoms = {reacting_atoms}
metal_atom = '{metal_atom}'
scaled1 = {scaled1}
scaled2 = {scaled2}
creation_dir = '{creation_dir}'
rxn = {rxn}
rxn_name = '{rxn_name}'
balsam_exe_settings = {balsam_exe_settings}
minima_dir = os.path.join(creation_dir, facetpath, 'minima')
ts_dir = 'TS_estimate'
path_to_ts_estimate = os.path.join(creation_dir, facetpath, rxn_name, ts_dir)


ts = TS(
    facetpath,
    slab,
    ts_dir,
    yamlfile,
    repeats,
    creation_dir)

ts.prepare_ts_estimate(
    rxn,
    scfactor,
    scfactor_surface,
    pytemplate_xtb,
    species_list,
    reacting_atoms,
    metal_atom,
    scaled1,
    scaled2)

dependancy_dict = IO().depends_on(facetpath, yamlfile, creation_dir)
jobs_to_be_finished = dependancy_dict[rxn_name]

dependency_workflow_name = facetpath + '_01_' + rxn_name
workflow_name = facetpath + '_02_' + rxn_name
dependent_workflow_name = facetpath + '_03_' + rxn_name

pending_simulations_dep = BalsamJob.objects.filter(
    workflow__contains=dependent_workflow_name
).exclude(state="JOB_FINISHED")

# for a given rxn_name, get all BalsamJob objects that it depends on
pending_simulations = []
for dep_job in jobs_to_be_finished:
    pending_simulations.append(BalsamJob.objects.filter(
        name=dep_job).exclude(state="JOB_FINISHED"))

# create BalsamJob objects
for py_script in Path(path_to_ts_estimate).glob('**/*.py'):
    job_dir, script_name = os.path.split(str(py_script))
    job_to_add = BalsamJob(
        name=script_name,
        workflow=workflow_name,
        application='python',
        args=str(py_script),
        input_files='',
        ranks_per_node=1,
        threads_per_rank=balsam_exe_settings['threads_per_rank'],
        node_packing_count={node_packing_count},
        user_workdir=job_dir,
    )
    job_to_add.save()

    # for menten
    # ranks_per_node=1,
    # node_packing_count={node_packing_count},
    # threads_per_rank not specified

    # all job_to_add_ are childs of 01 job, as from jobs_to_be_finished
    # nested for loop becouse BalsamJob.objects.filter(name=dep_job) returns
    # django.query object for a single dep_job, e.g. (H_00_relax.py)
    # no nested loop required if workflow__contains=dependent_workflow_name
    for job in pending_simulations:
        for sub_job in job:
            add_dependency(sub_job, job_to_add)  # parent, child
    # do not run 03 until all 02 for a given reaction are done
    for job in pending_simulations_dep:
        add_dependency(job_to_add, job)  # parent, child
