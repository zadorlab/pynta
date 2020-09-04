#!/usr/bin/env python3

import os
from pathlib import Path

from rmgcat_to_sella.ts import TS
from rmgcat_to_sella.io import IO

from balsam.launcher.dag import BalsamJob, add_dependency

slab = '{slab}'
repeats = {repeats}
yamlfile = '{yamlfile}'
facetpath = '{facetpath}'
rotAngle = {rotAngle}
scfactor = {scfactor}
scfactor_surface = {scfactor_surface}
pytemplate_xtb = '{pytemplate_xtb}'
species = {species_dict}
current_dir = os.path.dirname(os.getcwd())
minima_dir = os.path.join(facetpath, 'minima')
scaled1 = {scaled1}
scaled2 = {scaled2}
ts_dir = 'TS_estimate'
workflow_name = yamlfile+facetpath+'02'
dependency_workflow_name = yamlfile+facetpath+'01'
creation_dir = '{creation_dir}'

ts = TS(facetpath, slab, ts_dir, yamlfile, repeats, creation_dir)
# ts.copy_minimas_prev_calculated(current_dir, species, minima_dir)
ts.prepare_ts_estimate(scfactor, scfactor_surface, rotAngle,
                       pytemplate_xtb, species, scaled1, scaled2)
all_rxns = IO().get_list_all_rxns_names(yamlfile)

# dependent_workflow_name = yamlfile+facetpath+'03'
pending_simulations_dep = BalsamJob.objects.filter(
    workflow__contains=dependent_workflow_name
).exclude(state="JOB_FINISHED")

pending_simulations = BalsamJob.objects.filter(
    workflow__contains=dependency_workflow_name
).exclude(state="JOB_FINISHED")
cwd = Path.cwd().as_posix()

for rxn in all_rxns:
    dependent_workflow_name = yamlfile+facetpath+'03'+rxn
    path_to_ts_estimate = os.path.join(facetpath, rxn, 'TS_estimate')
    for py_script in Path(path_to_ts_estimate).glob('**/*.py'):
        print(py_script)
        job_dir, script_name = os.path.split(py_script)
        # job_dir = Path.cwd().as_posix() + '/' + '/'.join(
        #     py_script.strip().split('/')[:-1]
        # )
        # script_name = py_script.strip().split('/')[-1]
        # job_to_add = BalsamJob(
        #         name=script_name,
        #         workflow=workflow_name,
        #         application='python',
        #         args=cwd + '/' + py_script,
        #         input_files='',
        #         ranks_per_node=1,
        #         node_packing_count=64,
        #         user_workdir=job_dir,
        #         )
        # job_to_add.save()
        # for job in pending_simulations:
        #     add_dependency(job, job_to_add)  # parent, child
        # for job in pending_simulations_dep:
        #     add_dependency(job_to_add, job)  # parent, child

