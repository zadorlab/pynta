#!/usr/bin/env python3
import os

from pathlib import Path

from rmgcat_to_sella.irc import IRC

from balsam.launcher.dag import BalsamJob, add_dependency

slab = '{slab}'
repeats = {repeats}
yamlfile = '{yamlfile}'
facetpath = '{facetpath}'
pytemplate = '{pytemplate}'
ts_estimate_dir = 'TS_estimate_unique'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
rxn = {rxn}
rxn_name = '{rxn_name}'
cwd = Path.cwd().as_posix()
path_to_irc_opt = os.path.join(facetpath, rxn_name, 'IRC')

irc = IRC(
    facetpath,
    slab,
    repeats,
    ts_estimate_dir,
    yamlfile,
    pseudopotentials,
    pseudo_dir,
    balsam_exe_settings,
    calc_keywords,
    creation_dir)

irc.opt_after_IRC(
    rxn,
    pytemplate)

workflow_name = yamlfile + facetpath + '05' + rxn_name
dependency_workflow_name = yamlfile + facetpath + '04' + rxn_name

pending_simulations = BalsamJob.objects.filter(
    workflow__contains=dependency_workflow_name
).exclude(state="JOB_FINISHED")
cwd = Path.cwd().as_posix()

for py_script in Path(path_to_irc_opt).glob('*_opt.py'):
    job_dir, script_name = os.path.split(str(py_script))
    job_to_add = BalsamJob(
        name=script_name,
        workflow=workflow_name,
        application='python',
        args=cwd + '/' + str(py_script),
        input_files='',
        user_workdir=job_dir,
        node_packing_count=48,
        ranks_per_node=1,
    )
    job_to_add.save()
    for job in pending_simulations:
        add_dependency(job, job_to_add)  # parent, child
