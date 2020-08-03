#!/usr/bin/env python3

from rmgcat_to_sella.irc import IRC

slab             = '{slab}'
repeats          = {repeats}
yamlfile         = '{yamlfile}'
facetpath        = '{facetpath}'
pytemplate       = '{pytemplate}'
ts_dir           = 'TS_estimate_unique'
irc_dir          = 'IRC'
pseudopotentials = {pseudopotentials}
pseudo_dir       = '{pseudo_dir}'
workflow_name    = yamlfile+facetpath+'05'
dependency_workflow_name = yamlfile+facetpath+'04' 

irc = IRC(facetpath, slab, repeats, ts_dir, yamlfile,
          pseudopotentials, pseudo_dir)
irc.opt_after_IRC(irc_dir, pytemplate)

from glob import glob
from pathlib import Path

from balsam.launcher.dag import BalsamJob
from balsam.core.models import ApplicationDefinition
BalsamJob = BalsamJob
pending_simulations = BalsamJob.objects.filter(workflow__contains=dependency_workflow_name).exclude(state=“JOB_FINISHED”)
myPython, created= ApplicationDefinition.objects.get_or_create(
            name="Python",
            executable="python3")
myPython.save()
cwd=Path.cwd().as_posix()
for py_script in glob('{facetpath}/IRC/irc*/*.py'):
    #creation_dir=Path.cwd().as_posix()+'/'+'/'.join(py_script.strip().split('/')[:-1])
    job_to_add = BalsamJob(
            name = py_script,
            workflow = workflow_name,
            application = myPython.name,
            args = cwd+py_script,
            ranks_per_node = 1,
            )
    job_to_add.save()
    for job in pending_simulations:
        add_dependency(job,job_to_add) # parent, child

