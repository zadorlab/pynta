#!/usr/bin/env python3
import os
import sys
sys.path.append(os.getcwd())

from rmgcat_to_sella.ts import TS

slab             = '{slab}'
repeats          = {repeats}
yamlfile         = '{yamlfile}'
facetpath        = '{facetpath}'
rotAngle         = {rotAngle}
scfactor         = {scfactor}
scfactor_surface = {scfactor_surface}
pytemplate_xtb   = '{pytemplate_xtb}'
species          = ['{sp1}', '{sp2}']
current_dir      = os.path.dirname(os.getcwd())
minima_dir       = os.path.join(facetpath, 'minima')
scaled1          = {scaled1}
scaled2          = {scaled2}
ts_dir           = 'TS_estimate'
workflow_name    = yamlfile+facetpath+'02'
dependency_workflow_name    = yamlfile+facetpath+'01'

ts = TS(facetpath, slab, ts_dir, yamlfile, repeats)
ts.copy_minimas_prev_calculated(current_dir, species, minima_dir)
ts.prepare_ts_estimate(scfactor, scfactor_surface, rotAngle,
                       pytemplate_xtb, species, scaled1, scaled2)

from glob import glob
from pathlib import Path

from balsam.launcher.dag import BalsamJob
from balsam.core.models import ApplicationDefinition
myPython= ApplicationDefinition.objects.get_or_create(
            name="Python",
            executable="python")
myPython.save()
BalsamJob = BalsamJob
pending_simulations = BalsamJob.objects.filter(workflow__contains=dependency_workflow_name).exclude(state=“JOB_FINISHED”)
for py_script in glob('{facetpath}/TS_estimate/*/*.py'):
    creation_dir=Path.cwd().as_posix()+'/'+'/'.join(py_script.strip().split('/')[:-1])
    job_to_add = BalsamJob(
            name = py_script,
            workflow = workflow_name,
            application = myPython,
            args = py_script,
            ranks_per_node = 1,
#            data={"creation_dir": Path.cwd().as_posix()+'./{facetpath}/minima'}
            )
    job_to_add.save()
    for job in pending_simulations:
        add_dependency(job,job_to_add) # parent, child
