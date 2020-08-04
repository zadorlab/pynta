#!/usr/bin/env python3

from rmgcat_to_sella.adsorbates import Adsorbates

facetpath        = '{facetpath}'
slab             = '{slabopt}'
repeats          = {repeats}
yamlfile         = '{yamlfile}'
pytemplate       = '{pytemplate}'
pseudopotentials = {pseudopotentials}
pseudo_dir       = '{pseudo_dir}'
executable       = {executable}
balsam_exe_settings = {balsam_exe_settings}
calc_keywords    = {calc_keywords}

put_adsorbates = Adsorbates(facetpath, slab, repeats, yamlfile)
put_adsorbates.adjacency_to_3d()
put_adsorbates.create_relax_jobs(pytemplate, pseudopotentials, pseudo_dir,executable,balsam_exe_settings,calc_keywords)

from glob import glob
from pathlib import Path


from balsam.launcher.dag import BalsamJob
from balsam.core.models import ApplicationDefinition
myPython, created= ApplicationDefinition.objects.get_or_create(
            name="Python",
            executable="python3")
myPython.save()
cwd=Path.cwd().as_posix()
workflow_name    = yamlfile+facetpath+'01'
for py_script in glob('{facetpath}/minima/*.py'):
    job_to_add = BalsamJob(
            name = py_script,
            workflow = workflow_name,
            application = myPython.name,
            args = cwd+py_script,
            ranks_per_node = 1,
            )
    job_to_add.save()
