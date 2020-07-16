#!/usr/bin/env python3

import os
import sys
sys.path.append(os.getcwd())

from rmgcat_to_sella.adsorbates import Adsorbates

facetpath        = '{facetpath}'
slab             = '{slabopt}'
repeats          = {repeats}
yamlfile         = '{yamlfile}'
pytemplate       = '{pytemplate}'
pseudopotentials = {pseudopotentials}
pseudo_dir       = '{pseudo_dir}'
workflow_name    = yamlfile+facetpath+'01'

put_adsorbates = Adsorbates(facetpath, slab, repeats, yamlfile)
put_adsorbates.adjacency_to_3d()
put_adsorbates.create_relax_jobs(pytemplate, pseudopotentials, pseudo_dir)

from glob import glob
from pathlib import Path


from balsam.launcher.dag import BalsamJob
from balsam.core.models import ApplicationDefinition
myPython= ApplicationDefinition.objects.get_or_create(
            name="Python",
            executable="python")
myPython.save()
for py_script in glob('./{facetpath}/minima/*.py'):
    job_to_add = BalsamJob(
            name = py_script,
            workflow = workflow_name,
            application = myPython,
            args = py_script,
            ranks_per_node = 1,
            #data={"creation_dir": Path.cwd().as_posix()+'/{facetpath}/minima'}
            )
    job_to_add.save()
