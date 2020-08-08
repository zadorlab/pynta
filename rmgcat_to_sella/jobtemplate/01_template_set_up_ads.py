#!/usr/bin/env python3

from glob import glob
from pathlib import Path

from rmgcat_to_sella.adsorbates import Adsorbates

from balsam.launcher.dag import BalsamJob

facetpath = '{facetpath}'
slab = '{slabopt}'
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
    pytemplate, pseudopotentials, pseudo_dir, 
    balsam_exe_settings, calc_keywords
)

cwd = Path.cwd().as_posix()
workflow_name = yamlfile + facetpath + '01'
for py_script in glob('{facetpath}/minima/*.py'):
    job_dir = cwd + '/' + '/'.join(py_script.strip().split('/')[:-1])
    script_name = py_script.strip().split('/')[-1]
    job_to_add = BalsamJob(
            name=script_name,
            workflow=workflow_name,
            application='python',
            args=cwd + '/' + py_script,
            input_files='',
            user_workdir=job_dir,
            node_packing_count=64,
            ranks_per_node=1,
            )
    job_to_add.save()
