#!/usr/bin/env python3
import os
from pynta.vib import minimaVib
from pathlib import Path
from balsam.launcher.dag import BalsamJob

socket_calculator = '{socket_calculator}'
facetpath = '{facetpath}'
yamlfile = '{yamlfile}'
pytemplate = '{pytemplate}'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
node_packing_count = {node_packing_count}
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
yamlfile_path = os.path.join(creation_dir, yamlfile)
path_to_minima_vib = os.path.join(creation_dir, facetpath, 'minima_vib')

mv = minimaVib(facetpath, creation_dir)
mv.create_minima_vib_all(socket_calculator, facetpath, yamlfile_path,
                         pytemplate, balsam_exe_settings, pseudo_dir,
                         pseudopotentials, calc_keywords,
                         creation_dir)

workflow_name = facetpath + '_01_vib'
for py_script in Path(path_to_minima_vib).glob('**/*.py'):
    job_dir, script_name = os.path.split(str(py_script))
    job_to_add = BalsamJob(
        name=script_name,
        workflow=workflow_name,
        application='python',
        args=str(script_name),
        input_files='',
        user_workdir=job_dir,
        node_packing_count=node_packing_count,
        ranks_per_node=1,
    )
    job_to_add.save()
