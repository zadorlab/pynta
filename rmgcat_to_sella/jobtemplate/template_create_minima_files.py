from rmgcat_to_sella.main import WorkFlow
import os
job_file_dir_name = '{job_file_dir_name}'
facetpath = '{facetpath}'
py_job_dir = os.path.join(job_file_dir_name, facetpath)
template_ads_vib = '{template_ads_vib}'
repeats = {repeats}
pytemplate_set_up_ads_vib = '{pytemplate_set_up_ads_vib}'
pseudopotentials = {pseudopotentials}
pseudo_dir = '{pseudo_dir}'
node_packing_count = {node_packing_count}
balsam_exe_settings = {balsam_exe_settings}
calc_keywords = {calc_keywords}
creation_dir = '{creation_dir}'
unique_adsorbates_prefixes = {unique_adsorbates_prefixes}

WorkFlow.set_up_ads_vib(
    template_ads_vib,
    py_job_dir,
    facetpath,
    repeats,
    pytemplate_set_up_ads_vib,
    pseudopotentials,
    pseudo_dir,
    node_packing_count,
    balsam_exe_settings,
    calc_keywords,
    creation_dir,
    unique_adsorbates_prefixes
)
minima_vib_py_script_list = WorkFlow.get_minima_vib_py_scripts(
    facetpath, creation_dir, unique_adsorbates_prefixes)
for minima_vib in minima_vib_py_script_list:
    WorkFlow().exe('', minima_vib, facetpath)
