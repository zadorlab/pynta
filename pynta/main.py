#!/usr/bin/env python3
from pynta.check_input import InputChecker
from pynta.restart import LowLevelRestart, HighLevelRestart
from pynta.io import IO

from typing import List, Dict, Tuple, Any, Optional
import os
import sys
import __main__
from pathlib import Path, PosixPath
from warnings import warn


# check which file calls this module and adjust working_dir path accordingly
calling_py = os.path.basename(__main__.__file__)
if any([calling_py == i for i in ['run_me.py', 'restart_me.py']]):
    working_dir = os.getcwd()
else:
    working_dir = os.path.dirname(os.path.dirname(os.getcwd()))

# check if all necessary input files are in your working directory
check_yaml = os.path.join(working_dir, 'reactions.yaml')
check_inputR2S = os.path.join(working_dir, 'inputR2S.py')
check_run_me_py = os.path.join(working_dir, 'run_me.py')
check_run_me_sh = os.path.join(working_dir, 'run_me.sh')

InputChecker(check_yaml, check_inputR2S, check_run_me_py,
             check_run_me_sh, working_dir).check_all()

# add working dir to system path
sys.path.insert(1, working_dir)

try:
    import inputR2S
    """
    User defined parameters

    Here we only read them. They are set up in inputR2S.py (submit directory)
    """
except ImportError:
    warn(
        'Missing input file. You cannot run calculations '
        'but will be able to use most of the workflow.'
    )

else:
    optimize_slab = inputR2S.optimize_slab
    surface_types_and_repeats = inputR2S.surface_types_and_repeats
    metal_atom = inputR2S.metal_atom
    a = inputR2S.a
    vacuum = inputR2S.vacuum
    pseudo_dir = inputR2S.pseudo_dir
    pseudopotentials = inputR2S.pseudopotentials
    yamlfile = inputR2S.yamlfile
    scfactor = inputR2S.scfactor
    scfactor_surface = inputR2S.scfactor_surface
    scaled1 = inputR2S.scaled1
    scaled2 = inputR2S.scaled2
    all_reacting_atoms = IO.get_all_reacting_atoms(check_yaml)
    species_dict = IO().get_species_dict(check_yaml)
    all_species = IO().get_all_unique_species(check_yaml)
    executable = inputR2S.executable
    node_packing_count = inputR2S.node_packing_count
    balsam_exe_settings = inputR2S.balsam_exe_settings
    calc_keywords = inputR2S.calc_keywords
    creation_dir = inputR2S.creation_dir
    surface_types = surface_types_and_repeats.keys()
    facetpaths = IO().get_facetpaths(metal_atom, surface_types)
    job_file_dir_name = 'job_files'

####################################################
#                    Scripts                       #
####################################################

# These template and pytemplate scripts can be modified by users
# (np intended, though) to tune them to given calculation setup,
# i.e. calculator, method, queue system,
# etc. The current version works for SLURM and Quantum Espresso.

path = os.path.abspath(__file__)
dir_path = os.path.dirname(path)
path_template = os.path.join(dir_path, 'jobtemplate/')
path_pytemplate = os.path.join(dir_path, 'pytemplate/')
template_slab_opt = os.path.join(
    path_template + '00_template_set_up_slab_opt.py')
template_big_slab_opt = os.path.join(
    path_template + '00_template_set_up_big_slab_opt.py')
template_ads = os.path.join(path_template + '01_template_set_up_ads.py')
template_ads_vib = os.path.join(
    path_template + '01_template_set_up_ads_vib.py')
template_set_up_ts_with_xtb = os.path.join(
    path_template + '02_template_set_up_ts_with_xtb.py')
template_set_up_ts = os.path.join(
    path_template + '03_template_checksym_xtb_runTS.py')
template_set_up_ts_vib = os.path.join(
    path_template + '04_template_set_up_TS_vib.py')
template_set_up_after_ts = os.path.join(
    path_template + '05_template_set_up_after_ts.py')
pytemplate_big_slab_opt = os.path.join(
    path_pytemplate + 'pytemplate_set_up_big_slab_opt.py')
pytemplate_relax_ads = os.path.join(
    path_pytemplate + 'pytemplate_set_up_ads_on_slab.py')
pytemplate_set_up_ads_vib = os.path.join(
    path_pytemplate + 'pytemplate_set_up_ads_vib.py')
pytemplate_xtb = os.path.join(path_pytemplate + 'pytemplate_set_up_xtb.py')
pytemplate_set_up_ts = os.path.join(
    path_pytemplate + 'pytemplate_set_up_ts.py')
pytemplate_set_up_ts_vib = os.path.join(
    path_pytemplate + 'pytemplate_set_up_ts_vib.py')
pytemplate_set_up_after_ts = os.path.join(
    path_pytemplate + 'pytemplate_set_up_opt_after_ts.py')
slab_opt = '00_set_up_slab_opt.py'

####################################################
#                    Initialize                    #
####################################################


class WorkFlow:

    def __init__(self):
        ''' Setup the balsam application for this workflow run.

            Once we start using QE will want one app for QE,
            one for xtb most likely
        '''
        print('Checking Balsam DB...')
        try:
            from balsam.core.models import ApplicationDefinition

            self.myPython, _ = ApplicationDefinition.objects.get_or_create(
                name="python",
                executable=sys.executable
            )
            self.myPython.save()
            self.slab_opt_job = ''

            # TODO: instead of directly importing EspressoBalsam, we should
            # write a function which returns the appropriate class from
            # balsamcalc.py based on the user-provided input file
            from pynta.balsamcalc import (
                EspressoBalsam, EspressoBalsamSocketIO
            )
            EspressoBalsam.exe = executable
            EspressoBalsamSocketIO.exe = executable
            EspressoBalsam.create_application()
            EspressoBalsamSocketIO.create_application()
        except SystemExit:
            print('---')
            print('Please create Balsam DB and/or activate it')
            print('---')

    @staticmethod
    def get_ts_xtb_py_script_list(
            facetpath: str) -> List[str]:
        '''Get a list with all 02 job scripts

        Parameters
        ----------
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        Returns
        -------
        ts_with_xtb_py_script_list : list(str)
            a list with all ts_with_xtb_py jobs for a given facetpath

        '''
        reactions = IO().open_yaml_file(yamlfile)
        ts_with_xtb_py_script_list = []
        for rxn in reactions:
            rxn_name = IO.get_rxn_name(rxn)
            fname = '02_{}_set_up_TS_with_xtb_{}.py'.format(
                facetpath, rxn_name)
            ts_with_xtb_py_script_list.append(fname)
        return ts_with_xtb_py_script_list

    @staticmethod
    def get_ts_estimate_unique_list(
            facetpath: str) -> List[str]:
        ''' Get a list with all 03 job scripts

        Parameters:
        ----------
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        Returns:
        -------
        ts_sella_py_script_list : list(str)
            a list with all ts_sella_py jobs for a given facetpath

        '''
        reactions = IO().open_yaml_file(yamlfile)
        ts_sella_py_script_list = []
        for rxn in reactions:
            rxn_name = IO.get_rxn_name(rxn)
            fname = '03_{}_run_TS_{}.py'.format(facetpath, rxn_name)
            ts_sella_py_script_list.append(fname)
        return ts_sella_py_script_list

    @staticmethod
    def get_ts_vib_list(
            facetpath: str) -> List[str]:
        ''' Get a list with all 04 job scripts

        Parameters:
        ----------
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        Returns:
        -------
        ts_vib_py_scripts_list : list(str)
            a list with all ts_vib_py jobs for a given facetpath

        '''
        reactions = IO().open_yaml_file(yamlfile)
        ts_vib_py_scripts_list = []
        for rxn in reactions:
            rxn_name = IO.get_rxn_name(rxn)
            fname = '04_{}_set_up_TS_vib_{}.py'.format(facetpath, rxn_name)
            ts_vib_py_scripts_list.append(fname)
        return ts_vib_py_scripts_list

    @staticmethod
    def get_after_ts_py_scripts(
            facetpath: str) -> List[str]:
        ''' Get a list with all 05 job scripts

        Parameters:
        ----------
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        Returns:
        -------
        after_ts_py_scripts_list : list(str)
            a list with all after_ts_py jobs for a given facetpath

        '''
        reactions = IO().open_yaml_file(yamlfile)
        after_ts_py_scripts_list = []
        for rxn in reactions:
            rxn_name = IO.get_rxn_name(rxn)
            fname = '05_{}_set_up_after_TS_{}.py'.format(facetpath, rxn_name)
            after_ts_py_scripts_list.append(fname)
        return after_ts_py_scripts_list

    @staticmethod
    def create_big_slab_pyjob(
            pytemplate: str,
            facetpath: str,
            slab_name: str,
            repeats: Tuple[int, int, int],
            creation_dir: PosixPath) -> None:
        ''' Create a ase pyjob for big slab optimization

        Parameters
        ----------
        pytemplate : python file
            a template to prepare submission scripts
            for big_slab minimization
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        slab_name : str
            name of the slab file (no extintion) that is about to be created
            e.g.
            slab_name = 'Cu_111_slab_opt'
        repeats : tuple(int, int, int)
            how to replicate unit cell in (x, y, z) direction,
            e.g. (3, 3, 1)
        creation_dir : posix
            a posix path to the working directory

        '''
        with open(pytemplate, 'r') as r:
            pytemplate_text = r.read()
            py_job_name = '{}_big_slab_opt_job.py'.format(facetpath)
            py_job = os.path.join(
                creation_dir, job_file_dir_name, facetpath, py_job_name)
            with open(py_job, 'w') as c:
                c.write(pytemplate_text.format(
                    facetpath=facetpath,
                    slab_name=slab_name,
                    repeats=repeats,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    creation_dir=creation_dir
                ))

    @staticmethod
    def create_job_files() -> None:
        ''' For each surface type and for each reaction
            generate submit scripts for 6 stages of the workflow

        '''
        for surface_type, reps in surface_types_and_repeats.items():
            # get a facetpath
            facetpath = IO().get_facetpath(metal_atom, surface_type)
            # Create a directory to store six (00-05) main *py job files
            py_job_dir = os.path.join(job_file_dir_name, facetpath)
            os.makedirs(py_job_dir, exist_ok=True)

            # naming convetion for the slab file
            slab_name = facetpath + '_slab_opt'
            slab = slab_name + '.xyz'

            # define repeats and repeats_surface
            repeats, repeats_surface = reps

            # create all input files for each main job 00-05
            WorkFlow.set_up_slab(
                template_slab_opt,
                py_job_dir,
                facetpath,
                surface_type,
                repeats_surface,
                slab_name
            )

            WorkFlow.set_up_big_slab(
                template_big_slab_opt,
                py_job_dir,
                facetpath,
                slab_name,
                repeats,
                pytemplate_big_slab_opt
            )

            WorkFlow.set_up_ads(
                template_ads,
                py_job_dir,
                facetpath,
                slab,
                repeats,
                yamlfile,
                pytemplate_relax_ads
            )

            WorkFlow.set_up_ads_vib(
                template_ads_vib,
                py_job_dir,
                facetpath,
                pytemplate_set_up_ads_vib
            )

            reactions = IO().open_yaml_file(yamlfile)
            # the rest of jobs depends on reaction and reacting_atoms
            for rxn, react_at in zip(reactions, all_reacting_atoms.values()):
                WorkFlow.set_up_TS_with_xtb(
                    rxn,
                    react_at,
                    template_set_up_ts_with_xtb,
                    py_job_dir,
                    slab,
                    repeats,
                    facetpath,
                    scfactor,
                    scfactor_surface,
                    pytemplate_xtb
                )

                WorkFlow.set_up_run_TS(
                    rxn,
                    template_set_up_ts,
                    py_job_dir,
                    facetpath,
                    slab,
                    repeats,
                    pytemplate_set_up_ts
                )

                WorkFlow.set_up_TS_vib(
                    rxn,
                    template_set_up_ts_vib,
                    py_job_dir,
                    facetpath,
                    slab,
                    repeats,
                    pytemplate_set_up_ts_vib
                )

                WorkFlow.set_up_opt_after_TS(
                    rxn,
                    template_set_up_after_ts,
                    py_job_dir,
                    facetpath,
                    slab,
                    repeats,
                    pytemplate_set_up_after_ts
                )

###########################
#   Create submit files   #
###########################
    @staticmethod
    def set_up_slab(
            template: str,
            py_job_dir: str,
            facetpath: str,
            surface_type: str,
            repeats_surface: Tuple[int, int, int],
            slab_name: str) -> None:
        ''' Create 00_{facetpath}_set_up_slab_opt.py file

        Parameters
        ----------
        template : py file
            a template for 00 job (slab optimization)
        py_job_dir : str
            a path to where all job files 00-05 are about to be created, e.g.
            {'current_dir'}/job_files/Cu_111
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        surface_type : str
            type of the surface. Available options are:
            fcc111, fcc211, fcc100, bcc111, bcc110, hcp0001, diamond111,
            diamond100
        metal_atom : str
            atomic metal_atom of the studied metal surface
            e.g. 'Cu'
        a : float
            a lattice constant
        repeats_surface : tuple(int, int, int)
            surface multiplication in (x, y, z) direction
             eg. (1, 1, 4)
        vacuum : float
            amout of empty space in z direction (Angstrem)
        slab_name : str
            name of the slab file (no extintion) that is about to be created
            e.g.
            slab_name = 'Cu_111_slab_opt'
        pseudopotentials : dict{str:str}
            a dictionary with QE pseudopotentials for all species.
            e.g.
            dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                H='H.pbe-kjpaw_psl.1.0.0.UPF',
                O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                C='C.pbe-n-kjpaw_psl.1.0.0.UPF',
                )
        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.
            '/home/mgierad/espresso/pseudo'
        balsam_exe_settings : dict{str:int}
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.
            balsam_exe_settings = {'num_nodes': 1,
                                   'ranks_per_node': 48,
                                   'threads_per_rank': 1}
        calc_keywords : dict{str:str}
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            calc_keywords = {'kpts': (3, 3, 1), 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}
        creation_dir : posix
            a posix path to the working directory

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            py_job_fname = os.path.join(
                py_job_dir, '00_{}_set_up_slab_opt.py'.format(facetpath))
            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    surface_type=surface_type,
                    metal_atom=metal_atom,
                    a=a,
                    repeats_surface=repeats_surface,
                    vacuum=vacuum,
                    slab_name=slab_name,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=creation_dir
                ))

    @staticmethod
    def set_up_big_slab(
            template: str,
            py_job_dir: str,
            facetpath: str,
            slab_name: str,
            repeats: Tuple[int, int, int],
            pytemplate: str) -> None:
        ''' Set up a big_slab_opt ase py job generator

        Parameters
        ----------
        template : py file
            a template to set up 01 job
        py_job_dir : str
            a path to where all job files 00-05 are about to be created, e.g.
            {'current_dir'}/job_files/Cu_111
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        slab_name : str
            name of the slab file (no extintion) that is about to be created
            e.g.
            slab_name = 'Cu_111_slab_opt'
        repeats : tuple(int, int, int)
            how to replicate unit cell in (x, y, z) direction,
            e.g. (3, 3, 1)
        pytemplate : python file
            a template to prepare submission scripts
            for big_slab minimization

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            py_job_fname = os.path.join(
                py_job_dir, '00_{}_set_up_big_slab_opt.py'.format(facetpath))

            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    facetpath=facetpath,
                    slab_name=slab_name,
                    repeats=repeats,
                    pytemplate=pytemplate,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=creation_dir,
                    job_file_dir_name=os.path.join(
                        creation_dir, job_file_dir_name, facetpath),
                    node_packing_count=node_packing_count
                ))

    @staticmethod
    def set_up_ads(
            template,
            py_job_dir,
            facetpath,
            slab,
            repeats,
            yamlfile,
            pytemplate):
        ''' Create 01_{facetpath}_set_up_ads_on_slab_{rxn_name}.pyfile

        Parameters:
        ___________
        template : py file
            a template to set up 01 job
        py_job_dir : str
            a path to where all job files 00-05 are about to be created, e.g.
            {'current_dir'}/job_files/Cu_111
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        slab : str
            a path to .xyz file with optimized slab
        repeats : tuple(int, int, int)
            how to replicate unit cell in (x, y, z) direction,
            e.g. (3, 3, 1)
        yamlfile : str
            a name of the .yaml file with a reaction list
        pytemplate : python file
            a template to prepare submission scripts
            for adsorbate+surface minimization

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            py_job_fname = os.path.join(
                py_job_dir, '01_{}_set_up_ads_on_slab.py'.format(facetpath))

            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    facetpath=facetpath,
                    slab=slab,
                    yamlfile=yamlfile,
                    repeats=repeats,
                    pytemplate=pytemplate,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    node_packing_count=node_packing_count,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=creation_dir
                ))

    @staticmethod
    def set_up_ads_vib(
            template: str,
            py_job_dir: str,
            facetpath: str,
            pytemplate: str) -> None:
        '''[summary]

        Parameters
        ----------
        template : str
            a template file for ads_vib job
        py_job_dir : str
            a path to where all job files 00-05 are about to be created, e.g.
            {'current_dir'}/job_files/Cu_111
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        pytemplate : str
            a template to prepare submission scripts for vibfrational
            frequencies calculations of symmetry distinct minima

        '''

        with open(template, 'r') as f:
            template_txt = f.read()
            py_job_fname = os.path.join(
                creation_dir,
                py_job_dir,
                '{}_set_up_ads_vib.py'.format(facetpath))
            with open(py_job_fname, 'w') as c:
                c.write(template_txt.format(
                    facetpath=facetpath,
                    yamlfile=yamlfile,
                    pytemplate=pytemplate,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    node_packing_count=node_packing_count,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=creation_dir
                ))

    @staticmethod
    def set_up_TS_with_xtb(
            rxn: Dict[str, str],
            reacting_atoms: List[str],
            template: str,
            py_job_dir: str,
            slab: str,
            repeats: Tuple[int, int, int],
            facetpath: str,
            scfactor: float,
            scfactor_surface: float,
            pytemplate_xtb: str) -> None:
        ''' Create 02_{facetpath}_set_up_TS_with_xtb_{rxn_name}.py files

        Parameters:
        ___________

        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        template : py file
            a template to set up 02 job for a particular reaction
        py_job_dir : str
            a path to where all job files 00-05 are about to be created, e.g.
            {'current_dir'}/job_files/Cu_111
        slab : str
            a '.xyz' file name with the optimized slab
            e.g.
            'Cu_111_slab_opt.xyz'
        repeats : tuple(int, int, int)
            how to replicate unit cell in (x, y, z) direction,
            e.g. (3, 3, 1)
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        scfator : float
            a scaling factor to scale a bond distance between
            atoms taking part in the reaction
            e.g. 1.4
        scfactor_surface : float
            a scaling factor to scale the target bond distance, i.e.
            the average distance between adsorbed atom and the nearest
            surface atom. Helpful e.g. when H is far away form the surface
            in TS, whereas for minima it is close to the surface
            e.g. 1.0
        pytemplate_xtb : python script
            a template file for penalty function minimization job

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            rxn_name = IO.get_rxn_name(rxn)
            rxn_no = rxn['index']

            py_job_fname = os.path.join(
                py_job_dir, '02_{}_set_up_TS_with_xtb_{}.py'.format(
                    facetpath, rxn_name))

            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    facetpath=facetpath,
                    slab=slab,
                    repeats=repeats,
                    yamlfile=yamlfile,
                    scfactor=scfactor,
                    scfactor_surface=scfactor_surface,
                    pytemplate_xtb=pytemplate_xtb,
                    species_list=species_dict['rxn' + str(rxn_no)],
                    metal_atom=metal_atom,
                    scaled1=scaled1,
                    scaled2=scaled2,
                    reacting_atoms=reacting_atoms,
                    creation_dir=creation_dir,
                    rxn=rxn,
                    rxn_name=rxn_name,
                    node_packing_count=node_packing_count
                ))

    @staticmethod
    def set_up_run_TS(
            rxn: Dict[str, str],
            template: str,
            py_job_dir: str,
            facetpath: str,
            slab: str,
            repeats: Tuple[int, int, int],
            pytemplate: str) -> None:
        ''' Create 03_{facetpath}_set_up_run_TS_{rxn_name}.py file

        Parameters:
        ___________

        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        template : py file
            a template to set up 03 job for a particular reaction
        py_job_dir : str
            a path to where all job files 00-05 are about to be created, e.g.
            {'current_dir'}/job_files/Cu_111
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        slab : str
            a '.xyz' file name with the optimized slab
            e.g.
            'Cu_111_slab_opt.xyz'
        repeats : tuple(int, int, int)
            how to replicate unit cell in (x, y, z) direction,
            e.g. (3, 3, 1)
        pytemplate : python script
            a template file for ts optimization with sella

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            rxn_name = IO.get_rxn_name(rxn)
            py_job_fname = os.path.join(
                py_job_dir, '03_{}_run_TS_{}.py'.format(
                    facetpath, rxn_name))

            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    facetpath=facetpath,
                    slab=slab,
                    repeats=repeats,
                    yamlfile=yamlfile,
                    pytemplate=pytemplate,
                    pseudo_dir=pseudo_dir,
                    pseudopotentials=pseudopotentials,
                    node_packing_count=node_packing_count,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=creation_dir,
                    rxn=rxn,
                    rxn_name=rxn_name
                ))

    @staticmethod
    def set_up_TS_vib(
            rxn: Dict[str, str],
            template: str,
            py_job_dir: str,
            facetpath: str,
            slab: str,
            repeats: Tuple[int, int, int],
            pytemplate: str) -> None:
        ''' Create '04_{facetpath}_set_up_TS_vib_{rxn_name}.py file

        Parameters:
        ___________

        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        template : py file
            a template to set up 03 job for a particular reaction
        py_job_dir : str
            a path to where all job files 00-05 are about to be created, e.g.
            {'current_dir'}/job_files/Cu_111
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        slab : str
            a '.xyz' file name with the optimized slab
            e.g.
            'Cu_111_slab_opt.xyz'
        repeats : tuple(int, int, int)
            how to replicate unit cell in (x, y, z) direction,
            e.g. (3, 3, 1)
        pytemplate : python script
            a template file for setting up frequencies calculations

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            rxn_name = IO.get_rxn_name(rxn)
            py_job_fname = os.path.join(
                py_job_dir, '04_{}_set_up_TS_vib_{}.py'.format(
                    facetpath, rxn_name))

            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    facetpath=facetpath,
                    slab=slab,
                    repeats=repeats,
                    yamlfile=yamlfile,
                    pytemplate=pytemplate,
                    pseudo_dir=pseudo_dir,
                    pseudopotentials=pseudopotentials,
                    node_packing_count=node_packing_count,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=creation_dir,
                    rxn=rxn,
                    rxn_name=rxn_name
                ))

    @staticmethod
    def set_up_opt_after_TS(
            rxn: Dict[str, str],
            template: str,
            py_job_dir: str,
            facetpath: str,
            slab: str,
            repeats: Tuple[int, int, int],
            pytemplate: str) -> None:
        ''' Create 05_{facetpath}_set_up_after_TS_{rxn_name}.py file

        Parameters:
        ___________

        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        template : py file
            a template to set up 03 job for a particular reaction
        py_job_dir : str
            a path to where all job files 00-05 are about to be created, e.g.
            {'current_dir'}/job_files/Cu_111
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        slab : str
            a '.xyz' file name with the optimized slab
            e.g.
            'Cu_111_slab_opt.xyz'
        repeats : tuple(int, int, int)
            how to replicate unit cell in (x, y, z) direction,
            e.g. (3, 3, 1)
        pytemplate : python script
            a template file for setting up an alternative to IRC (minimization
            of displaced structures following the imaginary mode of
            oscillation)

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            rxn_name = IO.get_rxn_name(rxn)
            py_job_fname = os.path.join(
                py_job_dir, '05_{}_set_up_after_TS_{}.py'.format(
                    facetpath, rxn_name))

            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    facetpath=facetpath,
                    slab=slab,
                    repeats=repeats,
                    yamlfile=yamlfile,
                    pytemplate=pytemplate,
                    pseudo_dir=pseudo_dir,
                    pseudopotentials=pseudopotentials,
                    node_packing_count=node_packing_count,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=creation_dir,
                    rxn=rxn,
                    rxn_name=rxn_name
                ))

##############################
# Submit jobs and execute it #
##############################

    def exe(
            self,
            parent_job: str,
            job_script: str,
            facetpath: str,
            cores: int = 1) -> Any:
        ''' Execute a py script

        Parameters:
        ___________

        parent_job : str
            a parent job on which subbmited jobs depends on. Formatted as:
            00 or 01 or 02 or 03 or 04, depending on the job
        job_script : str
            a script that is about to be submitted
        cores : int
            number of cores for exe job
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        Returns:
        ________

        job_to_add : BalsamJob
            job that will be submitted to balsam database

        '''
        from balsam.launcher.dag import BalsamJob

        job_files_path = os.path.join(
            creation_dir, job_file_dir_name, facetpath)

        # get rxn_name from job_script by spliting and joining job_script name
        # exeption for two first jobs
        slab_opt = '00_{}_set_up_slab_opt.py'.format(facetpath)
        big_slab_opt = '00_{}_set_up_big_slab_opt.py'.format(facetpath)
        ads_opt = '01_{}_set_up_ads_on_slab.py'.format(facetpath)
        ads_vib = '{}_set_up_ads_vib.py'.format(facetpath)

        # specify rxn_name TODO a better method required!
        if job_script in [ads_opt, slab_opt, ads_vib]:
            rxn_name = ''
        else:
            rxn_name = '_'.join(job_script.split('_')[-2:])[:-3]

        # specify workflow name
        if job_script == big_slab_opt:
            # special case for big_slab_opt jobs
            workflow_name = facetpath + '_big_slab_opt'
        elif job_script == ads_vib:
            # special case for ads_vib jobs
            workflow_name = facetpath + '_vib'
        else:
            # default name for the workflow
            try:
                workflow_name = facetpath + '_' + \
                    job_script[0:2] + '_' + rxn_name
            except ValueError:
                workflow_name = facetpath + '_error'

        job = os.path.join(job_files_path, job_script)
        job_to_add = BalsamJob(
            name=job_script,
            workflow=workflow_name,
            application=self.myPython.name,
            args=job,
            ranks_per_node=cores,
            input_files='',
            node_packing_count=node_packing_count,
            user_workdir=job_files_path
        )
        job_to_add.save()

        # if there is a parent job, specify dependency
        if parent_job != '':
            from balsam.launcher.dag import add_dependency
            try:
                add_dependency(parent_job, job_to_add)  # parent, child
            except ValueError:
                dependency = str(parent_job[0:2])

                if parent_job == '01':
                    # a special case for 01 where there is on job script
                    # for all reactions for a given facet
                    dependency_workflow_name = os.path.join(
                        facetpath + '_' + dependency + '_')
                else:
                    dependency_workflow_name = os.path.join(
                        facetpath + '_' + dependency + '_' + rxn_name)

                BalsamJob = BalsamJob
                pending_simulations = BalsamJob.objects.filter(
                    workflow__contains=dependency_workflow_name
                ).exclude(state='JOB_FINISHED')

                for job in pending_simulations:
                    add_dependency(job, job_to_add)  # parent, child

        return job_to_add

    def run_slab_optimization(
            self,
            facetpath: str) -> None:
        ''' Submit slab_optimization_job

        Parameters:
        ___________

        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        slab_opt = '00_{}_set_up_slab_opt.py'.format(facetpath)
        self.slab_opt_job = self.exe('', slab_opt, facetpath, cores=1)

    def run_big_slab_opt(
            self,
            facetpath: str) -> None:
        ''' Submit big slab optimization

        Parameters
        ----------
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        big_slab_opt = '00_{}_set_up_big_slab_opt.py'.format(facetpath)
        return self.exe('', big_slab_opt, facetpath, cores=1)

    def run_opt_surf_and_adsorbate(
            self,
            facetpath: str) -> None:
        ''' Run optmization of adsorbates on the surface

        Parameters:
        ___________

        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        ads_surf_opt_script = '01_{}_set_up_ads_on_slab.py'.format(facetpath)
        slab_opt = '00_{}_set_up_slab_opt.py'.format(facetpath)
        return self.exe(slab_opt, ads_surf_opt_script, facetpath)

    def run_opt_surf_and_adsorbate_no_depend(
            self,
            facetpath: str) -> None:
        ''' Run optmization of adsorbates on the surface
            if there is no dependency on other jobs

        Parameters:
        ___________

        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        ads_surf_opt_script = '01_{}_set_up_ads_on_slab.py'.format(facetpath)
        return self.exe('', ads_surf_opt_script, facetpath)

    def run_minima_vib(
            self,
            dependent_job: str,
            facetpath: str) -> None:
        ''' Run vibrational frequnecies calculations for the lowest energy
            conformer of a given adsorbate

        Parameters
        ----------
        dependant_job : str
            a prefix of the the job (minimum: two integers,
            max: whole job script name) that the current run depends on,
            e.g. '01'
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        ads_vib_script = '{}_set_up_ads_vib.py'.format(facetpath)
        return self.exe(dependent_job, ads_vib_script, facetpath)

    def run_minima_vib_no_depend(
            self,
            facetpath: str) -> None:
        ads_vib_script = '{}_set_up_ads_vib.py'.format(facetpath)
        return self.exe('', ads_vib_script, facetpath)

    def run_ts_estimate(
            self,
            dependent_job: str,
            facetpath: str) -> None:
        ''' Run TS estimation calculations

        Parameters:
        ___________

        dependant_job : str
            a prefix of the the job (minimum: two integers,
            max: whole job script name) that the current run depends on,
            e.g. '01'
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        ts_xtb_py_script_list = WorkFlow.get_ts_xtb_py_script_list(facetpath)
        for ts_xtb in ts_xtb_py_script_list:
            self.exe(dependent_job, ts_xtb, facetpath)

    def run_ts_estimate_no_depend(
            self,
            facetpath: str) -> None:
        ''' Run TS estimate calculations if there is no dependency on the
            other jobs

        Parameters:
        ___________

        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        ts_xtb_py_script_list = WorkFlow.get_ts_xtb_py_script_list(facetpath)
        for ts_xtb in ts_xtb_py_script_list:
            self.exe('', ts_xtb, facetpath)

    def run_ts_with_sella(
            self,
            dependant_job: str,
            facetpath: str) -> None:
        ''' Run TS minimization with Sella

        Parameters:
        ___________

        dependant_job : str
            a prefix of the the job (minimum: two integers,
            max: whole job script name) that the current depends on,
            e.g. '02'
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        ts_sella_py_script_list = WorkFlow.get_ts_estimate_unique_list(
            facetpath)
        for ts_sella in ts_sella_py_script_list:
            self.exe(dependant_job, ts_sella, facetpath)

    def run_ts_vib(
            self,
            dependant_job: str,
            facetpath: str) -> None:
        ''' Run frequency calculations for TS

        Parameters:
        ___________

        dependant_job : str
            a prefix of the the job (minimum: two integers,
            max: whole job script name) that the current depends on,
            e.g. '03'
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        ts_vib_py_script_list = WorkFlow.get_ts_vib_list(facetpath)
        for ts_vib in ts_vib_py_script_list:
            self.exe(dependant_job, ts_vib, facetpath)

    def run_opt_after_ts(
            self,
            dependant_job: str,
            facetpath: str) -> None:
        ''' Run minimization of minima obtained as nudging TS structure
            towards imaginary mode of oscilation

        Parameters:
        ___________

        dependant_job : str
            a prefix of the the job (minimum: two integers,
            max: whole job script name) that the current depends on,
            e.g. '04'
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        after_irc_py_scripts = WorkFlow.get_after_ts_py_scripts(facetpath)
        for after_irc in after_irc_py_scripts:
            self.exe(dependant_job, after_irc, facetpath)

    @staticmethod
    def check_all_species(
            yamlfile: str,
            facetpath: str) -> Dict[str, bool]:
        ''' Check all species(all reactions) to find whether
            there are previous calculation the code can use

        Parameters:
        ___________

        yamlfile: str
            a name of the .yaml file with a reaction list
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        Return:
        _______
        checked_species: dict(str: bool)
            a dictionary with True/False values for every species(keys)
            True if there are previous calculation, otherwise False
            e.g. {'C': True, 'H': True, 'O': False, 'OH': True, 'CH': True}

        '''
        checked_species = {}
        all_species = IO().get_all_unique_species(yamlfile)
        for species in all_species:
            checked_species[species] = WorkFlow.is_minima_out_files(
                species, facetpath)
        return checked_species

    @staticmethod
    def is_minima_dir(
            species: str,
            facetpath: str) -> bool:
        ''' Return True if directory for a given species exists in
            {facepath}/minima. Otherwise, False

        Parameters:
        ___________

        species: str
            a species metal_atom
            e.g. 'H' or 'CO'
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        species_minima_dir = os.path.join(
            creation_dir, facetpath, 'minima', species)
        if os.path.isdir(species_minima_dir):
            return True
        return False

    @staticmethod
    def is_minima_out_files(
            species: str,
            facetpath: str) -> bool:
        ''' Check for the previously calculated * relax.out files for a given
            species. Return True if there are previous calculations. Otherwise,
            False.

        Parameters:
        ___________

        species: str
            a species metal_atom
            e.g. 'H' or 'CO'
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        minima_dir = os.path.join(creation_dir, facetpath, 'minima')
        # if minima dir exists, check for outfiles
        if WorkFlow.is_minima_dir(species, facetpath):
            keyphrase = '{}_{}_*relax.py.out'.format(facetpath, species)
            search_for_outfiles = Path(minima_dir).glob(keyphrase)
            outfiles = []
            for outfile in search_for_outfiles:
                outfiles.append(str(outfile))
            # empty list
            if not outfiles:
                return False
            return True
        return False

    @staticmethod
    def is_slab(
            facetpath: str) -> Tuple[bool, Optional[str]]:
        ''' Check whether slab has been already optimized

        Parameters:
        ___________

        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        Returns:
        ________

        tuple(bool, str=None):
            True if there are previous calculations
                (True, path_to_prev_calc)
            False otherwise
                (False, )

        '''
        slab_opt_path_str = []
        # the code will look for anything like Cu_111*.xyz starting from the
        # facetpath directory including all subdirectories.
        keyphrase = os.path.join(facetpath + '_slab_opt.xyz')
        slab_opt_path_posix = Path(str(os.getcwd())).glob(keyphrase)
        for slab_opt_path in slab_opt_path_posix:
            slab_opt_path_str.append(slab_opt_path)
        if len(slab_opt_path_str) >= 1:
            return True, slab_opt_path_str[0]
        return (False, )

    @staticmethod
    def is_big_slab(
            facetpath: str) -> bool:
        ''' Check for big_slab calculations. True if there is a big_slab file,
            False if not. If multiple matches found, print all and raise error.

        big_slab is a slab multiplied by repeats.

        Parameters:
        ___________
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        big_slab_list = []
        keyphrase = '{}_big_slab_opt*.xyz'.format(facetpath)
        big_slab_paths = Path(creation_dir).glob(keyphrase)

        for big_slab_path in big_slab_paths:
            big_slab_list.append(big_slab_path)

        if len(big_slab_list) > 1:
            print('Multiple matches found. Please check the work dir')
            print(big_slab_list)
            exit()
        elif len(big_slab_list) == 1:
            return True
        else:
            print('No matches found for {} '
                  'Big slab optimization required'.format(keyphrase))
            return False

    @staticmethod
    def check_if_slab_opt_exists(
            facetpath: str) -> Tuple[bool, Optional[str]]:
        ''' Check whether slab has been already optimized

        Parameters:
        ___________

        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        Returns:
        ________

        tuple(bool, str=None):
            True if there are previous calculations
                (True, path_to_prev_calc)
            False otherwise
                (False, )

        '''
        slab_opt_path_str = []
        # the code will look for anything like Cu_111*.xyz starting from the
        # facetpath directory including all subdirectories.
        keyphrase = str(facetpath) + '*.xyz'
        slab_opt_path_posix = Path(str(os.getcwd())).glob(keyphrase)
        for slab_opt_path in slab_opt_path_posix:
            slab_opt_path_str.append(slab_opt_path)
        if len(slab_opt_path_str) >= 1:
            return True, slab_opt_path_str[0]
        return (False, )

    def execute(
            self,
            facetpath: str) -> None:
        ''' The main execute method for a given surface

        Parameters:
        ___________

        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'

        '''
        if optimize_slab:
            # if slab found in previous calculation, do nothing
            if WorkFlow.is_slab(facetpath)[0] is False:
                # If the code cannot locate optimized slab .xyz file,
                # a slab optimization will be launched.
                self.run_slab_optimization(facetpath)
            if WorkFlow.is_big_slab(facetpath) is False:
                self.run_big_slab_opt(facetpath)
            # check if species were already calculated
            # TODO there is a bug here as it counts CO as C
            if all(WorkFlow.check_all_species(yamlfile, facetpath).values()):
                # If all are True, start by generating TS guesses and run
                # the penalty function minimization
                self.run_minima_vib_no_depend(facetpath)
                self.run_ts_estimate_no_depend(facetpath)
            else:
                # If any of sp_check_list is False
                # run optimization of surface + reactants; surface + products
                try:
                    self.run_opt_surf_and_adsorbate(facetpath)
                except NameError:
                    self.run_opt_surf_and_adsorbate_no_depend(facetpath)
                self.run_minima_vib('01', facetpath)
                self.run_ts_estimate('01', facetpath)
        else:
            # this is executed if user provide .xyz with the optimized slab
            # and explicitly define oiptimize_slab = False
            if WorkFlow.check_if_slab_opt_exists(facetpath)[0]:
                pass
            else:
                raise FileNotFoundError(
                    'It appears there is no slab_opt.xyz file'
                )
            if WorkFlow.is_big_slab(facetpath) is False:
                self.run_big_slab_opt(facetpath)
            if all(WorkFlow.check_all_species(yamlfile, facetpath).values()):
                # If all minima were calculated some time age pynta
                # will use that calculations. Start from 02 step
                self.run_ts_estimate_no_depend(facetpath)
            else:
                # run optimization of surface + reactants; surface + products
                # May need to put a post process on surface adsorbate
                # to call the next step
                # wait until optimization of surface + reactants; surface
                # + products finish and submit calculations to get TS guesses
                try:
                    self.run_opt_surf_and_adsorbate(facetpath)
                except NameError:
                    self.run_opt_surf_and_adsorbate_no_depend(facetpath)
                self.run_minima_vib('01', facetpath)
                self.run_ts_estimate('01', facetpath)
        # search for the 1st order saddle point
        self.run_ts_with_sella('02', facetpath)
        # run frequencies calculations for all TSs
        self.run_ts_vib('03', facetpath)
        # for each distinct TS, nudge towards imaginary frequency and
        # optimize to minima
        self.run_opt_after_ts('04', facetpath)

    def execute_all(self) -> None:
        ''' Main execute method for the entire workflow '''
        for facetpath in facetpaths:
            self.execute(facetpath)

    @ staticmethod
    def restart() -> None:
        LowLevelRestart().restart()
        HighLevelRestart().restart()
