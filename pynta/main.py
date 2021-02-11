#!/usr/bin/env python
from pynta.check_input import InputChecker
from pynta.restart import LowLevelRestart, HighLevelRestart
from pynta.io import IO
import json

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
check_input = os.path.join(working_dir, 'input.json')
check_run_me_py = os.path.join(working_dir, 'run_me.py')
check_run_me_sh = os.path.join(working_dir, 'run_me.sh')


# add working dir to system path
sys.path.insert(1, working_dir)
try:
    InputChecker(check_yaml, check_input, check_run_me_py,
                 check_run_me_sh, working_dir).check_all()
except ImportError:
    warn(
        'Missing input file. You cannot run calculations '
        'but will be able to use most of the workflow.'
    )

else:
    input_file = os.path.join(working_dir, 'input.json')
    with open(input_file, 'r') as f:
        input_json = json.load(f)

    quantum_chemistry = input_json['quantum_chemistry']
    calculator, socket_calculator = IO.get_calculators(quantum_chemistry)
    optimize_slab = input_json['optimize_slab']
    surface_types_and_repeats = input_json['surface_types_and_repeats']
    metal_atom = input_json['metal_atom']
    a = input_json['a']
    vacuum = input_json['vacuum']
    pseudo_dir = input_json['pseudo_dir']
    pseudopotentials = input_json['pseudopotentials']
    yamlfile = input_json['yamlfile']
    scfactor = input_json['scfactor']
    scfactor_surface = input_json['scfactor_surface']
    scaled1 = input_json['scaled1']
    scaled2 = input_json['scaled2']
    all_reacting_atoms = IO.get_all_reacting_atoms(check_yaml)
    species_dict = IO().get_species_dict(check_yaml)
    all_species = IO().get_all_unique_species_symbols(check_yaml)
    executable = input_json['executable']
    node_packing_count = input_json['node_packing_count']
    balsam_exe_settings = input_json['balsam_exe_settings']
    calc_keywords = input_json['calc_keywords']
    creation_dir = Path.cwd().as_posix()
    surface_types = surface_types_and_repeats.keys()
    facetpaths = IO().get_facetpaths(metal_atom, surface_types)
    job_file_dir_name = 'job_files'

# ###################################################
#                     Scripts                       #
# ###################################################

# These template and pytemplate scripts can be modified by users
# (if intended, though) to tune them to given calculation setup

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
        ''' Setup the balsam application.

        Once Quantum Chemistry calculations starts there will be one app for
        desired QE (e.g. `Quantum Espresso`) package and one for ``xTB``.

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

            IO.set_calculators(executable, calculator, socket_calculator)

        except SystemExit:
            print('---')
            print('Please create Balsam DB and/or activate it')
            print('---')

    @staticmethod
    def get_ts_xtb_py_scripts_list(
            facetpath: str) -> List[str]:
        '''Get a list with all 02 job scripts

        Parameters
        ----------
        facetpath : str
            a path to the workflow's main dir
            e.g. ``'Cu_111'``

        Returns
        -------
        ts_with_xtb_py_script_list : List[str]
            a list with all ``xTB`` jobs for a given facetpath

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

        Parameters
        ----------
        facetpath : str
            a path to the workflow's main dir
            e.g. ``'Cu_111'``

        Returns
        -------
        ts_sella_py_script_list : List[str]
            a list with all ``ts_sella_py`` jobs for a given facetpath

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

        Parameters
        ----------
        facetpath : str
            a path to the workflow's main dir
            e.g. ``'Cu_111'``

        Returns
        -------
        ts_vib_py_scripts_list : List[str]
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

        Parameters
        ----------
        facetpath : str
            a path to the workflow's main dir
            e.g. ``'Cu_111'``

        Returns:
        -------
        after_ts_py_scripts_list : List[str]
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
            socket_calculator: str,
            pytemplate: str,
            facetpath: str,
            slab_name: str,
            repeats: Tuple[int, int, int],
            creation_dir: PosixPath) -> None:
        ''' Create a python ASE job file for a big slab optimization

        Parameters
        ----------
        pytemplate : str
            a name of a template to prepare submission scripts
            for big_slab minimization
        facetpath : str
            a path to the workflow's main dir
            e.g.

            >>> facetpath = 'Cu_111'

        slab_name : str
            name of the slab file (no extintion) that is about to be created
            e.g.

            >>> slab_name = 'Cu_111_slab_opt'

        repeats : Tuple[int, int, int]
            how to replicate unit cell in (x, y, z) direction,
            e.g.

            >>> repeats = (3, 3, 1)

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
                    socket_calculator=socket_calculator,
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
        ''' For each surface type and for each reaction generate submit scripts
        for 6 stages of the workflow

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
                socket_calculator,
                template_slab_opt,
                py_job_dir,
                facetpath,
                surface_type,
                repeats_surface,
                slab_name
            )

            WorkFlow.set_up_big_slab(
                socket_calculator,
                template_big_slab_opt,
                py_job_dir,
                facetpath,
                slab_name,
                repeats,
                pytemplate_big_slab_opt
            )

            WorkFlow.set_up_ads(
                socket_calculator,
                template_ads,
                py_job_dir,
                facetpath,
                slab,
                repeats,
                yamlfile,
                pytemplate_relax_ads
            )

            WorkFlow.set_up_ads_vib(
                socket_calculator,
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
                    socket_calculator,
                    rxn,
                    template_set_up_ts,
                    py_job_dir,
                    facetpath,
                    slab,
                    repeats,
                    pytemplate_set_up_ts
                )

                WorkFlow.set_up_TS_vib(
                    socket_calculator,
                    rxn,
                    template_set_up_ts_vib,
                    py_job_dir,
                    facetpath,
                    slab,
                    repeats,
                    pytemplate_set_up_ts_vib
                )

                WorkFlow.set_up_opt_after_TS(
                    socket_calculator,
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
            socket_calculator: str,
            template: str,
            py_job_dir: str,
            facetpath: str,
            surface_type: str,
            repeats_surface: Tuple[int, int, int],
            slab_name: str) -> None:
        ''' Create ``'00_{facetpath}_set_up_slab_opt.py'`` file

        Parameters
        ----------
        template : str
            a name of a template for ``00`` job (slab optimization)
        py_job_dir : str
            a path to where all job files ``00``-``05`` are about to be
            created, e.g.

            >>> py_job_dir = {'current_dir'}/job_files/Cu_111

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        surface_type : str
            type of the surface. Available options are:
            fcc111, fcc211, fcc100, bcc111, bcc110, hcp0001, diamond111,
            diamond100
        metal_atom : str
            atomic metal_atom of the studied metal surface
            e.g. ``'Cu'``
        a : float
            a lattice constant
        repeats_surface : Tuple[int, int, int]
            surface multiplication in (x, y, z) direction
            eg.

            >>> repeats_surface = (1, 1, 4)

        vacuum : float
            amout of empty space in ``z`` direction (Angstrem)
        slab_name : str
            name of the slab file (no extintion) that is about to be created
            e.g.

            >>> slab_name = 'Cu_111_slab_opt'

        pseudopotentials : Dict[str,str]
            a dictionary with QE pseudopotentials for all species.
            e.g.

            >>> pseudopotentials = dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',
                                H='H.pbe-kjpaw_psl.1.0.0.UPF',
                                O='O.pbe-n-kjpaw_psl.1.0.0.UPF',
                                C='C.pbe-n-kjpaw_psl.1.0.0.UPF')

        pseudo_dir : str
            a path to the QE's pseudopotentials main directory
            e.g.

            >>> pseudo_dir = '/home/mgierad/espresso/pseudo'

        balsam_exe_settings : Dict[str,int]
            a dictionary with balsam execute parameters (cores, nodes, etc.),
            e.g.

            >>> balsam_exe_settings = {'num_nodes': 1,
                                    'ranks_per_node': 48,
                                    'threads_per_rank': 1}

        calc_keywords : Dict[str,str]
            a dictionary with parameters to run DFT package. Quantum Espresso
            is used as default, e.g.

            >>> calc_keywords = {'kpts': (3, 3, 1),
                                'occupations': 'smearing',
                                'smearing':  'marzari-vanderbilt',
                                'degauss': 0.01,
                                'ecutwfc': 40,
                                'nosym': True,
                                'conv_thr': 1e-11,
                                'mixing_mode': 'local-TF'}

        creation_dir : posix
            a posix path to the working directory

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            py_job_fname = os.path.join(
                py_job_dir, '00_{}_set_up_slab_opt.py'.format(facetpath))
            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    socket_calculator=socket_calculator,
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
            socket_calculator: str,
            template: str,
            py_job_dir: str,
            facetpath: str,
            slab_name: str,
            repeats: Tuple[int, int, int],
            pytemplate: str) -> None:
        ''' Set up a big_slab_opt ase py job generator

        Parameters
        ----------
        template : str
            a name of a template to set up ``01`` job
        py_job_dir : str
            a path to where all job files ``00``-``05`` are about to be
            created, e.g.

            >>> py_job_dir = {'current_dir'}/job_files/Cu_111

        facetpath : str
            a path to the workflow's main dir
            e.g.

            >>> facetpath = 'Cu_111'

        slab_name : str
            name of the slab file (no extintion) that is about to be created
            e.g.

            >>> slab_name = 'Cu_111_slab_opt'

        repeats : Tuple[int, int, int]
            how to replicate unit cell in (x, y, z) direction,
            e.g.

            >>> repeats = (3, 3, 1)

        pytemplate : str
            a template to prepare submission scripts for big_slab minimization

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            py_job_fname = os.path.join(
                py_job_dir, '00_{}_set_up_big_slab_opt.py'.format(facetpath))

            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    socket_calculator=socket_calculator,
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
            set_up_ads,
            template,
            py_job_dir,
            facetpath,
            slab,
            repeats,
            yamlfile,
            pytemplate):
        ''' Create ``'01_{facetpath}_set_up_ads_on_slab_{rxn_name}.py'`` file

        Parameters
        ----------

        template : str
            a name of a template to set up `01` job
        py_job_dir : str
            a path to where all job files ``00``-``05`` are about to be
            created, e.g.

            >>> py_job_dir = {'current_dir'}/job_files/Cu_111

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        slab : str
            a path to :literal:`*.xyz` file with optimized slab
        repeats : Tuple[int, int, int]
            how to replicate unit cell in (x, y, z) direction,
            e.g.

            >>> repeats = (3, 3, 1)

        yamlfile : str
            a name of the :literal:`*.yaml` file with a reaction list
        pytemplate : str
            a name of a template to prepare submission scripts for
            adsorbate+surface minimization

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            py_job_fname = os.path.join(
                py_job_dir, '01_{}_set_up_ads_on_slab.py'.format(facetpath))

            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    socket_calculator=socket_calculator,
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
            socket_calculator: str,
            template: str,
            py_job_dir: str,
            facetpath: str,
            pytemplate: str) -> None:
        ''' Create ``'{facetpath}_set_up_ads_vib.py'`` file

        Parameters
        ----------
        template : str
            a template file for ads_vib job
        py_job_dir : str
            a path to where all job files ``00``-``05`` are about to be
            created, e.g.

            >>> py_job_dir = {'current_dir'}/job_files/Cu_111

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        pytemplate : str
            a name of a template to prepare submission scripts for vibfrational
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
                    socket_calculator=socket_calculator,
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
        ''' Create ``'02_{facetpath}_set_up_TS_with_xtb_{rxn_name}.py'`` files

        Parameters
        ----------

        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction :literal:`*.yaml` file to a
            single reaction :literal:`*.yaml` file
        template : str
            a template to set up ``02`` job for a particular reaction
        py_job_dir : str
            a path to where all job files ``00``-``05`` are about to be
            created, e.g.

            >>> py_job_dir = {'current_dir'}/job_files/Cu_111

        slab : str
            a :literal:`*.xyz` file name with the optimized slab
            e.g.

            >>> slab = 'Cu_111_slab_opt.xyz'

        repeats : Tuple[int, int, int]
            how to replicate unit cell in (x, y, z) direction,
            e.g.

            >>> repeats = (3, 3, 1)

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        scfator : float
            a scaling factor to scale a bond distance between atoms taking
            part in the reaction
            e.g.

            >>> scfator = 1.4

        scfactor_surface : float
            a scaling factor to scale the target bond distance, i.e.
            the average distance between adsorbed atom and the nearest
            surface atom. Helpful e.g. when H is far away form the surface
            in TS, whereas for minima it is close to the surface
            e.g.

            >>> scfactor_surface = 1.0

        pytemplate_xtb : str
            a name of a template file for penalty function minimization job

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
                    node_packing_count=node_packing_count,
                    balsam_exe_settings=balsam_exe_settings,
                ))

    @staticmethod
    def set_up_run_TS(
            socket_calculator: str,
            rxn: Dict[str, str],
            template: str,
            py_job_dir: str,
            facetpath: str,
            slab: str,
            repeats: Tuple[int, int, int],
            pytemplate: str) -> None:
        ''' Create ``'03_{facetpath}_set_up_run_TS_{rxn_name}.py'`` file

        Parameters
        ----------

        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction :literal:`*.yaml` file to a
            single reaction :literal:`*.yaml` file
        template : str
            a template to set up 03 job for a particular reaction
        py_job_dir : str
            a path to where all job files ``00``-``05`` are about to be
            created, e.g.

            >>> py_job_dir = {'current_dir'}/job_files/Cu_111

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        slab : str
            a :literal:`*.xyz` file name with the optimized slab
            e.g.

            >>> slab = 'Cu_111_slab_opt.xyz'

        repeats : Tuple[int, int, int]
            how to replicate unit cell in (x, y, z) direction,
            e.g.

            >>> repeats = (3, 3, 1)

        pytemplate : str
            a name of a template file for ts optimization with sella

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            rxn_name = IO.get_rxn_name(rxn)
            py_job_fname = os.path.join(
                py_job_dir, '03_{}_run_TS_{}.py'.format(
                    facetpath, rxn_name))

            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    socket_calculator=socket_calculator,
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
            socket_calculator: str,
            rxn: Dict[str, str],
            template: str,
            py_job_dir: str,
            facetpath: str,
            slab: str,
            repeats: Tuple[int, int, int],
            pytemplate: str) -> None:
        ''' Create ``''04_{facetpath}_set_up_TS_vib_{rxn_name}.py'`` file

        Parameters
        ----------

        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction :literal:`*.yaml` file to a
            single reaction :literal:`*.yaml` file
        template : str
            a template to set up `04` job for a particular reaction
        py_job_dir : str
            a path to where all job files ``00``-``05`` are about to be
            created, e.g.

            >>> py_job_dir = {'current_dir'}/job_files/Cu_111

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        slab : str
            a :literal:`*.xyz` file name with the optimized slab
            e.g.
            >>> slab = 'Cu_111_slab_opt.xyz'

        repeats : Tuple[int, int, int]
            how to replicate unit cell in (x, y, z) direction,
            e.g.

            >>> repeats = (3, 3, 1)

        pytemplate : str
            a name of a template file for setting up frequencies calculations

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            rxn_name = IO.get_rxn_name(rxn)
            py_job_fname = os.path.join(
                py_job_dir, '04_{}_set_up_TS_vib_{}.py'.format(
                    facetpath, rxn_name))

            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    socket_calculator=socket_calculator,
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
            socket_calculator: str,
            rxn: Dict[str, str],
            template: str,
            py_job_dir: str,
            facetpath: str,
            slab: str,
            repeats: Tuple[int, int, int],
            pytemplate: str) -> None:
        ''' Create ``'05_{facetpath}_set_up_after_TS_{rxn_name}.py'`` file

        Parameters
        ----------

        rxn : Dict[str, str]
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction :literal:`*.yaml` file to a
            single reaction :literal:`*.yaml` file
        template : str
            a template to set up `05` job for a particular reaction
        py_job_dir : str
            a path to where all job files ``00``-``05`` are about to be
            created, e.g.

            >>> py_job_dir = {'current_dir'}/job_files/Cu_111

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        slab : str
            a :literal:`*.xyz` file name with the optimized slab
            e.g.

            >>> slab = 'Cu_111_slab_opt.xyz'

        repeats : Tuple[int, int, int]
            how to replicate unit cell in (x, y, z) direction,
            e.g.

            >>> repeats = (3, 3, 1)

        pytemplate : str
            a name of a template file for setting up an alternative to IRC
            (minimization of displaced structures following the imaginary mode
            of oscillation)

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            rxn_name = IO.get_rxn_name(rxn)
            py_job_fname = os.path.join(
                py_job_dir, '05_{}_set_up_after_TS_{}.py'.format(
                    facetpath, rxn_name))

            with open(py_job_fname, 'w') as c:
                c.write(template_text.format(
                    socket_calculator=socket_calculator,
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
        ''' Add python job file to balsam DB

        Parameters
        ----------

        parent_job : str
            a parent job on which subbmited jobs depends on. Formatted as:
            ``00`` or ``01`` or ``02`` or ``03`` or ``04`` or ``05``,
            depending on the job
        job_script : str
            a script that is about to be submitted
        cores : int
            number of cores for exe job
        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Returns
        -------

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

                pending_simulations = BalsamJob.objects.filter(
                    workflow__contains=dependency_workflow_name
                ).exclude(state='JOB_FINISHED')

                for job in pending_simulations:
                    add_dependency(job, job_to_add)  # parent, child

        return job_to_add

    def run_slab_optimization(
            self,
            facetpath: str) -> None:
        ''' Submit ``slab_optimization_job``

        Parameters
        ----------

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Returns
        -------

        job_to_add : BalsamJob
            job that will be submitted to balsam database

        '''
        slab_opt = '00_{}_set_up_slab_opt.py'.format(facetpath)
        self.slab_opt_job = self.exe('', slab_opt, facetpath, cores=1)

    def run_big_slab_opt(
            self,
            facetpath: str) -> None:
        ''' Submit ``big_slab_optimization`` job

        Parameters
        ----------

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Returns
        -------

        job_to_add : BalsamJob
            job that will be submitted to balsam database

        '''
        big_slab_opt = '00_{}_set_up_big_slab_opt.py'.format(facetpath)
        return self.exe('', big_slab_opt, facetpath, cores=1)

    def run_opt_surf_and_adsorbate(
            self,
            facetpath: str) -> None:
        ''' Run optmization of adsorbates on the surface

        Parameters
        ----------

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Returns
        -------

        job_to_add : BalsamJob
            job that will be submitted to balsam database

        '''
        ads_surf_opt_script = '01_{}_set_up_ads_on_slab.py'.format(facetpath)
        slab_opt = '00_{}_set_up_slab_opt.py'.format(facetpath)
        return self.exe(slab_opt, ads_surf_opt_script, facetpath)

    def run_opt_surf_and_adsorbate_no_depend(
            self,
            facetpath: str) -> None:
        ''' Run optmization of adsorbates on the surface if there is no
        dependency on other jobs

        Parameters
        ----------

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Returns
        -------

        job_to_add : BalsamJob
            job that will be submitted to balsam database

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
            e.g.
            >>> dependant_job = '01'
        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Returns
        -------

        job_to_add : BalsamJob
            job that will be submitted to balsam database

        '''
        ads_vib_script = '{}_set_up_ads_vib.py'.format(facetpath)
        return self.exe(dependent_job, ads_vib_script, facetpath)

    def run_minima_vib_no_depend(
            self,
            facetpath: str) -> None:
        ''' Run minima vibrational frequencies job

        Parameters
        ----------
        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Returns
        -------

        job_to_add : BalsamJob
            job that will be submitted to balsam database

        '''
        ads_vib_script = '{}_set_up_ads_vib.py'.format(facetpath)
        return self.exe('', ads_vib_script, facetpath)

    def run_ts_estimate(
            self,
            dependent_job: str,
            facetpath: str) -> None:
        ''' Run TS estimation calculations

        Parameters
        ----------

        dependant_job : str
            a prefix of the the job (minimum: two integers,
            max: whole job script name) that the current run depends on,
            e.g.

            >>> dependant_job = '01'

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        '''
        ts_xtb_py_script_list = WorkFlow.get_ts_xtb_py_scripts_list(facetpath)
        for ts_xtb in ts_xtb_py_script_list:
            self.exe(dependent_job, ts_xtb, facetpath)

    def run_ts_estimate_no_depend(
            self,
            facetpath: str) -> None:
        ''' Run TS estimate calculations if there is no dependency on the
        other jobs

        Parameters
        ----------

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        '''
        ts_xtb_py_script_list = WorkFlow.get_ts_xtb_py_scripts_list(facetpath)
        for ts_xtb in ts_xtb_py_script_list:
            self.exe('', ts_xtb, facetpath)

    def run_ts_with_sella(
            self,
            dependant_job: str,
            facetpath: str) -> None:
        ''' Run TS minimization with Sella

        Parameters
        ----------

        dependant_job : str
            a prefix of the the job (minimum: two integers,
            max: whole job script name) that the current depends on,
            e.g.

            >>> dependant_job = '02'

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

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

        Parameters
        ----------

        dependant_job : str
            a prefix of the the job (minimum: two integers,
            max: whole job script name) that the current depends on,
            e.g.

            >>> dependant_job = '03'

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

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

        Parameters
        ----------

        dependant_job : str
            a prefix of the the job (minimum: two integers,
            max: whole job script name) that the current depends on,
            e.g.

            >>> dependant_job = '04'

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        '''
        after_irc_py_scripts = WorkFlow.get_after_ts_py_scripts(facetpath)
        for after_irc in after_irc_py_scripts:
            self.exe(dependant_job, after_irc, facetpath)

    @staticmethod
    def check_all_species(
            yamlfile: str,
            facetpath: str) -> Dict[str, bool]:
        ''' Check all species(all reactions) to find whether there are
        previous calculation the code can use

        Parameters
        ----------

        yamlfile: str
            a name of the :literal:`*.yaml` file with a reaction list
        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Return
        ------

        checked_species: dict(str: bool)
            a dictionary with True/False values for every species(keys)
            True if there are previous calculation, otherwise False
            e.g. {'C': True, 'H': True, 'O': False, 'OH': True, 'CH': True}

        '''
        checked_species = {}
        all_species = IO().get_all_unique_species_symbols(yamlfile)
        for species in all_species:
            checked_species[species] = WorkFlow.is_minima_out_files(
                species, facetpath)
        return checked_species

    @staticmethod
    def is_minima_dir(
            species: str,
            facetpath: str) -> bool:
        ''' Check if directory with all minima calculations exists

        Parameters
        ----------

        species: str
            a species metal_atom
            e.g. ``'H'`` or ``'CO'``
        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Return
        ------
            ``True`` if directory for a given species exists in
            ``{facepath}/minima``. Otherwise, ``False``

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
        ''' Check for the previously calculated :literal:`*.relax.out` files
        for a given species.

        Parameters
        ----------

        species: str
            a species metal_atom
            e.g. ``'H'`` or ``'CO'``
        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Return
            ``True`` if there are previous calculations. Otherwise, ``False``.

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

        Parameters
        ----------

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Returns
        -------

        Tuple[bool, str=None]
            ``True`` if there are previous calculations

                >>> (True, path_to_prev_calc)

            `False` otherwise

                >>> (False, )

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
        ''' Check for ``big_slab`` calculations.


        Parameters
        ----------

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Return
        ------
        bool:
            ``True`` if there is a big_slab file, ``False`` otherwise.
            If multiple matches found, print all and raise error.

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

        Parameters
        ----------

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        Returns
        -------

        Tuple[bool, str=None]:
            ``True`` if there are previous calculations

                >>> (True, path_to_prev_calc)

            `False` otherwise

                >>> (False, )

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

        Parameters
        ----------

        facetpath : str
            a path to the workflow's main dir

            >>> facetpath = 'Cu_111'

        '''
        print('Starting calculations...')
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
                self.run_minima_vib_no_depend(facetpath)
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
        print('Running!')

    def execute_all(self) -> None:
        ''' Main execute method for the entire workflow '''
        for facetpath in facetpaths:
            self.execute(facetpath)

    @staticmethod
    def restart() -> None:
        LowLevelRestart().restart()
        HighLevelRestart().restart()
