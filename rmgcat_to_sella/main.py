#!/usr/bin/env python3
import os
import sys
import shutil
from pathlib import Path
from warnings import warn
from rmgcat_to_sella.io import IO
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
    slab_name = inputR2S.slab_name
    surface_types = inputR2S.surface_types
    symbol = inputR2S.symbol
    a = inputR2S.a
    vacuum = inputR2S.vacuum
    pseudo_dir = inputR2S.pseudo_dir
    pseudopotentials = inputR2S.pseudopotentials
    slabopt = inputR2S.slabopt
    yamlfile = inputR2S.yamlfile
    repeats = inputR2S.repeats
    repeats_surface = inputR2S.repeats_surface
    rotAngle = inputR2S.rotAngle
    scfactor = inputR2S.scfactor
    scfactor_surface = inputR2S.scfactor_surface
    scaled1 = inputR2S.scaled1
    scaled2 = inputR2S.scaled2
    species_dict = inputR2S.species_dict
    executable = inputR2S.executable
    balsam_exe_settings = inputR2S.balsam_exe_settings
    calc_keywords = inputR2S.calc_keywords
    creation_dir = inputR2S.creation_dir
    facetpaths = IO().get_facetpaths(symbol, surface_types)
    slab_names = [facetpath + '_slab_opt' for facetpath in facetpaths]
    slabopt = [slab_name + '.xyz' for slab_name in slab_names]

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
template_ads = os.path.join(path_template + '01_template_set_up_ads.py')
template_set_up_ts_with_xtb = os.path.join(
    path_template + '02_template_set_up_ts_with_xtb.py')
template_set_up_ts = os.path.join(
    path_template + '03_template_checksym_xtb_runTS.py')
template_set_up_ts_vib = os.path.join(
    path_template + '04_template_set_up_TS_vib.py')
template_set_up_after_ts = os.path.join(
    path_template + '05_template_set_up_after_ts.py')
pytemplate_relax_ads = os.path.join(
    path_pytemplate + 'pytemplate_relax_Cu_111_ads.py')
pytemplate_xtb = os.path.join(path_pytemplate + 'pytemplate_set_up_xtb.py')
pytemplate_set_up_ts = os.path.join(
    path_pytemplate + 'pytemplate_set_up_ts.py')
pytemplate_set_up_ts_vib = os.path.join(
    path_pytemplate + 'pytemplate_set_up_ts_vib.py')
pytemplate_set_up_after_ts = os.path.join(
    path_pytemplate + 'pytemplate_set_up_opt_after_ts.py')
slab_opt = '00_set_up_slab_opt.py'
ads_surf_opt_script = '01_set_up_ads.py'

####################################################
#                    Initialize                    #
####################################################


class WorkFlow:

    # def __init__(self):

        #     """Setup the balsam application for this workflow run.

        #     Once we start using QE will want one app for QE,
        #     one for xtb most likely
        #     """
        #     from balsam.core.models import ApplicationDefinition
        #     self.myPython, _ = ApplicationDefinition.objects.get_or_create(
        #         name="python",
        #         executable=sys.executable
        #     )
        #     self.myPython.save()
        #     self.slab_opt_job = ''

        #     # TODO: instead of directly importing EspressoBalsam, we should
        #     # write a function which returns the appropriate class from
        #     # balsamcalc.py based on the user-provided input file
        #     from rmgcat_to_sella.balsamcalc import (
        #         EspressoBalsam, EspressoBalsamSocketIO
        #     )
        #     EspressoBalsam.exe = executable
        #     EspressoBalsamSocketIO.exe = executable
        #     EspressoBalsam.create_application()
        #     EspressoBalsamSocketIO.create_application()

    def get_ts_xtb_py_script_list(self):
        ''' Get a list with all 02 job scripts '''
        reactions = IO().open_yaml_file(yamlfile)
        ts_with_xtb_py_script_list = []
        for rxn in reactions:
            rxn_name = IO().get_rxn_name(rxn)
            fname = '02_set_up_TS_with_xtb_{}.py'.format(rxn_name)
            ts_with_xtb_py_script_list.append(fname)
        return ts_with_xtb_py_script_list

    def get_ts_estimate_unique_list(self):
        ''' Get a list with all 03 job scripts '''
        reactions = IO().open_yaml_file(yamlfile)
        ts_sella_py_script_list = []
        for rxn in reactions:
            rxn_name = IO().get_rxn_name(rxn)
            fname = '03_checksym_xtb_run_TS_{}.py'.format(rxn_name)
            ts_sella_py_script_list.append(fname)
        return ts_sella_py_script_list

    def get_ts_vib_list(self):
        ''' Get a list with all 04 job scripts '''
        reactions = IO().open_yaml_file(yamlfile)
        ts_vib_py_scripts_list = []
        for rxn in reactions:
            rxn_name = IO().get_rxn_name(rxn)
            fname = '04_set_up_TS_vib_{}.py'.format(rxn_name)
            ts_vib_py_scripts_list.append(fname)
        return ts_vib_py_scripts_list

    def get_after_ts_py_scripts(self):
        ''' Get a list with all 05 job scripts '''
        reactions = IO().open_yaml_file(yamlfile)
        after_ts_py_scripts_list = []
        for rxn in reactions:
            rxn_name = IO().get_rxn_name(rxn)
            fname = '05_set_up_after_TS_{}.py'.format(rxn_name)
            after_ts_py_scripts_list.append(fname)
        return after_ts_py_scripts_list

    def gen_job_files(self):
        ''' Generate submit scripts for 6 stages of the workflow '''
        self.set_up_slab(
            template_slab_opt,
            surface_type,
            symbol,
            a,
            repeats_surface,
            vacuum,
            slab_name,
            pseudopotentials,
            pseudo_dir,
            balsam_exe_settings,
            calc_keywords,
            creation_dir
        )
        self.set_up_ads(
            template_ads,
            facetpath,
            slabopt,
            repeats,
            yamlfile,
            pytemplate_relax_ads,
            pseudopotentials,
            pseudo_dir,
            balsam_exe_settings,
            calc_keywords,
            creation_dir
        )

        reactions = IO().open_yaml_file(yamlfile)
        for rxn in reactions:
            self.set_up_TS_with_xtb(
                rxn,
                template_set_up_ts_with_xtb,
                slabopt,
                repeats,
                yamlfile,
                facetpath,
                rotAngle,
                scfactor,
                scfactor_surface,
                pytemplate_xtb,
                species_dict,
                creation_dir
            )

            self.set_up_run_TS(
                rxn,
                template_set_up_ts,
                facetpath,
                slabopt,
                repeats,
                yamlfile,
                pytemplate_set_up_ts,
                pseudopotentials,
                pseudo_dir,
                balsam_exe_settings,
                calc_keywords,
                creation_dir
            )

            self.set_up_TS_vib(
                rxn,
                template_set_up_ts_vib,
                facetpath,
                slabopt,
                repeats,
                yamlfile,
                pytemplate_set_up_ts_vib,
                pseudopotentials,
                pseudo_dir,
                balsam_exe_settings,
                calc_keywords,
                creation_dir
            )

            self.set_up_opt_after_TS(
                rxn,
                template_set_up_after_ts,
                facetpath,
                slabopt,
                repeats,
                yamlfile,
                pytemplate_set_up_after_ts,
                pseudopotentials,
                pseudo_dir,
                balsam_exe_settings,
                calc_keywords,
                creation_dir
            )

###########################
#   Create submit files   #
###########################

    def set_up_slab(
            self,
            template,
            surface_type,
            symbol,
            a,
            repeats_surface,
            vacuum,
            slab_name,
            pseudopotentials,
            pseudo_dir,
            balsam_exe_settings,
            calc_keywords,
            creation_dir):
        ''' Create 00_set_up_slab_opt.py file

        Parameters
        ----------
        template : py file
            a template for 00 job (slab optimization)
        surface_type : str
            type of the surface. Available options are:
            fcc111, fcc211, fcc100, bcc111, bcc110, hcp0001, diamond111,
            diamond100
        symbol : str
            atomic symbol of the studied metal surface
            e.g. 'Cu'
        a : float
            a lattice constant
        repeats_surface : tuple
            specify reapeats in (x, y, z) direction,
            eg. (3, 3, 1)
        vacuum : float
            amount of vacuum in the z direction (Units: Angstrem)
        slab_name : str
            a user defined name of the slab
        repeats_surface : tuple(int, int, int)
            surface multiplication in (x, y, z) direction
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
            with open('00_set_up_slab_opt.py', 'w') as c:
                c.write(template_text.format(
                    surface_type=surface_type,
                    symbol=symbol, a=a,
                    repeats_surface=repeats_surface,
                    vacuum=vacuum, slab_name=slab_name,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords, creation_dir=creation_dir
                ))

    def set_up_ads(
        self,
        template,
        facetpath,
        slabopt,
        repeats,
        yamlfile,
        pytemplate,
        pseudopotentials,
        pseudo_dir,
        balsam_exe_settings,
        calc_keywords,
        creation_dir
    ):
        ''' Create 01_set_up_ads.py file

        Parameters:
        ___________
        template : py file
            a template to set up 01 job
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        slabopt : str
            a path to .xyz file with optimized slab
        repeats : tuple(int, int, int)
            how to replicate unit cell in (x, y, z) direction
        yamlfile : str
            a name of the .yaml file with a reaction list
        pytemplate : python file
            a template to prepare submission scripts
            for adsorbate+surface minimization
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
            with open('01_set_up_ads.py', 'w') as c:
                c.write(template_text.format(
                    facetpath=facetpath, slabopt=slabopt,
                    yamlfile=yamlfile, repeats=repeats,
                    pytemplate=pytemplate,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=creation_dir
                ))
            c.close()
        r.close()

    def set_up_TS_with_xtb(
        self,
        rxn,
        template,
        slab,
        repeats,
        yamlfile,
        facetpath,
        rotAngle,
        scfactor,
        scfactor_surface,
        pytemplate_xtb,
        species_dict,
        creation_dir
    ):
        ''' Create 02_set_up_TS_with_xtb_{rxn_name}.py files

        Parameters:
        ___________

        rxn : dict(yaml[str:str])
            a dictionary with info about the paricular reaction. This can be
            view as a splitted many reaction .yaml file to a single reaction
            .yaml file
        template : py file
            a template to set up 02 job for a particular reaction
        slab : str
            a '.xyz' file name with the optimized slab
            e.g.
            'Cu_111_slab_opt.xyz'
        repeats : tuple(int, int, int)
            how to replicate unit cell in (x, y, z) direction
        yamlfile : str
            a name of the .yaml file with a reaction list
        facetpath : str
            a path to the workflow's main dir
            e.g. 'Cu_111'
        rotAngle : float
            an angle (deg) of rotation  of the TS guess adduct on the surface
            e.g. 60.0 (to be removed - not really necessary)
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
        species_dict : dict{str:list[str]}
            a dictionary holding info about particular reaction and key species
            for that reaction
            e.g. {'rxn1': ['O', 'H'], 'rxn2': ['C', 'H']}
        creation_dir : posix
            a posix path to the working directory

        '''
        with open(template, 'r') as r:
            template_text = r.read()
            rxn_name = IO().get_rxn_name(rxn)
            rxn_no = rxn['index'] + 1
            fname = '02_set_up_TS_with_xtb_{}.py'.format(rxn_name)
            with open(fname, 'w') as c:
                c.write(template_text.format(
                    facetpath=facetpath,
                    slab=slab,
                    repeats=repeats,
                    yamlfile=yamlfile,
                    rotAngle=rotAngle,
                    scfactor=scfactor,
                    scfactor_surface=scfactor_surface,
                    pytemplate_xtb=pytemplate_xtb,
                    species_list=species_dict['rxn' + str(rxn_no)],
                    scaled1=scaled1,
                    scaled2=scaled2,
                    creation_dir=creation_dir,
                    rxn=rxn,
                    rxn_name=rxn_name
                ))

    def set_up_run_TS(
        self,
        rxn,
        template,
        facetpath,
        slab,
        repeats,
        yamlfile,
        pytemplate,
        pseudopotentials,
        pseudo_dir,
        balsam_exe_settings,
        calc_keywords,
        creation_dir
    ):
        ''' Create 03_checksym_xtb_run_TS.py file '''
        with open(template, 'r') as r:
            template_text = r.read()
            rxn_name = IO().get_rxn_name(rxn)
            fname = '03_checksym_xtb_run_TS_{}.py'.format(rxn_name)
            with open(fname, 'w') as c:
                c.write(template_text.format(
                    facetpath=facetpath,
                    slab=slab,
                    repeats=repeats,
                    yamlfile=yamlfile,
                    pytemplate=pytemplate,
                    pseudo_dir=pseudo_dir,
                    pseudopotentials=pseudopotentials,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=creation_dir,
                    rxn=rxn,
                    rxn_name=rxn_name
                ))

    def set_up_TS_vib(
        self,
        rxn,
        template,
        facetpath,
        slab,
        repeats,
        yamlfile,
        pytemplate,
        pseudopotentials,
        pseudo_dir,
        balsam_exe_settings,
        calc_keywords,
        creation_dir
    ):
        ''' Create '04_set_up_TS_vib_{rxn_name}.py file '''
        with open(template, 'r') as r:
            template_text = r.read()
            rxn_name = IO().get_rxn_name(rxn)
            fname = '04_set_up_TS_vib_{}.py'.format(rxn_name)
            with open(fname, 'w') as c:
                c.write(template_text.format(
                    facetpath=facetpath,
                    slab=slab,
                    repeats=repeats,
                    yamlfile=yamlfile,
                    pytemplate=pytemplate,
                    pseudo_dir=pseudo_dir,
                    pseudopotentials=pseudopotentials,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=creation_dir,
                    rxn=rxn,
                    rxn_name=rxn_name
                ))

    def set_up_opt_after_TS(
        self,
        rxn,
        template,
        facetpath,
        slab,
        repeats,
        yamlfile,
        pytemplate,
        pseudopotentials,
        pseudo_dir,
        balsam_exe_settings,
        calc_keywords,
        creation_dir
    ):
        ''' Create '05_set_up_after_TS_{rxn_name}.py file '''
        with open(template, 'r') as r:
            template_text = r.read()
            rxn_name = IO().get_rxn_name(rxn)
            fname = '05_set_up_after_TS_{}.py'.format(rxn_name)
            with open(fname, 'w') as c:
                c.write(template_text.format(
                    facetpath=facetpath,
                    slab=slab,
                    repeats=repeats,
                    yamlfile=yamlfile,
                    pytemplate=pytemplate,
                    pseudo_dir=pseudo_dir,
                    pseudopotentials=pseudopotentials,
                    balsam_exe_settings=balsam_exe_settings,
                    calc_keywords=calc_keywords,
                    creation_dir=creation_dir,
                    rxn=rxn,
                    rxn_name=rxn_name
                ))

##############################
# Submit jobs and execute it #
##############################

    def exe(self,
            parent_job,
            job_script,
            cores=1):
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

        Returns:
        ________

        job_to_add : balsam job
            job that will be submitted to balsam queue/database

        '''
        from balsam.launcher.dag import BalsamJob
        from os import getcwd
        cwd = getcwd()

        # get rxn_name from job_script by spliting and joining job_script name
        # exeption for two first jobs
        if job_script in [ads_surf_opt_script, slab_opt]:
            rxn_name = ''
        else:
            rxn_name = '_'.join(job_script.split('_')[-2:])[:-3]
        try:
            workflow_name = facetpath + '_' + job_script[0:2] + '_' + rxn_name
        except ValueError:
            workflow_name = facetpath + '_error'

        job_to_add = BalsamJob(
            name=job_script,
            workflow=workflow_name,
            application=self.myPython.name,
            args=cwd + '/' + job_script,
            ranks_per_node=cores,
            input_files='',
            node_packing_count=48,
            user_workdir=cwd
        )
        job_to_add.save()

        # if there is a parent job, specify dependency
        if parent_job != '':
            from balsam.launcher.dag import add_dependency

            try:
                add_dependency(parent_job, job_to_add)  # parent, child
            except ValueError:
                dependency = str(parent_job[0:2])

                # a special case for 01 where there is on job script for all
                # reactions
                if parent_job == '01':
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

    def run_slab_optimization(self):
        ''' Submit slab_optimization_job '''
        self.slab_opt_job = self.exe('', slab_opt, cores=1)
        # submit slab_optimization_job 1 task probably, was 48 originally

    def run_opt_surf_and_adsorbate(self):
        ''' Run optmization of adsorbates on the surface '''
        return self.exe(self.slab_opt_job, ads_surf_opt_script)

    def run_opt_surf_and_adsorbate_no_depend(self):
        ''' Run optmization of adsorbates on the surface
            if there is no dependency on other jobs '''
        return self.exe('', ads_surf_opt_script)

    def run_ts_estimate(self, dependent_job):
        ''' Run TS estimation calculations '''
        ts_xtb_py_script_list = self.get_ts_xtb_py_script_list()
        for ts_xtb in ts_xtb_py_script_list:
            self.exe(dependent_job, ts_xtb)

    def run_ts_estimate_no_depend(self):
        ''' Run TS estimate calculations if there is
            no dependency on other jobs '''
        ts_xtb_py_script_list = self.get_ts_xtb_py_script_list()
        for ts_xtb in ts_xtb_py_script_list:
            self.exe('', ts_xtb)

    def run_ts_with_sella(self, dependant_job):
        ''' Run TS minimization with Sella '''
        ts_sella_py_script_list = self.get_ts_estimate_unique_list()
        for ts_sella in ts_sella_py_script_list:
            self.exe(dependant_job, ts_sella)

    def run_ts_vib(self, dependant_job):
        ''' Run frequency calculations for TS '''
        ts_vib_py_script_list = self.get_ts_vib_list()
        for ts_vib in ts_vib_py_script_list:
            self.exe(dependant_job, ts_vib)

    def run_opt_after_ts(self, dependant_job):
        ''' Run minimization of minima obtained as nudging TS structure
        towards imaginary mode of oscilation'''
        after_irc_py_scripts = self.get_after_ts_py_scripts()
        for after_irc in after_irc_py_scripts:
            self.exe(dependant_job, after_irc)

    def check_all_species(self, yamlfile):
        ''' Check all species (all reactions) to find whether
            there are previous calculation the code can use

        Parameters:
        ___________
        yamlfile : str
            a name of the .yaml file with a reaction list

        Return:
        _______
        checked_species : dict(str:bool)
            a dictionary with True/False values for every species (keys)
            True if there are previous calculation, otherwise False
            e.g. {'C': True, 'H': True, 'O': False, 'OH': True, 'CH': True}

        '''
        io = IO()
        checked_species = {}
        all_species = io.get_all_species(yamlfile)
        for species in all_species:
            # a bug to be resolved - why it inverts the name?
            if species == 'OH':
                species = 'HO'
            checked_species[species] = self.check_for_minima_outfiles(species)
        return checked_species

    def check_for_minima_dir(self, species):
        ''' Return True if directory for a given species exists in
            {facepath}/minima. Otherwise, False

        Parameters:
        ___________
        species : str
            a species symbol
            e.g. 'H' or 'CO'

        '''
        species_minima_dir = os.path.join(facetpath, 'minima', species)
        if os.path.isdir(species_minima_dir):
            return True
        return False

    def check_for_minima_outfiles(self, species):
        ''' Check for the previously calculated *relax.out files for a given
            species. Return True if there are previous calculations. Otherwise,
            False.

        Parameters:
        ___________
        species : str
            a species symbol
            e.g. 'H' or 'CO'

        '''
        minima_dir = os.path.join(facetpath, 'minima')
        # if minima dir exists, check for outfiles
        if self.check_for_minima_dir(species):
            search_for_outfiles = Path(minima_dir).glob(species + '_??_*out')
            outfiles = []
            for outfile in search_for_outfiles:
                outfiles.append(str(outfile))
            # empty list
            if not outfiles:
                return False
            return True
        return False

    def check_if_slab_opt_exists(self):
        ''' Check whether slab has been already optimized

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

    def copy_slab_opt_file(self):
        ''' Copy .xyz of previously optimized slab '''
        self.slab_exists = self.check_if_slab_opt_exists()
        if self.slab_exists[0]:
            src = self.slab_exists[1]
            dst = os.getcwd()
            try:
                shutil.copy2(src, dst)
                self.slab_opt_job = ''
            except shutil.SameFileError:
                pass

    def execute(self):
        ''' The main executable '''

        if optimize_slab:
            # if slab found in previous calculation, do nothing
            if self.check_if_slab_opt_exists()[0]:
                pass
                # self.copy_slab_opt_file()
            else:
                # If the code cannot locate optimized slab .xyz file,
                # a slab optimization will be launched.
                self.run_slab_optimization()
            # check if species were already calculated
            if all(self.check_all_species(yamlfile).values()):
                # If all are True, start by generating TS guesses and run
                # the penalty function minimization
                self.run_ts_estimate_no_depend()
            else:
                # If any of sp_check_list is False
                # run optimization of surface + reactants; surface + products
                try:
                    self.run_opt_surf_and_adsorbate()
                except NameError:
                    self.run_opt_surf_and_adsorbate_no_depend()
                self.run_ts_estimate('01')
        else:
            # this is executed if user provide .xyz with the optimized slab
            # and explicitly define oiptimize_slab = False
            if self.check_if_slab_opt_exists()[0]:
                pass
            else:
                raise FileNotFoundError(
                    'It appears that there is no slab_opt.xyz file'
                )
            if all(self.check_all_species(yamlfile).values()):
                # If all minima were calculated some time age rmgcat_to_sella
                # will use that calculations. Start from 02 step
                self.run_ts_estimate_no_depend()
            else:
                # run optimization of surface + reactants; surface + products
                # May need to put a post process on surface adsorbate
                # to call the next step
                # wait until optimization of surface + reactants; surface
                # + products finish and submit calculations to get TS guesses
                try:
                    self.run_opt_surf_and_adsorbate()
                except NameError:
                    self.run_opt_surf_and_adsorbate_no_depend()
                self.run_ts_estimate('01')
        # search for the 1st order saddle point
        self.run_ts_with_sella('02')
        # run frequencies calculations for all TSs
        self.run_ts_vib('03')
        # for each distinct TS, nudge towards imaginary frequency and
        # optimize to minima
        self.run_opt_after_ts('04')
