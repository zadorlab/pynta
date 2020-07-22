import os
import time
import shutil
from pathlib import Path
try:
    import inputR2S
    '''
    User defined parameters. Here we only read them. They are set up in inputR2S.py (submit directory)
    '''
    facetpath = inputR2S.facetpath
    slab_name = inputR2S.slab_name
    surface_type = inputR2S.surface_type
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
    sp1 = inputR2S.sp1
    sp2 = inputR2S.sp2
    scaled1 = inputR2S.scaled1
    scaled2 = inputR2S.scaled2
    species_list = [sp1, sp2]
    slab_opt = inputR2S.slab_opt_script
    SurfaceAdsorbate = inputR2S.SurfaceAdsorbateScript
    TSxtb = inputR2S.TSxtbScript
    TS = inputR2S.TSScript
    IRC = inputR2S.IRCScript
    IRCopt = inputR2S.IRCoptScript
#    from pathlib import Path
#    creation_dir = Path.cwd().as_posix() 

except ImportError:
    print('Missing input file. You cannot run calculations but will be able to use most of the workflow.')

# These template and pytemplate scripts can be modified by users to tune
# them to given calculation setup, i.e. calculator, method, queue menager,
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
template_set_up_IRC = os.path.join(path_template + '04_template_set_up_irc.py')
template_set_up_optIRC = os.path.join(
    path_template + '05_template_set_up_opt_after_irc.py')
pytemplate_relax_ads = os.path.join(
    path_pytemplate + 'pytemplate_relax_Cu_111_ads.py')
pytemplate_xtb = os.path.join(path_pytemplate + 'pytemplate_set_up_xtb.py')
pytemplate_set_up_ts = os.path.join(
    path_pytemplate + 'pytemplate_set_up_ts.py')
pytemplate_f = os.path.join(path_pytemplate + 'pytemplate_set_up_irc_f.py')
pytemplate_r = os.path.join(path_pytemplate + 'pytemplate_set_up_irc_r.py')
pytemplate_optIRC = os.path.join(
    path_pytemplate + 'pytemplate_set_up_opt_irc.py')

slab_opt = inputR2S.slab_opt_script
SurfaceAdsorbate = inputR2S.SurfaceAdsorbateScript
TSxtb = inputR2S.TSxtbScript
TS = inputR2S.TSScript
IRC = inputR2S.IRCScript
IRCopt = inputR2S.IRCoptScript
##
currentDir = os.path.dirname(os.getcwd())
# sp1 = inputR2S.sp1
# sp2 = inputR2S.sp2
# facetpath = inputR2S.facetpath
optimize_slab = inputR2S.optimize_slab

####################################################
#################### Initialize ####################
####################################################


class WorkFlow:

    def __init__(self):
        ''' Setup the balsam application for this workflow run, once we start using QE will want one app for QE, one for xtb most likely '''
        from balsam.core.models import ApplicationDefinition
        self.myPython, self.app_created = ApplicationDefinition.objects.get_or_create(name="Python", executable="python")
        # envscript="/path/to/setup-envs.sh",
        #postprocess="python /path/to/post.py"
    # def __init__(self, facetpath):
    #     self.facetpath = facetpath

    def gen_job_files(self):
        ''' Generate submt scripts for 6 stages of the workflow '''
        self.set_up_slab(template_slab_opt, surface_type, symbol, a,
                             repeats_surface, vacuum, slab_name,
                             pseudopotentials, pseudo_dir)
        self.set_up_ads(template_ads, facetpath, slabopt,
                            repeats, yamlfile, pytemplate_relax_ads,
                            pseudopotentials, pseudo_dir)
        self.set_up_TS_with_xtb(template_set_up_ts_with_xtb, slabopt,
                                    repeats, yamlfile, facetpath, rotAngle,
                                    scfactor, scfactor_surface, pytemplate_xtb,
                                    sp1, sp2)
        self.set_up_run_TS(template_set_up_ts, facetpath, slabopt,
                               repeats, yamlfile, pytemplate_set_up_ts,
                               pseudopotentials, pseudo_dir)
        self.set_up_run_IRC(template_set_up_IRC, facetpath, slabopt,
                                repeats, pytemplate_f, pytemplate_r, yamlfile,
                                pseudopotentials, pseudo_dir)
        self.set_up_opt_IRC(template_set_up_optIRC,
                                facetpath, slabopt, repeats,
                                pytemplate_optIRC,
                                pseudopotentials, pseudo_dir)

###########################
#   Create submit files   #
###########################

    def set_up_slab(self, template, surface_type, symbol, a, repeats_surface,
                    vacuum, slab_name, pseudopotentials, pseudo_dir):
        ''' Create 00_set_up_slab_opt.py file '''
        with open(template, 'r') as r:
            template_text = r.read()
            with open('00_set_up_slab_opt.py', 'w') as c:
                c.write(template_text.format(surface_type=surface_type,
                                        symbol=symbol, a=a,
                                        repeats_surface=repeats_surface,
                                        vacuum=vacuum, slab_name=slab_name,
                                        pseudopotentials=pseudopotentials,
                                        pseudo_dir=pseudo_dir))
            c.close()
        r.close()

    def set_up_ads(self, template, facetpath, slabopt, repeats, yamlfile,
                   pytemplate, pseudopotentials, pseudo_dir):
        ''' Create 01_set_up_ads.py file '''
        with open(template, 'r') as r:
            template_text = r.read()
            with open('01_set_up_ads.py', 'w') as c:
                c.write(template_text.format(facetpath=facetpath, slabopt=slabopt,
                    yamlfile=yamlfile, repeats=repeats,
                    pytemplate=pytemplate,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir))
            c.close()
        r.close()

    def set_up_TS_with_xtb(self, template, slab,
                           repeats, yamlfile, facetpath, rotAngle,
                           scfactor, scfactor_surface,
                           pytemplate_xtb, sp1, sp2):
        ''' Create 02_set_up_TS_with_xtb.py file'''
        with open(template, 'r') as r:
            template_text = r.read()
            with open('02_set_up_TS_with_xtb.py', 'w') as c:
                c.write(template_text.format(facetpath=facetpath, slab=slab,
                                        repeats=repeats, yamlfile=yamlfile,
                                        rotAngle=rotAngle, scfactor=scfactor,
                                        scfactor_surface=scfactor_surface,
                                        pytemplate_xtb=pytemplate_xtb, sp1=sp1,
                                        sp2=sp2, scaled1=scaled1,
                                        scaled2=scaled2))
            c.close()
        r.close()

    def set_up_run_TS(self, template, facetpath, slab, repeats, yamlfile,
                      pytemplate, pseudopotentials, pseudo_dir):
        ''' Create 03_checksym_xtb_runTS.py file '''
        with open(template, 'r') as r:
            template_text = r.read()
            with open('03_checksym_xtb_runTS.py', 'w') as c:
                c.write(template_text.format(facetpath=facetpath, slab=slab,
                                        repeats=repeats, yamlfile=yamlfile,
                                        pytemplate=pytemplate,
                                        pseudo_dir=pseudo_dir,
                                        pseudopotentials=pseudopotentials))
            c.close()
        r.close()

    def set_up_run_IRC(self, template, facetpath, slab, repeats,
                       pytemplate_f, pytemplate_r, yamlfile,
                       pseudopotentials, pseudo_dir):
        ''' Create 04_set_up_irc.py file '''
        with open(template, 'r') as r:
            template_text = r.read()
            with open('04_set_up_irc.py', 'w') as c:
                c.write(template_text.format(facetpath=facetpath,
                                        slab=slab,
                                        repeats=repeats,
                                        pytemplate_f=pytemplate_f,
                                        pytemplate_r=pytemplate_r,
                                        yamlfile=yamlfile,
                                        pseudo_dir=pseudo_dir,
                                        pseudopotentials=pseudopotentials))
            c.close()
        r.close()

    def set_up_opt_IRC(self, template, facetpath, slab, repeats, pytemplate,
                       pseudopotentials, pseudo_dir):
        ''' Create 05_set_up_opt_after_irc.py file'''
        with open(template, 'r') as r:
            template_text = r.read()
            with open('05_set_up_opt_after_irc.py', 'w') as c:
                c.write(template_text.format(facetpath=facetpath,
                                        slab=slab,
                                        repeats=repeats,
                                        pytemplate=pytemplate,
                                        yamlfile=yamlfile,
                                        pseudo_dir=pseudo_dir,
                                        pseudopotentials=pseudopotentials))
            c.close()
        r.close()

##############################
# Submit jobs and execute it #
##############################
    def exe(self, parent_job, job_script,cores=1):
        from balsam.launcher.dag import BalsamJob
        from os import getcwd
        cwd = getcwd()
        job_to_add = BalsamJob(
                name = job_script,
                workflow = "test",
                application = self.myPython,
                args = cwd+'/'+job_script,
                ranks_per_node = cores,
                user_workdir = cwd
                #working_directory = cwd
                )
        job_to_add.save()
        if parent_job!='':
            from balsam.launcher.dag import add_dependency
            try:
                add_dependency(parent_job,job_to_add) # parent, child
            except ValueError:
                dependency=str(int(parent_job[0:1]))
                dependency_workflow_name = yamlfile+facetpath+dependency
                BalsamJob = BalsamJob
                pending_simulations = BalsamJob.objects.filter(workflow__contains=dependency_workflow_name).exclude(state='JOB_FINISHED')
                for job in pending_simulations:
                    add_dependency(job,job_to_add) # parent, child
        return job_to_add

    def check_if_path_to_mimina_exists(self, WorkFlowDir, species):
        ''' Check for the paths to previously calculated minima and return 
            a list with all valid paths '''

        pathlist = Path(WorkFlowDir).glob('**/minima/' + species)
        paths = []
        for path in pathlist:
            # path = str(path)
            paths.append(str(path))
            return paths[0]
        if IndexError:
            return None

    def check_if_minima_already_calculated(self, currentDir, species,
                                           facetpath):
        ''' Check for previously calculated minima '''
        WorkFlowDirs = []
        uniqueMinimaDirs = []
        if facetpath == 'Cu_211':
            WorkFlowDirsList = Path(str(currentDir)).glob('*_Cu_211_methanol*')
        elif facetpath == 'Cu_100':
            WorkFlowDirsList = Path(str(currentDir)).glob('*_Cu_100_methanol*')
        elif facetpath == 'Cu_111':
            WorkFlowDirsList = Path(str(currentDir)).glob('*_Cu_methanol*')
        else:
            WorkFlowDirsList = Path(str(currentDir)).glob('*_Cu_methanol*')
        # transforming posix path to regular string
        for WorkFlowDir in WorkFlowDirsList:
            WorkFlowDirs.append(WorkFlowDir)
        # expected e.g. -> WorkFlowDirs = ['00_Cu_methanol_CO+O_CO2',
        # '01_Cu_methanol_OH_O+H', '02_Cu_methanol_CO+H_HCO']
        for WorkFlowDir in WorkFlowDirs:
            minimaDir = WorkFlow.check_if_path_to_mimina_exists(WorkFlowDir, species)
            if minimaDir is not None:
                uniqueMinimaDirs.append(minimaDir)
        if len(uniqueMinimaDirs) >= 1:
            print('More than one possible path were found for the species {}. Choosing the following path: {}'.format(
                species, uniqueMinimaDirs[0]))
            return True, uniqueMinimaDirs[0]
        elif IndexError:
            print('Species {} was not yet calculated. Setting up new calculations.'.format(
                species))
            return (False, )

    def run_slab_optimization(self):
        ''' Submit slab_optimization_job '''
        self.slab_opt_job=self.exe('', slab_opt,cores=48)
        # submit slab_optimization_job

    def run_opt_surf_and_adsorbate(self):
        ''' Run optmization of adsorbates on the surface '''
        return self.exe(self.slab_opt_job, SurfaceAdsorbate)

    def run_opt_surf_and_adsorbate_no_depend(self):
        ''' Run optmization of adsorbates on the surface
            if there is no dependency on other jobs '''
        return self.exe('', SurfaceAdsorbate)

    def run_ts_estimate(self, dependent_job):
        ''' Run TS estimation calculations '''
        TSxtb = inputR2S.TSxtbScript
        return self.exe(dependent_job, TSxtb)

    def run_ts_estimate_no_depend(self):
        ''' Run TS estimate calculations if there is
            no dependency on other jobs '''
        TSxtb = inputR2S.TSxtbScript
        return self.exe('', TSxtb)

    def check_all_species(self):
        ''' Check all species to find whether there are previous calculation
            the code can use

        Return:
        _______
        all_sp_checked : list(tuple(bool, str=None))
            a list of tuples with info whether a species were calculated
            (True, path_to_prev_calc)
            or not
            (False, )
        '''
        all_sp_checked = []
        for species in species_list:
            check_sp = self.check_if_minima_already_calculated(currentDir, species, facetpath)
            all_sp_checked.append(check_sp)
        return all_sp_checked

    def check_if_slab_opt_exists(self):
        ''' Check whether slab has been already optimized

        Returns : tuple(bool, str=None):
            True if there are previous calculations
                (True, path_to_prev_calc)
            False otherwise
                (False, )

        '''
        slab_opt_path_str = []
        # the code will look for anything like Cu_111*.xyz starting from the
        # facetpath directory including all subdirectories.
        keyphrase = '**/*' + str(facetpath) + '*.xyz'
        slab_opt_path_posix = Path(str(currentDir)).glob(keyphrase)
        for slab_opt_path in slab_opt_path_posix:
            slab_opt_path_str.append(slab_opt_path)
        if len(slab_opt_path_str) >= 1:
            return True, slab_opt_path_str[0]
        else:
            return (False, )

    def copy_slab_opt_file(self):
        ''' Copy .xyz of previously optimized slab '''
        self.slab_exists= self.check_if_slab_opt_exists()
        print(self.slab_exists)
        if self.slab_exists[0]:
            src = self.slab_exists[1]
            dst = os.getcwd()
            shutil.copy2(src, dst)
            self.slab_opt_job=''

    def execute(self):
        ''' The main executable '''
        checksp1, checksp2 = self.check_all_species()

        if optimize_slab is True:
            # If the code cannot locate optimized slab .xyz file,
            # a slab optimization will be launched.
            # a = self.check_if_slab_opt_exists()
            # print(a)
            if self.check_if_slab_opt_exists():
                self.run_slab_optimization()
                # wait a bit in case the file write process is too slow
            else:
                self.copy_slab_opt_file()
            # check whether sp1 and sp2 was already cacluated
            if checksp1[0] is False or checksp2[0] is False:
                # If any of these is False
                # run optimization of surface + reactants; surface + products
                try:
                    self.run_opt_surf_and_adsorbate()
                except NameError:
                    self.run_opt_surf_and_adsorbate_no_depend()
                # run calculations to get TS guesses
                self.run_ts_estimate('01')
            else:
                # If both are True, start by generating TS guesses and run
                # the penalty function minimization
                self.run_ts_estimate_no_depend()
                # self.run_ts_estimate('00')
        else:
            # this is executed if user provide .xyz with the optimized slab
            # check whether sp1 and sp2 was already cacluated
            if self.check_if_slab_opt_exists():
                pass
            else:
                raise FileNotFoundError(
                    'It appears that there is no slab_opt.xyz file')

            if checksp1[0] is False or checksp2[0] is False:
                # run optimization of surface + reactants; surface + products
                self.exe('',SurfaceAdsorbate)
                # wait a bit in case the file write process is too slow
                """ May need to put a post process on surface adsorbate to call the next step """
                # wait until optimization of surface + reactants; surface + products
                # finish and submit calculations to get TS guesses
                self.exe('01', TSxtb)
            else:
                # If all minimas were calculated some time age for other reaction,
                # rmgcat_to_sella will use that calculations.
                self.exe('', TSxtb)
        # search for the 1st order saddle point
        self.exe('02', TS)
        # for each distinct TS, run IRC calculations
        self.exe('03', IRC)
        # run optimizataion of both IRC (forward, reverse) trajectory
        self.exe('04', IRCopt)
