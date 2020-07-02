import os
import time
from pathlib import Path
import sys
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
sp1 = inputR2S.sp1
sp2 = inputR2S.sp2
facetpath = inputR2S.facetpath
optimize_slab = inputR2S.optimize_slab

####################################################
#################### Initialize ####################
####################################################

# TODO: It would be great to hava a class here


class WorkFlow:
    # def __init__(self,):
    #     self.surface_type = surface_type

    def genJobFiles(self):
        ''' Generate submt scripts for 6 stages of the workflow '''
        WorkFlow.set_up_slab(self, template_slab_opt, surface_type, symbol, a,
                             repeats_surface, vacuum, slab_name,
                             pseudopotentials)
        WorkFlow.set_up_ads(self, template_ads, facetpath, slabopt, yamlfile,
                            repeats, pytemplate_relax_ads)
        WorkFlow.set_up_TS_with_xtb(self, template_set_up_ts_with_xtb, slabopt,
                                    repeats, yamlfile, facetpath, rotAngle,
                                    scfactor, scfactor_surface, pytemplate_xtb,
                                    sp1, sp2)
        WorkFlow.set_up_run_TS(self, template_set_up_ts,
                               facetpath, pytemplate_set_up_ts)
        WorkFlow.set_up_run_IRC(self, template_set_up_IRC, facetpath,
                                pytemplate_f, pytemplate_r, yamlfile)
        WorkFlow.set_up_opt_IRC(self, template_set_up_optIRC,
                                facetpath, pytemplate_optIRC)

###########################
#   Create submit files   #
###########################

    def set_up_slab(self, template, surface_type, symbol, a, repeats_surface,
                    vacuum, slab_name, pseudopotentials):
        ''' Create 00_set_up_slab_opt.py file '''
        with open(template, 'r') as r:
            template = r.read()
            with open('00_set_up_slab_opt.py', 'w') as c:
                c.write(template.format(surface_type=surface_type,
                                        symbol=symbol, a=a,
                                        repeats_surface=repeats_surface,
                                        vacuum=vacuum, slab_name=slab_name,
                                        pseudopotentials=pseudopotentials))
            c.close()
        r.close()

    def set_up_ads(self, template, facetpath, slabopt, yamlfile, repeats,
                   pytemplate):
        ''' Create 01_set_up_ads.py file '''
        with open(template, 'r') as r:
            template = r.read()
            with open('01_set_up_ads.py', 'w') as c:
                c.write(template.format(facetpath=facetpath, slabopt=slabopt,
                                        yamlfile=yamlfile, repeats=repeats,
                                        pytemplate=pytemplate))
            c.close()
        r.close()

    def set_up_TS_with_xtb(self, template, slab,
                           repeats, yamlfile, facetpath, rotAngle,
                           scfactor, scfactor_surface,
                           pytemplate_xtb, sp1, sp2):
        ''' Create 02_set_up_TS_with_xtb.py file'''
        with open(template, 'r') as r:
            template = r.read()
            with open('02_set_up_TS_with_xtb.py', 'w') as c:
                c.write(template.format(facetpath=facetpath, slab=slab,
                                        repeats=repeats, yamlfile=yamlfile,
                                        rotAngle=rotAngle, scfactor=scfactor,
                                        scfactor_surface=scfactor_surface,
                                        pytemplate_xtb=pytemplate_xtb, sp1=sp1,
                                        sp2=sp2, scaled1=scaled1,
                                        scaled2=scaled2))
            c.close()
        r.close()

    def set_up_run_TS(self, template, facetpath, pytemplate):
        ''' Create 03_checksym_xtb_runTS.py file '''
        with open(template, 'r') as r:
            template = r.read()
            with open('03_checksym_xtb_runTS.py', 'w') as c:
                c.write(template.format(facetpath=facetpath,
                                        pytemplate=pytemplate))
            c.close()
        r.close()

    def set_up_run_IRC(self, template, facetpath,
                       pytemplate_f, pytemplate_r, yamlfile):
        ''' Create 04_set_up_irc.py file '''
        with open(template, 'r') as r:
            template = r.read()
            with open('04_set_up_irc.py', 'w') as c:
                c.write(template.format(facetpath=facetpath,
                                        pytemplate_f=pytemplate_f,
                                        pytemplate_r=pytemplate_r,
                                        yamlfile=yamlfile))
            c.close()
        r.close()

    def set_up_opt_IRC(self, template, facetpath, pytemplate):
        ''' Create 05_set_up_opt_after_irc.py file'''
        with open(template, 'r') as r:
            template = r.read()
            with open('05_set_up_opt_after_irc.py', 'w') as c:
                c.write(template.format(facetpath=facetpath,
                                        pytemplate=pytemplate,
                                        yamlfile=yamlfile))
            c.close()
        r.close()

##############################
# Submit jobs and execute it #
##############################

    def get_slurm_jobs_id(self, slurm_id_subm):
        ''' Get slurm IDs of just submitted jobs '''
        slurm_jobs_id = []
        with open(slurm_id_subm, 'r') as f:
            for line in f.readlines():
                line = line.split()[3]
                slurm_jobs_id.append(line)
        f.close()
        print(slurm_jobs_id)
        return slurm_jobs_id

    def gen_slurm_command(self, slurm_id_subm):
        ''' Prepare a bash command to submit jobs '''
        slurmID = WorkFlow.get_slurm_jobs_id(self, slurm_id_subm)
        slurmID = ",".join(["{}"] * len(slurmID)).format(*slurmID)
        # if not slurmID:
        #     print('Error')
        #     sys.exit('No submitted jobs, probably all files have been already generated')
        #     return command = False
        # else:
        command = os.path.join('sbatch --dependency=afterany:' + str(slurmID))
        print(command)
        return command

    def run(self, slurm_id_subm, job_script):
        ''' Submit slurm jobs '''
        command = WorkFlow.gen_slurm_command(self, slurm_id_subm)
        os.popen(str(os.path.join(command + ' ' + job_script)))

    def exe(self, prevSlurmID, job_script):
        ''' Check if the previous step of calculations terminated. If so, run the next step'''
        while not os.path.exists(prevSlurmID):
            time.sleep(60)
        WorkFlow.run(self, prevSlurmID, job_script)

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
        else:
            WorkFlowDirsList = Path(str(currentDir)).glob('*_Cu_methanol*')
        # transforming posix path to regular string
        for WorkFlowDir in WorkFlowDirsList:
            WorkFlowDirs.append(WorkFlowDir)
        # expected e.g. -> WorkFlowDirs = ['00_Cu_methanol_CO+O_CO2',
        # '01_Cu_methanol_OH_O+H', '02_Cu_methanol_CO+H_HCO']
        for WorkFlowDir in WorkFlowDirs:
            minimaDir = WorkFlow.check_if_path_to_mimina_exists(
                self, WorkFlowDir, species)
            if minimaDir is not None:
                uniqueMinimaDirs.append(minimaDir)

        if len(uniqueMinimaDirs) >= 1:
            print('More than one possible path were found for the species {}. Choosing the following path: {}'.format(
                species, uniqueMinimaDirs[0]))
            # print(uniqueMinimaDirs[0])
            return True, uniqueMinimaDirs[0]
        elif IndexError:
            print('Species {} was not yet calculated. Setting up new calculations.'.format(
                species))
            return False

    def run_slab_optimization(self):
        ''' Submit slab_optimization_job '''
        # submit slab_optimization_job
        run_slab_command = os.path.join(
            'sbatch ' + slab_opt + ' > submitted_00.txt')
        run_slab_opt = os.popen(run_slab_command)
        print(run_slab_opt.read())

    def run_opt_surf_and_adsorbate(self):
        return WorkFlow.exe(self, 'submitted_00.txt', SurfaceAdsorbate)

    def run_ts_estimate(self, submit_txt):
        return WorkFlow.exe(self, submit_txt, TSxtb)

    # def check_all_species(self):
    #     checksp1 = WorkFlow.check_if_minima_already_calculated(
    #             self, currentDir, sp1, facetpath)
    #     checksp2 = WorkFlow.check_if_minima_already_calculated(
    #             self, currentDir, sp2, facetpath)

    def execute(self):
        ''' The main executable '''
        if optimize_slab is True:
            WorkFlow.run_slab_optimization(self)
            # wait a bit in case the file write process is too slow
            while not os.path.exists('submitted_00.txt'):
                time.sleep(3)
                # check whether sp1 and sp2 was already cacluated
            checksp1 = WorkFlow.check_if_minima_already_calculated(
                self, currentDir, sp1, facetpath)
            checksp2 = WorkFlow.check_if_minima_already_calculated(
                self, currentDir, sp2, facetpath)
            if checksp1 is False and checksp2 is False:
                # run optimization of surface + reactants; surface + products
                WorkFlow.run_opt_surf_and_adsorbate(self)
                # run calculations to get TS guesses
                WorkFlow.run_ts_estimate(self, 'submitted_01.txt')
            else:
                # If not, start by generating TS guesses and use
                # penalty function minimization
                WorkFlow.run_ts_estimate(self, 'submitted_00.txt')
        else:
            # this is executed if user provide .xyz with the optimized slab
            # check whether sp1 and sp2 was already cacluated
            checksp1 = WorkFlow.check_if_minima_already_calculated(
                self, currentDir, sp1, facetpath)
            checksp2 = WorkFlow.check_if_minima_already_calculated(
                self, currentDir, sp2, facetpath)
            if checksp1 is False and checksp2 is False:
                # run optimization of surface + reactants; surface + products
                runSurfAds = os.popen(os.path.join(
                    'sbatch ' + SurfaceAdsorbate))
                print(runSurfAds.read())
                # wait a bit in case the file write process is too slow
                while not os.path.exists('submitted_01.txt'):
                    time.sleep(1)
                # wait until optimization of surface + reactants; surface + products
                # finish and submit calculations to get TS guesses
                WorkFlow.run(self, 'submitted_01.txt', TSxtb)
            else:
                # If all minimas were calculated some time age for other reaction,
                # rmgcat_to_sella will use that calculations.
                run_TSxtb = os.popen(os.path.join(
                    'sbatch ' + TSxtb))
                print(run_TSxtb.read())
        # search for the 1st order saddle point
        WorkFlow.exe(self, 'submitted_02.txt', TS)
        # for each distinct TS, run IRC calculations
        WorkFlow.exe(self, 'submitted_03.txt', IRC)
        # run optimizataion of both IRC (forward, reverse) trajectory
        WorkFlow.exe(self, 'submitted_04.txt', IRCopt)
