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

####################################################
#################### Initialize ####################
####################################################

# TODO: It would be great to hava a class here


def genJobFiles():
    ''' Generate submt scripts for 5 stages of the workflow '''
    set_up_slab(template_slab_opt, surface_type, symbol, a, repeats_surface,
                vacuum, slab_name, pseudopotentials)
    set_up_ads(template_ads, facetpath, slabopt,
               yamlfile, repeats, pytemplate_relax_ads)
    set_up_TS_with_xtb(template_set_up_ts_with_xtb,
                       slabopt, repeats, yamlfile, facetpath,
                       rotAngle, scfactor, scfactor_surface, pytemplate_xtb,
                       sp1, sp2)
    set_up_run_TS(template_set_up_ts, facetpath, pytemplate_set_up_ts)
    set_up_run_IRC(template_set_up_IRC, facetpath,
                   pytemplate_f, pytemplate_r, yamlfile)
    set_up_opt_IRC(template_set_up_optIRC, facetpath, pytemplate_optIRC)


###########################
#   Create submit files   #
###########################

def set_up_slab(template, surface_type, symbol, a, repeats_surface, vacuum,
                slab_name, pseudopotentials):
    ''' Create 00_set_up_slab_opt.py file '''
    with open(template, 'r') as r:
        template = r.read()
        with open('00_set_up_slab_opt.py', 'w') as c:
            c.write(template.format(surface_type=surface_type, symbol=symbol,
                                    a=a, repeats_surface=repeats_surface,
                                    vacuum=vacuum, slab_name=slab_name,
                                    pseudopotentials=pseudopotentials))
        c.close()
    r.close()


def set_up_ads(template, facetpath, slabopt, yamlfile, repeats, pytemplate):
    ''' Create 01_set_up_ads.py file '''
    with open(template, 'r') as r:
        template = r.read()
        with open('01_set_up_ads.py', 'w') as c:
            c.write(template.format(facetpath=facetpath, slabopt=slabopt,
                                    yamlfile=yamlfile, repeats=repeats,
                                    pytemplate=pytemplate))
        c.close()
    r.close()


def set_up_TS_with_xtb(template, slab,
                       repeats, yamlfile, facetpath, rotAngle,
                       scfactor, scfactor_surface,
                       pytemplate_xtb, sp1, sp2):
    ''' Create 02_set_up_TS_with_xtb.py file'''
    with open(template, 'r') as r:
        template = r.read()
        with open('02_set_up_TS_with_xtb.py', 'w') as c:
            c.write(template.format(facetpath=facetpath, slab=slab,
                                    repeats=repeats, yamlfile=yamlfile,
                                    rotAngle=rotAngle, scfactor=scfactor, scfactor_surface=scfactor_surface,
                                    pytemplate_xtb=pytemplate_xtb, sp1=sp1,
                                    sp2=sp2, scaled1=scaled1, scaled2=scaled2))
        c.close()
    r.close()


def set_up_run_TS(template, facetpath, pytemplate):
    ''' Create 03_checksym_xtb_runTS.py file '''
    with open(template, 'r') as r:
        template = r.read()
        with open('03_checksym_xtb_runTS.py', 'w') as c:
            c.write(template.format(facetpath=facetpath,
                                    pytemplate=pytemplate))
        c.close()
    r.close()


def set_up_run_IRC(template, facetpath, pytemplate_f, pytemplate_r, yamlfile):
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


def set_up_opt_IRC(template, facetpath, pytemplate):
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


def get_slurm_jobs_id(slurm_id_subm):
    ''' Get slurm IDs of just submitted jobs '''
    slurm_jobs_id = []
    with open(slurm_id_subm, 'r') as f:
        for line in f.readlines():
            line = line.split()[3]
            slurm_jobs_id.append(line)
    f.close()
    print(slurm_jobs_id)
    return slurm_jobs_id


def genSbatchCommand(slurm_id_subm):
    ''' Prepare a bash command to submit jobs '''
    slurmID = get_slurm_jobs_id(slurm_id_subm)
    slurmID = ",".join(["{}"] * len(slurmID)).format(*slurmID)
    # if not slurmID:
    #     print('Error')
    #     sys.exit('No submitted jobs, probably all files have been already generated')
    #     return command = False
    # else:
    command = os.path.join('sbatch --dependency=afterany:' + str(slurmID))
    print(command)
    return command


def run(slurm_id_subm, job_script):
    ''' Submit slurm jobs '''
    command = genSbatchCommand(slurm_id_subm)
    os.popen(str(os.path.join(command + ' ' + job_script)))


def exe(prevSlurmID, job_script):
    ''' Check if the previous step of calculations terminated. If so, run the next step'''
    while not os.path.exists(prevSlurmID):
        time.sleep(60)
    run(prevSlurmID, job_script)


def check_if_path_to_mimina_exists(mainDir, species):
    ''' Check for the paths to previously calculated minima and return 
        a list with all valid paths '''

    pathlist = Path(mainDir).glob('**/minima/' + species)
    paths = []
    for path in pathlist:
        # path = str(path)
        paths.append(str(path))
        return paths[0]
    if IndexError:
        return None


def check_if_minima_already_calculated(currentDir, species, facetpath):
    ''' Check for previously calculated minima '''
    mainDirs = []
    uniqueMinimaDirs = []
    if facetpath == 'Cu_211':
        mainDirsList = Path(str(currentDir)).glob('*_Cu_211_methanol*')
    elif facetpath == 'Cu_100':
        mainDirsList = Path(str(currentDir)).glob('*_Cu_100_methanol*')
    else:
        mainDirsList = Path(str(currentDir)).glob('*_Cu_methanol*')
    # transforming posix path to regular string
    for mainDir in mainDirsList:
        mainDirs.append(mainDir)
    # expected e.g. -> mainDirs = ['00_Cu_methanol_CO+O_CO2',
    # '01_Cu_methanol_OH_O+H', '02_Cu_methanol_CO+H_HCO']
    for mainDir in mainDirs:
        minimaDir = check_if_path_to_mimina_exists(mainDir, species)
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
