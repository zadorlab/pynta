#!/usr/bin/env python3
#SBATCH -J main_job_Cu_111
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p week-long-cpu
#SBATCH -t 7-00:00:00
#SBATCH -e %x.err
#SBATCH -o %x.out
'''
####################################################
###################### Set up ######################
####################################################
'''
import os
import sys
submitDir = os.environ['SLURM_SUBMIT_DIR']
os.chdir(submitDir)
sys.path.append(os.getcwd())
import time
import inputR2S
from rmgcat_to_sella.main import run
from rmgcat_to_sella.main import exe
from rmgcat_to_sella.main import genJobFiles
from rmgcat_to_sella.main import check_if_minima_already_calculated
'''
####################################################
#################### Initialize ####################
####################################################
'''
genJobFiles()
SurfaceOpt       = inputR2S.SurfaceScript
SurfaceAdsorbate = inputR2S.SurfaceAdsorbateScript
TSxtb            = inputR2S.TSxtbScript
TS               = inputR2S.TSScript
IRC              = inputR2S.IRCScript
IRCopt           = inputR2S.IRCoptScript
##
currentDir       = os.path.dirname(os.getcwd())
sp1              = inputR2S.sp1
sp2              = inputR2S.sp2
facetpath        = inputR2S.facetpath
'''
####################################################
##################### Workflow #####################
####################################################
'''
'''
Run initial optimization of reactantts and products.
Among others, this will generate the average bond distances.
Minimize penalty function + xTB total energy
'''

checksp1 = check_if_minima_already_calculated(currentDir, sp1, facetpath)
checksp2 = check_if_minima_already_calculated(currentDir, sp2, facetpath)
if checksp1 is False and checksp2 is False):
    runSurfAds = os.popen(os.path.join(

    'sbatch ' + SurfaceAdsorbate))
    print(runSurfAds.read())
    while not os.path.exists('submitted_01.txt'):
        time.sleep(1)
    run('submitted_01.txt', TSxtb)
else:
    # If all minimas were calculated some time age for other reaction,
    # rmgcat_to_sella will use that calculations.
    run_TSxtb = os.popen(os.path.join(
    'sbatch ' + TSxtb))
    print(run_TSxtb.read())

'''
Check for symmetry distinct sites,
Use them to set up saddle point optimization in Sella,
Run calculations.
'''
exe('submitted_02.txt', TS)
'''
Run IRC calculations
'''
exe('submitted_03.txt', IRC)
'''
Minimize forward and reverse IRC trajectories
'''
exe('submitted_04.txt', IRCopt)
####################################################
####################################################
####################################################
