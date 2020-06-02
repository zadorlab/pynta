import os
import sys
import re
from pathlib import Path
from ase.io import read, write


def whichToRestart(PathToSpecies):
    print(PathToSpecies)
    pathlist = Path(PathToSpecies).glob('**/*final.png')
    finishedCalc = []
    for path in pathlist:
        path = os.path.split(str(path))
        ''' prefix for finished calculations for given species'''
        prefix = path[1][:2]
        finishedCalc.append(prefix)
    # print(sorted(finishedCalc, key = str))
    # print(sorted(next(os.walk(PathToSpecies))[1], key = str))
    ''' All conformers - finished calc '''
    ToBeRestarted = set(next(os.walk(PathToSpecies))[1]) - set(finishedCalc)
    if checkAllFinished(PathToSpecies) == True:
        print('All calculations finished succesfully')
    elif not ToBeRestarted:
        print('Calculations in progress or all need to be restarted!')
        pass
        # else:
        #     # sys.exit()
        # pass
    else:
        print('Some of calculations require restart')
        ToBeRestarted = sorted(ToBeRestarted)
        for i_prefix in ToBeRestarted:
            print(i_prefix + ' needs restart')
        print('Restarting...')
        # print(ToBeRestarted[int(i_prefix)])
        # print('{} needs restart'.format(ToBeRestarted[int(i_prefix)]))

    return ToBeRestarted


def createRestartJobs(PathToSpecies):
    restarting_prefix_list = whichToRestart(PathToSpecies)
    for prefix in restarting_prefix_list:
        PathToConformer = os.path.join(PathToSpecies, prefix + '.traj')
        newXYZ = os.path.join(PathToSpecies, prefix + '.xyz')
        # write(newXYZ, read(PathToConformer))
    return restarting_prefix_list
    # print(PathToConformer)


def createRestartJobsTS(PathToSpecies):
    restarting_prefix_list = whichToRestart(PathToSpecies)
    pathlist = Path(PathToSpecies).glob('**/*traj')
    rxn_names = []
    for path in pathlist:
        rxn_name_tmp = os.path.split(str(path))[1]
        rxn_name_tmp = re.search('\_(.+?)\.', rxn_name_tmp).group(1)
        rxn_names.append(rxn_name_tmp)
    rxn = list(set(rxn_names))[0]

    for prefix in restarting_prefix_list:
        PathToConformer = os.path.join(
            PathToSpecies, prefix, prefix + '_' + rxn + '.traj')
        newXYZ = os.path.join(PathToSpecies, prefix,
                              prefix + '_' + rxn + '_ts.xyz')
        write(newXYZ, read(PathToConformer))
    return restarting_prefix_list, rxn


def checkAllFinished(PathToSpecies):
    ''' Get dirs '''
    dirlist = []
    finalpnglist = []
    for dir in os.listdir(PathToSpecies):
        dirpath = os.path.join(PathToSpecies, dir)
        if os.path.isdir(dirpath):
            dirlist.append(dir)
        elif dir.endswith('_final.png'):
            finalpnglist.append(dir)
    nDir = len(dirlist)  # number of dirs
    npng = len(finalpnglist)  # number of *final.png files
    if nDir == npng:
        return True
    else:
        # print('Some of calculations require restart')
        return False


def restart(species, PathToSpecies, facetpath):
    restarting_prefix_list = createRestartJobs(PathToSpecies)
    for prefix in restarting_prefix_list:
        # print(prefix)
        fname = species + '_' + prefix + '_relax.py'
        command = "cd ./{}/minima/; sbatch {} >> ../../submitted_01.txt; cd ../../".format(
            facetpath, fname)
        bashCommand = os.popen(command)
        print(bashCommand.read())


def restartTS(PathToSpecies, facetpath):
    restarting_prefix_list, rxn = createRestartJobsTS(PathToSpecies)
    # rxn = TSxyz.split('_')[1]
    for prefix in restarting_prefix_list:
        fname = prefix + '_' + rxn + '_ts.py'
        command = "cd ./{}/TS_estimate_unique/; sbatch {} >> ../../submitted_03.txt; cd ../../".format(
            facetpath, fname)
        # print(command)
        bashCommand = os.popen(command)
        print(bashCommand.read())
