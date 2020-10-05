import os
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
    # Getting unfinished calculations, i.e  all conformers - finished calc
    print(finishedCalc)
    ToBeRestarted = set(next(os.walk(PathToSpecies))[1]) - set(finishedCalc)
    if checkAllFinished(PathToSpecies):
        print('All calculations finished succesfully')
    elif not ToBeRestarted:
        print('Calculations in progress or all need to be restarted!')
        pass
    else:
        print('Some of calculations require restart')
        ToBeRestarted = sorted(ToBeRestarted)
        for i_prefix in ToBeRestarted:
            print(i_prefix + ' needs restart')

    return ToBeRestarted


def createRestartJobs(PathToSpecies):
    restarting_prefix_list = whichToRestart(PathToSpecies)
    for prefix in restarting_prefix_list:
        PathToConformer = os.path.join(PathToSpecies, prefix + '.traj')
        newXYZ = os.path.join(PathToSpecies, prefix + '.xyz')
        print(PathToConformer)
        print(newXYZ)
        # write(newXYZ, read(PathToConformer))
    print('Restarting...')
    return restarting_prefix_list


def createRestartJobsTS(PathToSpecies):
    restarting_prefix_list = whichToRestart(PathToSpecies)
    pathlist = Path(PathToSpecies).glob('**/*traj')
    rxn_names = []
    for path in pathlist:
        rxn_name_tmp = os.path.split(str(path))[1]
        rxn_name_tmp = re.search(r'\_(.+?)\.', rxn_name_tmp).group(1)
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
        return False


def restart(species, PathToSpecies, facetpath):
    restarting_prefix_list = createRestartJobs(PathToSpecies)
    from balsam.launcher.dag import BalsamJob
    yamlfile = '{yamlfile}'
    workflow_name = yamlfile + facetpath + '01'
    for prefix in restarting_prefix_list:
        # print(prefix)
        py_script = os.path.join(
            facetpath, 'minima', '_'.join([species, prefix, 'relax.py'])
        )
        job_to_add = BalsamJob(
            name=py_script,
            workflow=workflow_name,
            application='python',
            args=py_script,
            ranks_per_node=1,
        )
        job_to_add.save()


def restartTS(PathToSpecies, facetpath):
    restarting_prefix_list, rxn = createRestartJobsTS(PathToSpecies)
    # rxn = TSxyz.split('_')[1]
    from balsam.launcher.dag import BalsamJob
    yamlfile = '{yamlfile}'
    workflow_name = yamlfile + facetpath + '03'
    for prefix in restarting_prefix_list:
        # print(prefix)
        py_script = os.path.join(
            facetpath, 'TS_estimate_unique', '_'.join([prefix, rxn, 'ts.py'])
        )
        job_to_add = BalsamJob(
            name=py_script,
            workflow=workflow_name,
            application='python',
            args=py_script,
            ranks_per_node=1,
        )
        job_to_add.save()

# TODO: Do I really need that?


def whichToRestart_optIRC(PathToSpecies):
    print(PathToSpecies)
    pathlist_finished = Path(PathToSpecies).glob('**/*final.png')
    finishedCalc = []
    for path in pathlist_finished:
        path = os.path.split(str(path))
        finishedCalc.append(path[0])
    # These are all calculation that were started but not finished
    allCalc = []
    pathlist_all = Path(PathToSpecies).glob('**/*relax*')
    for path in pathlist_all:
        path = os.path.split(str(path))
        allCalc.append(path[0])

    ToBeRestarted = (sorted(set(allCalc) - set(finishedCalc)))
    return ToBeRestarted


def createRestartJobs_optIRC(PathToSpecies):
    restarting_list = whichToRestart_optIRC(PathToSpecies)
    pathlist = Path(PathToSpecies).glob('**/irc*/*traj')
    rxn_names = []
    for path in pathlist:
        rxn_name_tmp = os.path.split(str(path))[1]
        rxn_name_tmp = re.search(r'\_(.+?)\.', rxn_name_tmp).group(1)
        rxn_names.append(rxn_name_tmp)
    rxn = list(set(rxn_names))[0]
    # print(rxn_names)

    for ircPath in restarting_list:
        _, _, prefix, irc = ircPath.split('/')
        PathToConformer = os.path.join(
            ircPath, prefix + '_' + rxn + '.traj')
        newXYZ = os.path.join(ircPath, prefix + '_' +
                              rxn + '_' + irc[:-4] + '.xyz')
        write(newXYZ, read(PathToConformer))
    print('Restarting...')
    return restarting_list, rxn


def restart_optIRC(PathToSpecies):
    restarting_list, rxn = createRestartJobs_optIRC(PathToSpecies)
    from balsam.launcher.dag import BalsamJob
    yamlfile = '{yamlfile}'
    workflow_name = yamlfile + facetpath + '03'
    for ircPath in restarting_list:
        _, _, prefix, irc = ircPath.split('/')
        py_script = os.path.join(ircPath, '_'.join([prefix, rxn, irc + '.py']))
        """ This may be wrong, probably needs debugging"""
        job_to_add = BalsamJob(
            name=py_script,
            workflow=workflow_name,
            application='python',
            args=py_script,
            ranks_per_node=1,
        )
        job_to_add.save()
