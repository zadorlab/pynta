import os

import shutil

from ase.io import read, write

from rmgcat_to_sella.adjacency_to_3d import rmgcat_to_gratoms

def checkEner(path, treshHold):
    allSpecies = next(os.walk(path))[1]

    dict_1 = dict()
    dict_2 = dict()
    dict_3 = dict()

    for file in sorted(os.listdir(path), key=str):
        if file.endswith('.out'):
            if file.startswith(str(allSpecies[0]) + '_'):
                tmpSplitter = file.split('_')
                key = str(tmpSplitter[0] + '_' + tmpSplitter[1])
                with open(os.path.join(path, file)) as f:
                    data = f.readlines()
                    enerLine = data[-1]
                    enerVal = enerLine.split()
                    dict_1[key] = float(enerVal[3])                
            elif file.startswith(str(allSpecies[1]) + '_'):
                tmpSplitter = file.split('_')
                key = str(tmpSplitter[0] + '_' + tmpSplitter[1])
                with open(os.path.join(path, file)) as f:
                    data = f.readlines()
                    enerLine = data[-1]
                    enerVal = enerLine.split()
                    dict_2[key] = float(enerVal[3])  
            elif file.startswith(str(allSpecies[2]) + '_'):
                tmpSplitter = file.split('_')
                key = str(tmpSplitter[0] + '_' + tmpSplitter[1])
                with open(os.path.join(path, file)) as f:
                    data = f.readlines()
                    enerLine = data[-1]
                    enerVal = enerLine.split()
                    dict_3[key] = float(enerVal[3])

    for item in dict_1.items():
        comparator = item[1] - min(dict_1.values()) # calculate energy difference
        print(item[0], '     ', "{0:.2e}".format(comparator), 'eV     ', "{:8.3f}".format(comparator * 23.06035 * 4.184) + ' kJ/mol')
        if comparator > treshHold or comparator == 0:
            species, prefix = item[0].split('_')
            copyDestDir = os.path.join(path + '_unique_ener_based', species)
            os.makedirs(copyDestDir, exist_ok=True)
            copySrcFile = os.path.join(path, species, prefix + '.xyz')
            shutil.copy2(copySrcFile, copyDestDir)
            copySrcFile_traj = os.path.join(path, species, prefix + '.traj')
            write(os.path.join(copyDestDir, prefix + '.png'), read(copySrcFile_traj))
    print()

    for item in dict_2.items():
        comparator = item[1] - min(dict_2.values()) # calculate energy difference
        print(item[0], '     ', "{0:.2e}".format(comparator), 'eV     ', "{:8.3f}".format(comparator * 23.06035 * 4.184) + ' kJ/mol')
        if comparator > treshHold or comparator == 0:
            species, prefix = item[0].split('_')
            copyDestDir = os.path.join(path + '_unique_ener_based', species)
            os.makedirs(copyDestDir, exist_ok=True)
            copySrcFile = os.path.join(path, species, prefix + '.xyz')
            shutil.copy2(copySrcFile, copyDestDir)
            copySrcFile_traj = os.path.join(path, species, prefix + '.traj')
            write(os.path.join(copyDestDir, prefix + '.png'), read(copySrcFile_traj))
    print()

    for item in dict_3.items():
        comparator = item[1] - min(dict_3.values()) # calculate energy difference
        print(item[0], '     ', "{0:.2e}".format(comparator), 'eV     ', "{:8.3f}".format(comparator * 23.06035 * 4.184) + ' kJ/mol')
        if comparator > treshHold or comparator == 0:
            species, prefix = item[0].split('_')
            copyDestDir = os.path.join(path + '_unique_ener_based', species)
            os.makedirs(copyDestDir, exist_ok=True)
            copySrcFile = os.path.join(path, species, prefix + '.xyz')
            shutil.copy2(copySrcFile, copyDestDir)       
            copySrcFile_traj = os.path.join(path, species, prefix + '.traj')
            write(os.path.join(copyDestDir, prefix + '.png'), read(copySrcFile_traj))

    # def findAllNebsEner():

    checkEner(path, treshHold)

