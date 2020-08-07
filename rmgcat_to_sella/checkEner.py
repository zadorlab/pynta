import os

import shutil

from ase.io import read, write
from ase.units import kJ, mol


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

    for collection in (dict_1, dict_2, dict_3):
        e_min = min(collection.values())
        for name, energy in collection.items():
            delta_e = energy - e_min
            print('{}     {:.2e}eV     {:8.3f} kJ/mol'.format(
                energy,
                delta_e,
                delta_e * mol / kJ
            ))
            if delta_e > treshHold or delta_e == 0:
                species, prefix = name.split('_')
                copyDestDir = os.path.join(
                    path + '_unique_ener_based', species
                )
                os.makedirs(copyDestDir, exist_ok=True)
                copySrcFile = os.path.join(path, species, prefix + '.xyz')
                shutil.copy2(copySrcFile, copyDestDir)
                copySrcFile_traj = os.path.join(
                    path, species, prefix + '.traj'
                )
                write(
                    os.path.join(copyDestDir, prefix + '.png'),
                    read(copySrcFile_traj)
                )
        print()

    # def findAllNebsEner():

    checkEner(path, treshHold)
