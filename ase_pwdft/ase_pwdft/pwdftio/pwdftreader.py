import re
import json
import numpy as np

from ase import Atoms
from ase.units import Hartree, Bohr
from ase.calculators.singlepoint import SinglePointDFTCalculator


def read_pwdft_out(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    matches = []
    for idx, line in enumerate(lines):
        if re.findall(r'>>> job completed', line):
            matches.append(idx)

    last_match = matches[-1]
    line = lines[last_match + 2]

    if 'Next rtdbstr' in line:
        json_str = line.split('Next rtdbstr=')[-1]

        data = json.loads(json_str)

        pspw = data['pspw']
        energy = pspw['energy'] * Hartree

        geo = data['geometries']
        geo1 = geo['geometry']
        nion = geo1['nion']

        symbols = geo1['symbols']
        cell = np.array(geo1['unita']).reshape(3, 3)
        nwpw = data['nwpw']
        dipole = np.array(nwpw['dipole'])
        coors = np.array(geo1['coords']).reshape(nion, 3)

        atoms = Atoms(symbols, positions=coors, cell=cell)

        forces = np.zeros((nion, 3))

        if 'fion' in pspw:
            forces = np.array(pspw['fion']).reshape(nion, 3)
        forces *= Hartree / Bohr

        calc = SinglePointDFTCalculator(atoms=atoms,
                                        energy=energy,
                                        free_energy=energy,  # XXX Is this right?
                                        forces=forces,
                                        dipole=dipole,
                                        # quadrupole=quadrupole,
                                        )
        calc.kpts = 1
        atoms.calc = calc
        return atoms
    else:
        pass
