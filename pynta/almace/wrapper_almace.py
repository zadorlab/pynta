from ase.calculators.calculator import Calculator
from ase.io import write, read
from ase.units import Bohr
from pynta.utils import name_to_ase_software
import os
import json
import subprocess
import shutil
import numpy as np
import yaml

class WrapperALMACE(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)
        # self.python_bin = '/lus/eagle/projects/catalysis_aesp/raymundohe/maceFlow/mace_env311/bin/python'
        # self.almace = '/lus/eagle/projects/catalysis_aesp/raymundohe/testPyntaMultiNode/testOpt/almace.py'

        if 'host' in kwargs:
            self.host = kwargs['host']
            self.connect = True
        else:
            self.host = 'localhost'
            self.connect = False

        if 'opt_method' in kwargs:
            self.opt_method = kwargs['opt_method']
        else:
            self.opt_method = None

        if 'sub_software' in kwargs:
            self.sub_software = kwargs['sub_software']
        else:
            self.sub_software = 'Espresso'

        if 'sub_software_kwargs' in kwargs:
            self.sub_software_kwargs = kwargs['sub_software_kwargs']
        else:
            print(" You need a DFT software for ALMACE")

        if 'storage' in kwargs:
            self.storage = kwargs['storage']
        else:
            self.storage = None

        if 'debug' in kwargs:
            self.debug = kwargs['debug']
        else:
            self.debug = False

        if 'python_mace' in kwargs:
            self.python_bin = kwargs['python_mace']
        else:
            raise ValueError 
        
        if 'almace_path' in kwargs:
            self.almace = kwargs['almace_path']
        else:
            raise ValueError 
        
        if 'min_data_to_train' in kwargs:
            self.min_data_to_train = kwargs["min_data_to_train"]
        else:
            self.min_data_to_train = 0
            
    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        # Create the input
        atoms.write('input.xyz')

        cwd_path = os.getcwd()

        if self.sub_software_kwargs is not None:
            file_yaml = os.path.join(cwd_path,'sub_software.yaml')

            with open(file_yaml, 'w') as file:
                yaml.dump(self.sub_software_kwargs, file)

        Ndata = len(os.listdir(self.storage))
        if Ndata >= self.min_data_to_train:
            self.training = False
        else:
            self.training = True

        if self.training == False:
            import subprocess
            command = f'{self.python_bin} {self.almace} input.xyz {cwd_path} {self.storage} {self.sub_software} sub_software.yaml'
            subprocess.run(command, shell=True)

            with open("input.json", "r") as file:
                data_json = json.load(file)
            self.results['energy'] = data_json['energy']
            self.results['forces'] = np.array(data_json['forces'])
            atoms.info['energy'] = self.results['energy']
            atoms.arrays['forces'] = self.results['forces']
        else:
            atoms.calc = name_to_ase_software(self.sub_software)(**self.sub_software_kwargs)
            atoms.info['energy'] =  atoms.get_potential_energy()
            atoms.arrays['forces'] = atoms.get_forces()
            self.results['energy'] = atoms.info['energy']
            self.results['forces'] = atoms.arrays['forces']

            write('almaceDFT000.xyz', atoms)

            if self.storage is not None:
                if not os.path.exists(self.storage):
                    os.makedirs(self.storage)

                src = os.path.join(os.getcwd(), "almaceDFT000.xyz");
                dst = os.path.join(self.storage, "almaceDFT000.xyz");

                if os.path.exists(dst):
                    count = 1
                    while True:
                        new_name = f"almaceDFT{count:03d}.xyz"
                        new_path = os.path.join(self.storage, new_name)
                        if not os.path.exists(new_path):
                            dst = new_path
                            break
                        count +=1
                shutil.copy(src, dst)

        return atoms
