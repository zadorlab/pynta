from ase.calculators.calculator import Calculator
from ase.io import write, read
import subprocess

class WrapperALMACE(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

        if 'python_mace' in kwargs:
            self.python_mace = kwargs['python_mace']
        else:
            raise ValueError 
        
        self.kwargs = kwargs 
        
    def calculate(self, atoms, properties, system_changes):
        super().calculate(self, atoms, properties, system_changes)
        
        write(atoms,"init_temp.xyz")
        
        python_script = f"""from alCalc import AlMaceCalculator
        from ase.io import read, write
        atoms = read("init_temp.xyz")
        calc = AlMaceCalculator(**{self.kwargs})
        atoms.calc = calc
        atoms.get_potential_energy()
        atoms.get_forces()
        write(atoms,"out.xyz")
        """
        
        with open("temp_script.py",'w') as f:
            f.write(python_script)
        
        command = f'{self.python_mace} temp_script.py'
        subprocess.run(command, shell=True)

        out_atoms = read("out.xyz")
        self.results['energy'] = out_atoms.info['energy']
        self.results['forces'] = out_atoms.arrays['forces']
        
        return 
    
        