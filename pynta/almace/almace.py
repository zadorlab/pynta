# Define the base directory for the new location
import os
import sys
import shutil
import json
import yaml

from alCalc import AlMaceCalculator
from mace.calculators.mace import MACECalculator

from ase import units
from ase.io import read, write, Trajectory
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.calculators.calculator import Calculator, all_changes
from importlib import import_module


def name_software(software_name):
    """
    go from software_name to the associated
    ASE calculator constructor
    """
    if software_name == "XTB":
        module = import_module("xtb.ase.calculator")
        return getattr(module, software_name)
    elif software_name == "PWDFT":
        module = import_module("pynta.ase_pwdft.pwdft")
        return getattr(module, software_name)
    else:
        module = import_module("ase.calculators."+software_name.lower())
        return getattr(module, software_name)


def read_yaml(file_path):
    """Reads a YAML file and returns the contents as a dictionary."""
    with open(file_path, 'r') as yaml_file:
        return yaml.safe_load(yaml_file)


def set_ALMACE(file_xyz, work_dir, storage, sub_calc, sub_calc_kwargs_yaml):

    config = read_yaml("/lus/eagle/projects/catalysis_aesp/raymundohe/testPyntaMultiNode/testOpt/config.yaml")
    sub_calc_kwargs = read_yaml(sub_calc_kwargs_yaml)

    if storage == 'None':
        storage = None

    training_set = loadStructures(storage)

    origin_subdir = config['origin']
    local_arg_foundational  =  os.path.join(origin_subdir, config['args_foundational'])

    # 1. Read in configuration to run MD on
    #fname = 'train-it0.xyz'
    #at = read(fname,4) # original
    # ya no lo necesito at = read(file_xyz) # modified by RHE
    args = read_yaml(local_arg_foundational)


    # 2. Define DFT calculator
    # In this test case: take already trained model as DFT calculator

    calc = name_software(sub_calc)(**sub_calc_kwargs) ## Define the sub calculator

    dir = os.path.join(work_dir, 'current')
    src_dir = dir
    #base_new_dir = os.path.join(config['al_directory'],'old-')


    # 3. Change directory to a new directlry
    # Just for testing
    #new_dir = base_new_dir + str(1)
    #while os.path.exists(new_dir):
    #    new_dir = base_new_dir + str(int(new_dir.split('-')[-1]) + 1)

    # Move the entire directory tree to the new location
    #if os.path.exists(src_dir):
    #    shutil.move(src_dir, new_dir)
    #os.makedirs(src_dir, exist_ok=True)
    # move to directory

    # This is the Calculator
    almace = AlMaceCalculator(
        AL_dir=dir,
        dft_calculator=calc,
        mlff_parameters=args,
        force_al_threshold=config['force_threshold'], # eV/A force error
        rel_force_al_threshold=config['rel_force_threshold'], # relative force error
        num_committes=2,
        #initial=config['initial'],
        #initial=None,
        device=config['device'],
        default_dtype='float32',
        #initial_atom=None,
        initial_atom=training_set,
        history = True,
        storage = storage,
    )


    return almace



def run_molecule(calc, name_xyz):

    molecule = read(name_xyz)

    molecule.calc = calc

    energy = molecule.get_potential_energy()

    forces = molecule.get_forces()

    name_json = name_xyz.replace(".xyz", ".json")

    data = {'energy': energy,
            'forces': forces.tolist(),
            'symbols': molecule.get_chemical_symbols(),
            'positions': molecule.get_positions().tolist()}

    if hasattr(molecule, 'cell'):
        data['cell'] = molecule.cell.tolist()
    if hasattr(molecule, 'pbc'):
        data['pbc'] = molecule.pbc.tolist()


    with open(name_json, "w") as file:
        json.dump(data, file)

def loadStructures(subdirectory):
    structures = []

    if not os.path.isdir(subdirectory):
        print("The subdirectory does not exist")
        return structures

## Load the xyz files
    iter=0
    for fname in os.listdir(subdirectory):
        if fname.endswith(".xyz"):
            filepath = os.path.join(subdirectory, fname)
            mol = read(filepath)
            structures.append(mol)
            iter += 1

## For MDMin the workflow creates a traj file with pure QE Calculations
    for fname in os.listdir(subdirectory):
        if fname.endswith(".traj"):
            filepath = os.path.join(subdirectory, fname)
            traj = Trajectory(filepath)
            subiter=0
            for mol in traj:
                structures.append(mol)
                subiter += 1
            iter += 1

    return structures


if  __name__ == "__main__":

    print(sys.argv)
    if len(sys.argv) == 6:
        #work_dir = sys.argv[2]
        work_dir = sys.argv[3]
        storage = sys.argv[3]
        sub_calc = sys.argv[4]
        sub_calc_kwargs_yaml = sys.argv[5]

        print(sub_calc_kwargs_yaml)

        almace = set_ALMACE(sys.argv[1], work_dir, storage, sub_calc, sub_calc_kwargs_yaml)
        run_molecule(almace, sys.argv[1])
    else:
        print(" ERROR to invoque ALMACE")

