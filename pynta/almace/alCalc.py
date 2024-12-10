from ase.calculators.calculator import Calculator, all_changes
import datetime
from ase.io import read, write
import os
from mace.calculators.mace import MACECalculator
from pynta.utils import import_module

import argparse
import sys
import shutil

import numpy as np
import logging
from glob import glob
from uuid import uuid4

import mace.cli.run_train

def name_to_ase_software(software_name):
    """
    go from software_name to the associated
    ASE calculator constructor
    """
    if software_name == "XTB":
        module = import_module("xtb.ase.calculator")
        return getattr(module, software_name)
    else:
        module = import_module("ase.calculators."+software_name.lower())
        return getattr(module, software_name)

class AlMaceCalculator(MACECalculator):
    def __init__(
        self,
        AL_dir,
        dft_calculator_name,
        dft_calculator_kwargs,
        mlff_parameters,
        force_al_threshold,
        rel_force_al_threshold,
        num_committes=4,
        mlff_train_cmd=None,
        logger=None,
        initial=None, # contains a already trained MLFF with correct number of neighbours
        initial_atom=None, # If no initial MLFF is given, use this atom to start training -> evaluate DFT + retrain
        history=False,
        reuse_checkpoints=True,
        storage = None,
        min_data_to_train = 0, 
        **kwargs,
    ):
        """Initialize the AL-MACE calculator with given parameters."""
        self.kwargs = kwargs

        if initial is not None and not os.path.exists(AL_dir):
            print(initial)
            print(AL_dir)
            #self.logger.info('Copying files from initial to current directory')
            #self.logger.info(f'Initial: {initial}, AL_dir: {AL_dir}')
            shutil.copytree(initial, AL_dir)
            print(" Are we copying??")

        if not os.path.exists(AL_dir):
            os.mkdir(AL_dir)

        if logger is None:
            self.logger = logging.getLogger(__name__)
            # remove all handlers
            for handler in self.logger.handlers[:]:  # Get a copy of the list to avoid modification issues
                self.logger.removeHandler(handler)
            if True:#not self.logger.hasHandlers():
                print('Adding handlers')
                handler = logging.StreamHandler()
                # self.logger.addHandler(handler)
                # Create a FileHandler
                log_fname = os.path.join(AL_dir, f'logfile-{uuid4().hex}.log')
                handler = logging.FileHandler(log_fname)
                formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
                handler.setFormatter(formatter)
                self.logger.addHandler(handler)
                self.logger.info(f'Logging to {log_fname}')
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger = logger

        self.logger.info("Initializing AL-MACE calculator: {}".format(AL_dir))



        # Name directories
        self.dir_AL = AL_dir
        self.dir_train = os.path.join(self.dir_AL, "TRAIN")
        self.dir_mlff = os.path.join(self.dir_AL, "MLFF")
        self.storage = storage
        # os.chdir(self.dir_AL)

        self.history = history
        self.at_history = []
        self.reuse_checkpoint = reuse_checkpoints

        # make all folders
        os.makedirs(self.dir_AL, exist_ok=True)
        os.makedirs(self.dir_train, exist_ok=True)
        [os.makedirs(os.path.join(self.dir_train,i), exist_ok=True) for i in ['new', 'train']]
        os.makedirs(self.dir_mlff, exist_ok=True)
        self.dir_initial_working = os.getcwd()

        # other inputs
        self.dft_calculator = name_to_ase_software(dft_calculator_name)(**dft_calculator_kwargs)
        self.mlff_parameters = mlff_parameters
        self.mlff_train_cmd = mlff_train_cmd

        # AL parameters
        self.force_al_threshold = force_al_threshold
        self.rel_froce_al_threshold = rel_force_al_threshold
        self.num_committes = num_committes
        self.min_data_to_train = min_data_to_train
        
        self.time_format = "%Y%m%d%H%M%S"
        # set to early time
        self.timestamp_fail = datetime.datetime.min.strftime(self.time_format)
        self.timestamp_train = self.timestamp_fail
        
        print(self.timestamp_fail)

        if len(self.get_fname_mlffs()) < self.num_committes and initial_atom is None:
            self.logger.info(
                f'Not enough MLFFs to run committee: {len(self.get_fname_mlffs())} < {self.num_committes}'
            )
            self.logger.info(
                f'Make sure the initial model directory ({initial}) has sufficient amount of models.' + 'Or make sure you supply an initial traning set to initial_atom',
                )

        mace_fnames = self.update_mlffs()
        
        super().__init__(mace_fnames, **kwargs)

    def get_fname_mlffs(self, current=True):
        """Get filenames of MLFF models."""
        mlff_fname_pat = os.path.join(self.dir_mlff, "{0}_{1}/*_swa.model")
        mlff_fnames = np.array(sorted(glob(mlff_fname_pat.format("*", "*"))))
        if len(mlff_fnames) == 0:
            return []

        if current:
            timestamps = np.array([int(fname.split('/')[-2].split('_')[1]) for fname in mlff_fnames])

            mlff_fnames = mlff_fnames[timestamps >= int(self.timestamp_train)]

        return mlff_fnames

    def calculate(self, atoms=None, properties=["energy", "forces"], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        if self.history:
            self.at_history.append(atoms.copy())
            self.at_history = self.at_history[-10:] # keep memory low
        
        Ndata = len(os.listdir(self.storage))
        
        if Ndata < self.min_data_to_train:
            new_at = self.calculate_dft(atoms)
            self.logger.info(f'Running DFT')
            self.create_new_training(new_at)
            self.results["forces"] = new_at.arrays["forces"]
            self.results["energy"] = new_at.info["energy"]
            return 
        
        self.mace.calculate(atoms)
        self.results = self.mace.results

        print(f'      MACE energy {self.results["energy"]}')

        if 'forces_comm' in self.results.keys():
            force_var = np.var(self.results['forces_comm'], axis=0)
            print(" Forces_comm in self.results.keys()")
        else:
            force_var = np.zeros_like(self.results['forces'])
            print(" No forces_comm in self.results.keys()")
        force_std = np.sqrt(np.sum(force_var, axis=1))
        # relative force error
        reg = 0.2
        fnorm = np.linalg.norm(self.results['forces'], axis=1)
        fnorm[fnorm<reg] = reg

        self.logger.debug(f'Max force error: {np.round(np.max(force_std),4)} Rel_force: {np.round(np.max(force_std / fnorm), 4)}')

        self.logger.debug(f'Force std: {np.max(force_std)}')
        

        if np.max(force_std) > self.force_al_threshold or np.max(force_std / fnorm) > self.rel_froce_al_threshold:
            self.logger.info(f'AL threshold exceeded: {np.max(force_std)} > {self.force_al_threshold} ')
            self.logger.info(f'OR Relative threshold {np.max(force_std / fnorm)} > {self.rel_froce_al_threshold}')
            self.timestamp_fail = datetime.datetime.now().strftime(self.time_format)
            new_at = self.calculate_dft(atoms)
            self.logger.info(f'Running DFT')
            self.create_new_training(new_at)
            self.timestamp_train = datetime.datetime.now().strftime(self.time_format)
            self.logger.info(f'Training Force Field: {self.timestamp_train}')
            self.update_mlffs()
            if self.history:
                self.logger.info('Skipping-back 10 steps')
                atoms = self.at_history[0]
                self.atoms.positions = atoms.positions
            self.mace.calculate(atoms)


    def calculate_dft(self, atoms):
        """Calculate properties using DFT."""
        self.dft_calculator.calculate(atoms)
        self.results = self.dft_calculator.results
        atoms.info['dft_energy'] = self.results['energy']
        atoms.arrays['dft_forces'] = self.results['forces']
        atoms.arrays['forces'] = self.results['forces']
        atoms.info['energy'] = self.results['energy']

        write("almaceDFT000.xyz", atoms);

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

    def create_new_training(self, atoms):
        """Save atomic configurations for further training."""
        filename = os.path.join(
            self.dir_train,f"new/new-{self.timestamp_fail}.xyz"
        )
        write(filename, atoms, append=True)
        all_new = glob(os.path.join(self.dir_train,"new/new-*.xyz"))
        traj_train = []
        [traj_train.extend(read(fname, index=':')) for fname in all_new]
        self.logger.info(f'Number of new training configs: {len(traj_train)}')
        filename = os.path.join(self.dir_train, f"train/train-{self.timestamp_fail}.xyz")
        write(filename, traj_train)



        # for each committe create sepeerate train and valid split
        # TODO: cumulative new - for all new configurations created since failure
        # or since last training configuration that we start from.

        indices = np.arange(len(traj_train))
        # for i in number of committees
        for it in range(self.num_committes):
            np.random.seed(it)
            np.random.shuffle(indices)

            # Assuming a 90-10 split. Adjust if needed.
            split_at = int(0.9 * len(indices))
            train_idx, validation_idx = indices[:split_at], indices[split_at:]
            train = [traj_train[i] for i in train_idx]
            validation = [traj_train[i] for i in validation_idx]
            logging.info(f"Split {it+1} | Train indexes: {list(map(traj_train.index, train))} | "
                        f"Validation indexes: {list(map(traj_train.index, validation))}")
            if type(atoms) is list:
                write(filename.replace('.xyz', f'_{it}.xyz'), [*atoms, *traj_train])
                write(filename.replace('.xyz', f'_{it}_val.xyz'), [*atoms, *validation])
            else:
                write(filename.replace('.xyz', f'_{it}.xyz'), [atoms, *traj_train])
                write(filename.replace('.xyz', f'_{it}_val.xyz'), [atoms, *validation])
            self.current_train_fname = filename


    def save_train(self, atoms):
        """Save atomic configurations for further training."""
        filename = os.path.join(
            self.train_dir, "data", f"train/train-{self.timestamp_fail}.xyz"
        )
        write(filename, atoms, append=True)

    def update_mlffs(self):
        """Update MLFF models based on new training data."""
        for seed in range(self.num_committes):
            if len(self.get_fname_mlffs()) < self.num_committes:
                self.train_mace(
                    self.current_train_fname.replace('.xyz', f'_{seed}.xyz'),
                    self.current_train_fname.replace('.xyz', f'_{seed}_val.xyz'),
                    name=f'{self.timestamp_fail}_{self.timestamp_train}_{seed}', seed=seed)
            else:
                break
        assert (
            len(self.get_fname_mlffs()) >= self.num_committes
        ), "Not enough MLFFs to run committee"
        mace_fnames = sorted(self.get_fname_mlffs())[-self.num_committes:]
        self.logger.info("Using MLFFs: {}".format(mace_fnames))
        self.mace = MACECalculator(mace_fnames, **self.mlff_parameters)
        return mace_fnames

    def train_mace(self,train_fname, valid_fname, name, seed):
        # copy latest checkpoint
        if self.reuse_checkpoint:
            fpat_checkpoint = os.path.join(self.dir_mlff, f"*_*/checkpoints/*{seed}_epoch-*_swa.pt")
            possible_checkpoints = glob(fpat_checkpoint)
            print(len(possible_checkpoints))
            print(fpat_checkpoint)
            if len(possible_checkpoints) > 0:
            # /home/lls34/rds/hpc-work/Data/Pynta/AL-calc/Dev/old-9-almost/MLFF/20231109123707_20231109123707/checkpoints/20231109123707_20231109123707_1_run-1_epoch-12_swa.pt

                sel_check = sorted(possible_checkpoints)[-1]
                os.makedirs(os.path.join(self.dir_mlff, f'{self.timestamp_fail}_{self.timestamp_train}/checkpoints'), exist_ok=True)
                to_fname = os.path.join(self.dir_mlff, f'{self.timestamp_fail}_{self.timestamp_train}/checkpoints/{self.timestamp_fail}_{self.timestamp_train}_{seed}_run-{seed}_epoch-0.pt')
                self.logger.info(f'Copying checkpoint: {sel_check} to {to_fname}')

                shutil.copy(sel_check, to_fname)
        # Create an argparse.Namespace object with the arguments
        args = argparse.Namespace(train_file=train_fname,
                                  valid_file=valid_fname,
                                  name=name,
                                  seed=seed,
                                  wandb_name=name,
                                  **self.mlff_parameters)

        # Temporarily replace sys.argv
        original_argv = sys.argv
        sys.argv = ['script_name']  # replace 'script_name' with the actual script name
        for key, value in vars(args).items():
            if value is not None:
                if type(value) == bool and value is True:
                    sys.argv.extend([f'--{key}'])
                else:
                    sys.argv.extend([f'--{key}', str(value)])

        # Call the main function
        mlff_run_dir = os.path.join(self.dir_mlff, f'{self.timestamp_fail}_{self.timestamp_train}')
        os.makedirs(mlff_run_dir, exist_ok=True)
        self.logger.info('Running MACE training in {}'.format(mlff_run_dir))
        os.chdir(mlff_run_dir)
        try:
            # Create a new logger
            new_logger = logging.getLogger('mace')
            new_logger.handlers = []
            # Replace the root logger
            logging.root = new_logger
            # Call
            # sys.exit()
            mace.cli.run_train.main()
            # Restore the original root logger
            mace.cli.run_train.main()
        except Exception as e:
            os.chdir(self.dir_initial_working)
            raise e

        os.chdir(self.dir_initial_working)
        # Restore the original sys.argv
        sys.argv = original_argv

    def retrain_mlff(self):
        """Retrain the MLFF with new data. """
        train_configs = self.read_train_configs()
        run_dir = os.path.join(self.train_dir, "mlff")
        self.timestamp_train = datetime.datetime.now().strftime(self.time_format)
        mace_name = self.mace_name_pat.format(self.timestamp_train, self.timestamp_fail)
        ...

    def read_train_configs(self):
        """Read atomic configurations used for training."""
        data_dir = os.path.join(self.train_dir, "data")
        return [
            read(os.path.join(data_dir, fname))
            for fname in os.listdir(data_dir)
            if fname.endswith(".xyz")
        ]
