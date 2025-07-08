import os
import json
# ase general
from ase.io import read, write, Trajectory
from ase.build import add_adsorbate, bulk
from ase.constraints import FixAtoms
# visulization
from ase.visualize import view
import numpy as np
import matplotlib.pyplot as plt
# call pynta module
from pynta.utils import name_to_ase_software
from pynta.tasks import *
from pynta.calculator import get_lattice_parameters
from pynta.main import generate_slab

# fireworks
from fireworks import LaunchPad, Workflow
from fireworks.queue.queue_launcher import rapidfire as rapidfirequeue
from fireworks.utilities.fw_serializers import load_object_from_file
from fireworks.core.fworker import FWorker


def plot_results(cutoff_energy, kpts_list, energies):
    """
    Plot the optimized energy against the kinetic energy cutoff and k-points.

    Parameters:
    - results: dict, a dictionary with ecut values as keys and corresponding optimized energies as values.
    - kpts_list: list, the list of k-point grids tested.
    - energies: list, the corresponding energies for each k-point grid.
    """

    def plot_energy_vs_cutoff():
        """
        Plot optimized energy against kinetic energy cutoff and save the plot.
        """
        # Extract ecut values and corresponding energies
        ecut_values = list(cutoff_energy.keys())
        optimized_energies = list(cutoff_energy.values())

        # Create the plot
        plt.figure(figsize=(10, 6))
        plt.plot(ecut_values, optimized_energies, marker='o', linestyle='-', color='b')

        # Add labels and title
        plt.xlabel('Kinetic Energy Cutoff (eV)', fontsize=14)
        plt.ylabel('Optimized Energy (eV)', fontsize=14)
        plt.title('Optimized Energy vs Kinetic Energy Cutoff', fontsize=16)
        plt.grid(True)

        # Save the plot as a PNG file
        plt.tight_layout()
        plt.savefig('optimized_energy_vs_cutoff.png')
        plt.close()  # Close the plot to free memory

    def plot_energy_vs_kpoints():
        """
        Plot optimized energy against k-points and save the plot.
        """
        # Convert kpts_list to a suitable format for plotting
        kpts_values = [kpts[0] * kpts[1] * kpts[2] for kpts in kpts_list]  # Calculate total k-points

        # Create the plot
        plt.figure(figsize=(10, 6))
        plt.plot(kpts_values, energies, marker='o', linestyle='-', color='r')

        # Add labels and title
        plt.xlabel('K-points', fontsize=14)
        plt.ylabel('Optimized Energy (eV)', fontsize=14)
        plt.title('Optimized Energy vs K-points', fontsize=16)
        plt.grid(True)

        # Save the plot as a PNG file
        plt.tight_layout()
        plt.savefig('optimized_energy_vs_kpoints.png')
        plt.close()  # Close the plot to free memory

    # Call the sub-functions to create the plots
    plot_energy_vs_cutoff()
    plot_energy_vs_kpoints()

# create Prep class for running calculations
"""
Steps for preprocessing follows:
1. Initialize Prep class
2. Compute optimal lattice constant for the slab with given metal 
3. Generate slab with optimal lattice constant
4. Add an adsorbate on the slab
5. Run energy convergence calculations using gamma point
6. With energy cutoff range, do k-point convergence test
7. Evaluate sites for custom surface (Print out the output in json file)
"""
class Prep:
    def __init__(self, metal='Pt', surface_type='fcc111', a0=3.96, software='Espresso', fmax=0.05, adsorbate=None, position=None, vacuum=10, slab=None,
                 launchpad_path=None,fworker_path=None,reset_launchpad=False,queue_adapter_path=None,queue=False,njobs_queue=0,num_jobs=25, ecut_range=None,
                 software_kwargs={'kpts': (3, 3, 1), 'tprnfor': True, 'occupations': 'smearing',
                            'smearing':  'marzari-vanderbilt',
                            'degauss': 0.01, 'ecutwfc': 40, 'nosym': True,
                            'conv_thr': 1e-6, 'mixing_mode': 'local-TF',
                            "pseudopotentials": 
                            {"Cu": 'Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',"H": 'H.pbe-kjpaw_psl.1.0.0.UPF',"O": 'O.pbe-n-kjpaw_psl.1.0.0.UPF',"C": 'C.pbe-n-kjpaw_psl.1.0.0.UPF',"N": 'N.pbe-n-kjpaw_psl.1.0.0.UPF',},
                              },):
        """
        Initialize the Prep class with metal, surface type, lattice constant, software, and software parameters.

        Parameters:
        - metal: str, the type of metal (default is Cu').
        - surface_type: str, the type of surface (default is 'fcc111').
        - a0: float, the lattice constant (default is 4.04 Å).
        - software: str, the software to use for calculations (default is 'Espresso').
        - software_kwargs: dict, additional parameters for the software calculator (default is None).
        - fmax : float, the maximum force threshold for convergence (default is 0.05 eV/Å).
        - adsorbate : str, adsorbate to place on ase supported/built slab
        - position : str, positions supported on given slab. Top, bridge, fcc, hcp
        - vacuum : float, additional parameter for building slab
        - slab : str, slab file path. /path/to/your/slab.xyz

        """
        # Fireworks keywords
        self.qadapter = None
        self.launchpad = launchpad  
        if launchpad_path:
            launchpad = LaunchPad.from_file(launchpad_path)
        else:
            launchpad = LaunchPad()
        if reset_launchpad:
            launchpad.reset('', require_password=False)

        if fworker_path:
            self.fworker = FWorker.from_file(fworker_path)
        else:
            self.fworker = FWorker()

        self.queue = queue
        if queue:
            self.qadapter = load_object_from_file(queue_adapter_path)

        self.njobs_queue = njobs_queue
        self.num_jobs = num_jobs

        # DFT, ASE, slab calculation related keywords
        self.metal = metal
        self.surface_type = surface_type
        self.a0 = a0
        self.software = software  # Add software as an attribute
        self.software_kwargs = software_kwargs if software_kwargs is not None else {}  # Initialize software_kwargs
        self.fmax = fmax
        self.adsorbate = adsorbate
        self.position = position
        self.vacuum = vacuum
        self.slab = slab
        self.ecut_range = ecut_range
    

    def opt_cutoff_energy(self, ecut_range=None, skip_launch=False):
        """
        Optimize bulk aluminum by varying the kinetic energy cutoff.

        Parameters:
        - slab: ASE Atoms object or None, the slab to optimize (default is None).
        - ecut_range: tuple, the range of kinetic energy cutoffs to test.
        - fmax: float, the maximum force threshold for convergence (default is 0.05 eV/Å).

        Returns:
        - cutoff_energy: list, a list of tuples with (ecut, optimized_energy).
        """
        soft = name_to_ase_software(self.software)(**self.software_kwargs)  # Initialize the calculator

        # Use provided slab or create a bulk structure of given metal
        if self.slab is None:
            self.slab = bulk(self.metal, self.surface_type[:3], a=self.a0)
            self.slab.calc = soft
            self.slab.pbc = (True, True, False)
        
        # Construct bulk slab for convergence calculation
        write(os.path.join(self.path,"test_bulk.xyz"),slab)

        cutoff_energy = {}  # Initialize results list

        # Loop over kinetic energy cutoffs
        for ecut in range(*ecut_range):
            # Update configuration based on the selected calculator
            if self.software == 'gpaw':
                self.software_kwargs['ecut'] = ecut  # Update ecut in configuration

            elif self.software == 'espresso':
                self.software_kwargs['ecutwfc'] = ecut  # Update ecutwfc in configuration

            elif self.software in ['nwchem', 'pwdft']:
                self.software_kwargs['cutoff'] = ecut  # Update cutoff in configuration

            # Initialize the optimizer : commented out since we may not need to optimize the bulk structure
            #optimizer = BFGS(slab)
            # Optimize the structure
            #optimizer.run(fmax=fmax)  # Stop when forces are below fmax

#energy_firework(xyz,software,label,software_kwargs={},parents=[],out_path=None,ignore_errors=False):

            # If skip_launch is true, do not use Fireworks. Calculations will be done locally.
            if skip_launch:
                optimized_energy = self.slab.get_potential_energy()
                cutoff_energy.append((ecut, optimized_energy))  # Append a tuple of (ecut, optimized_energy)
                print(f'Calculator: {self.software}, Cutoff: {ecut} eV, Optimized Energy: {optimized_energy:.6f} eV')

            else:  # Get the optimized energy via Fireworks
                fwenergy = energy_firework(os.path.join(self.path, "slab_bulk.xyz"), self.software, "energy_convergence_test", 
                                           software_kwargs=self.software_kwargs, out_path=os.path.join(self.path, "energy_convergence_test_energy.json"))
                wfenergy = Workflow([fwenergy], name=self.label)
                self.launchpad.add_wf(wfenergy)  # Only add once
                # Wait for the Firework to complete and retrieve the results
                while not os.path.exists(os.path.join(self.path, "energy_convergence_test_energy.json")):
                    time.sleep(1)  # Wait for the output file to be created
                # Read the optimized energy from the output file
                with open(os.path.join(self.path, "energy_convergence_test_energy.json"), 'r') as file:
                    optimized_energy = json.load(file)

                cutoff_energy.append((ecut, optimized_energy))  # Append the results
                print(f'Calculator: {self.software}, Cutoff: {ecut} eV, Optimized Energy: {optimized_energy:.6f} eV')

        # Save results to a JSON file
        with open('cutoff_energy.json', 'w') as f:
            json.dump(cutoff_energy, f, indent=4)  # Write the cutoff_energy list to the JSON file

        return cutoff_energy

    def opt_kpoints(self, slab=None, kpts_range=None, skip_launch=False):
        """
        Optimize k-points for a slab structure and save the results.

        Parameters:
        - slab: ASE Atoms object or None, the slab to optimize (default is None).
        - kpts_range: list or array, the range of k-point grids to test.
        - output_file: str, the name of the output file for logging (default is 'kpts.json').
        - size: tuple, the size of the slab.
        - vacuum: float, the vacuum thickness to add to the slab.
        - skip_launch: bool, whether to skip launching Fireworks and calculate locally.

        Returns:
        - kpts_list: list, the list of k-point grids tested.
        - energies: list, the corresponding energies for each k-point grid.
        """

        # Initialize the calculator
        soft = name_to_ase_software(self.software)(**self.software_kwargs)

        energies = []
        kpts_list = []

        # Check if slab is provided; if not, create it using make_test_slab
        if slab is None:
            slab = self.make_test_slab(size=self.size, vacuum=self.vacuum, adsorbate=self.adsorbate, position=self.position)  # Call the make_test_slab function
            write(os.path.join(self.path, "test_slab.xyz"), slab)  # Save the slab to a file

        # Dictionary to store results
        results = {}

        # Open the output file for writing results
        for kpts in kpts_range:
            self.conf_qe['kpts'] = kpts  # Update k-points in configuration

            # If skip_launch is true, calculate locally
            if skip_launch:
                slab.calc = soft  # Assign the calculator
                slab.pbc = (True, True, False)  # Set periodic boundary conditions

                # Calculate the potential energy
                eslab = slab.get_potential_energy()
                kpts_list.append(kpts)  # Store k-points
                energies.append(eslab)  # Store corresponding energy

                # Store the result in the dictionary with tuple as string key
                results[''.join(map(str, kpts))] = eslab  # Use concatenated string representation of the tuple

            else:  # Get the optimized energy via Fireworks
                kpts_name = ''.join(map(str, kpts))  # Create a name like '111', '221', etc.
                fwkpt = energy_firework(os.path.join(self.path, "test_slab.xyz"), self.software, f"kpoints_{kpts_name}", software_kwargs=self.software_kwargs)
                wfkpt = Workflow([fwkpt], name=self.label)
                self.launchpad.add_wf(wfkpt)  # Add the workflow to the launchpad

                # Wait for the Firework to complete and retrieve the results
                while not os.path.exists(os.path.join(self.path, f"kpoints_{kpts}_energy.json")):
                    time.sleep(1)  # Wait for the output file to be created

                # Read the optimized energy from the output file
                with open(os.path.join(self.path, f"kpoints_{kpts}_energy.json"), 'r') as file:
                    eslab = json.load(file)

                kpts_list.append(kpts)  # Store k-points
                energies.append(eslab)  # Store corresponding energy

                # Store the result in the dictionary with tuple as string key
                results[str(kpts)] = eslab  # Use string representation of the tuple

        # Write the results to the output JSON file
        with open('kpts.json', 'w') as f:
            json.dump(results, f, indent=4)  # Write the results dictionary to the JSON file

        return kpts_list, energies

    def generate_slab(self):
        """
        Create a slab structure for a given metal facet and size, and apply constraints.

        Parameters:
        - size: tuple, the size of the slab (e.g., (2, 2, 3)).
        - vacuum: float, the amount of vacuum to add (default is 10.0 Å).
        - adsorbate: str or None, the type of adsorbate to add (default is None).
        - position: str, ASE have pre-determined positions. 'top', 'bridge', 'fcc', 'hcp' , 'shortbridge', 'longbridge'
        - filename: str or None, the name of the file to save the slab (default is None).

        Returns:
        - slab: ASE Atoms object representing the slab.
        """

        slab_type = getattr(ase.build,self.surface_type)
        #optimize the lattice constant
        if self.a is None:
            a = get_lattice_parameters(self.metal,self.surface_type,self.software,self.lattice_opt_software_kwargs)
            print("computed lattice constants of: {} Angstroms".format(a))
            if isinstance(a,float):
                self.a = a
            else:
                self.a = a[0]
                self.c = a[1]
        
        logger.info('Construct slab with optimal lattice constant')
        #construct slab with optimial lattice constant
        if self.c:
            slab = slab_type(symbol=self.metal,size=self.repeats,a=self.a,vacuum=self.vacuum,c=self.c)
        else:
            slab = slab_type(symbol=self.metal,size=self.repeats,a=self.a,vacuum=self.vacuum)
        slab.pbc = self.pbc
        write(os.path.join(self.path,"slab_init.xyz"),slab)
        self.slab_path = os.path.join(self.path,"slab.xyz")
        if self.software != "XTB":
            fwslab = optimize_firework(os.path.join(self.path,"slab_init.xyz"),self.software,"slab",
                opt_method="BFGSLineSearch",socket=self.socket,software_kwargs=self.software_kwargs,
                run_kwargs={"fmax" : self.fmaxopt},out_path=os.path.join(self.path,"slab.xyz"),constraints=["freeze up to {}".format(self.freeze_ind)],priority=1000)
            wfslab = Workflow([fwslab], name=self.label+"_slab")
            self.launchpad.add_wf(wfslab)
            if skip_launch:
                return
            self.launch(single_job=True)
            while not os.path.exists(self.slab_path): #wait until slab optimizes, this is required anyway and makes the rest of the code simpler
                time.sleep(1)
            self.slab = read(self.slab_path)
        else: #testing
            self.slab = slab
            write(self.slab_path,slab)

        return test_slab
