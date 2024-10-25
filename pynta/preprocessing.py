import os
# ase general
from ase.io import read, write, Trajectory
from ase.calculators import calculator
from ase.build import add_adsorbate, bulk
from ase.constraints import FixAtoms
# visulization
from ase.visualize import view
import numpy as np
import matplotlib.pyplot as plt
# call pynta module
from pynta.utils import name_to_ase_software

# Analysis of the results: this does not require making prep class
def analyze_lattice_parameter(a0, c0, optimized_a, surface_type):
    if surface_type != "hcp0001":
        # Check if a0 is equal or greater than optimized_a by 0.05
        if a0 >= optimized_a + 0.05:
            print("Your lattice constant may not be optimal for Pynta calculation. Please generate slab using pynta. You can do it by pyn.generate_slab()")
    else:
        # Check if a0/c0 is equal or greater than optimized_a/a0 by 0.05
        if (a0 / c0) >= (optimized_a / c0) + 0.05:
            print("Your c and a may not be optimal for Pynta calculation. Please generate slab using pynta. You can do it by pyn.generate_slab()")


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
    #plot_energy_vs_cutoff()
    #plot_energy_vs_kpoints()

# create Prep class for running calculations
class Prep:
    def __init__(self, metal='Pt', surface_type='fcc111', a0=3.96, software='Espresso', fmax=0.05,
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
        - metal: str, the type of metal (default is 'Al').
        - surface_type: str, the type of surface (default is 'fcc111').
        - a0: float, the lattice constant (default is 4.04 Å).
        - software: str, the software to use for calculations (default is 'Espresso').
        - software_kwargs: dict, additional parameters for the software calculator (default is None).

        """
        self.metal = metal
        self.surface_type = surface_type
        self.a0 = a0
        self.software = software  # Add software as an attribute
        self.software_kwargs = software_kwargs if software_kwargs is not None else {}  # Initialize software_kwargs
        self.fmax = fmax

    def opt_cutoff_energy(self, slab=None, ecut_range=(200, 501, 50), output_file='cutoff_energy.txt'):
        """
        Optimize bulk aluminum by varying the kinetic energy cutoff.

        Parameters:
        - slab: ASE Atoms object or None, the slab to optimize (default is None).
        - ecut_range: tuple, the range of kinetic energy cutoffs to test.
        - fmax: float, the maximum force threshold for convergence (default is 0.05 eV/Å).
        - output_file: str, the name of the file to save results (default is 'cutoff_energy.txt').

        Returns:
        - cutoff_energy: list, a list of tuples with (ecut, optimized_energy).
        """
        soft = name_to_ase_software(self.software)(**self.software_kwargs)  # Initialize the calculator

        # Use provided slab or create a bulk structure of aluminum
        if slab is None:
            slab = bulk(self.metal, self.surface_type[:3], a=self.a0)
            slab.calc = soft
            slab.pbc = (True, True, True)

        cutoff_energy = []  # Initialize results list

        # Loop over kinetic energy cutoffs
        for ecut in range(*ecut_range):
            # Update configuration based on the selected calculator
            if self.software == 'gpaw':
                self.software_kwargs['ecut'] = ecut  # Update ecut in configuration

            elif self.software == 'espresso':
                self.software_kwargs['ecutwfc'] = ecut  # Update ecutwfc in configuration

            elif self.software in ['nwchem', 'pwdft']:
                self.software_kwargs['cutoff'] = ecut  # Update cutoff in configuration

            # Initialize the optimizer
            #optimizer = BFGS(slab)
            # Optimize the structure
            #optimizer.run(fmax=fmax)  # Stop when forces are below fmax

            # Get the optimized energy
            optimized_energy = slab.get_potential_energy()
            cutoff_energy.append((ecut, optimized_energy))  # Append a tuple of (ecut, optimized_energy)
            print(f'Calculator: {self.software}, Cutoff: {ecut} eV, Optimized Energy: {optimized_energy:.6f} eV')

        # Save results to a file
        with open(output_file, 'w') as f:
            f.write("# Cutoff Energy (eV), Optimized Energy (eV)\n")
            for ecut, energy in cutoff_energy:
                f.write(f"{ecut}, {energy:.6f}\n")

        return cutoff_energy

    def opt_kpoints(self, slab=None, kpts_range=None, output_file='kpts.txt', size=(2, 2, 3), vacuum=10):
        """
        Optimize k-points for a slab structure and save the results.

        Parameters:
        - slab: ASE Atoms object or None, the slab to optimize (default is None).
        - kpts_range: list or array, the range of k-point grids to test.
        - output_file: str, the name of the output file for logging (default is 'kpts.txt').

        Returns:
        - kpts_list: list, the list of k-point grids tested.
        - energies: list, the corresponding energies for each k-point grid.
        """
        # Get the slab type based on surface_type
        slab_type = getattr(ase.build, self.surface_type)
        # Initialize the calculator
        soft = name_to_ase_software(self.software)(**self.software_kwargs)

        energies = []
        kpts_list = []

        # Check if slab is provided; if not, create it
        if slab is None:
            slab = slab_type(self.metal, size=size)  
            slab.center(vacuum=vacuum.0, axis=2)  # Center the slab with vacuum

        # Open the output file for writing results
        with open(output_file, "w") as file:
            file.write(f"# Kpoints  Energy\n")  # Write header

            for kpts in kpts_range:
                self.conf_qe['kpts'] = kpts  # Update k-points in configuration
                
                slab.calc = soft  # Assign the calculator
                slab.pbc = (True, True, False)  # Set periodic boundary conditions

                # Calculate the potential energy
                eslab = slab.get_potential_energy()
                kpts_list.append(kpts)  # Store k-points
                energies.append(eslab)  # Store corresponding energy
                
                # Write the k-point and energy to the file
                file.write(f"{kpts}  {eslab:.6f}\n")

        return kpts_list, energies

    def make_test_slab(self, size=(2, 2, 3), vacuum=10.0, adsorbate=None, position=None):
        """
        Create a slab structure for a given metal facet and size, and apply constraints.

        Parameters:
        - size: tuple, the size of the slab (e.g., (2, 2, 3)).
        - vacuum: float, the amount of vacuum to add (default is 10.0 Å).
        - adsorbate: str or None, the type of adsorbate to add (default is None).
        - position: str, ASE have pre-determined positions. 'top', 'bridge', 'fcc', 'hcp' , 'shortbridge', 'longbridge'

        Returns:
        - slab: ASE Atoms object representing the slab.
        """
        # Get the slab type based on surface_type
        slab_type = getattr(ase.build, self.surface_type)
        # Initialize the calculator

        # Create the slab structure
        slab = slab_type(self.metal, size=size)
        slab.center(vacuum=vacuum, axis=2)

        # Create a mask to fix the bottom two layers of the slab
        mask_layers = [atom.tag > 2 for atom in slab]
        fix_bottom_two = FixAtoms(mask=mask_layers)
        slab.set_constraint(fix_bottom_two)

        # Optionally add an adsorbate
        if adsorbate is not None:
            add_adsorbate(slab, adsorbate, position=position)

        # Visualize the slab
        view(slab)

        return slab
