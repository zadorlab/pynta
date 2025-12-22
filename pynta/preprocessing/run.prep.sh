#!/bin/bash -l

#SBATCH --time=30-00:00:00
#SBATCH --partition=month-long-cpu
#SBATCH --account=shikim
#SBATCH --job-name=prep
#SBATCH --output=pynta_output.out
#SBATCH --error=pynta_output.error
#SBATCH --exclude=node7

python runprep.py

# CommonAdapter (SLURM) completed writing Template
