#!/bin/bash -l
#SBATCH --job-name=prep
#SBATCH --partition=month-long-cpu
#SBATCH --account=shikim
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=30-00:00:00
#SBATCH --output=pynta_output.out
#SBATCH --error=pynta_output.error
#SBATCH --exclude=node7

source /home/shikim/anaconda3/etc/profile.d/conda.sh
conda activate pynta_env

python runprep.py