#!/bin/bash
#SBATCH -J pynta_run_04
#SBATCH -p plgrid-short
#SBATCH -A plgpynta
#SBATCH --time=0:20:00
#SBATCH -N 3
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=5GB
#SBATCH -e %x.err
#SBATCH -o %x.out


export SLURM_HOSTS=$(scontrol show hostname)
source /net/scratch/people/plgveronalx/01_venvs/pynta/bin/activate
module load plgrid/apps/espresso/6.4.1-impi
yes | balsam init /net/scratch/people/plgveronalx/balsamDBs/run_04
yes | source balsamactivate /net/scratch/people/plgveronalx/balsamDBs/run_04

source balsamactivate /net/scratch/people/plgveronalx/balsamDBs/run_04
echo "executing python"
python3 $PWD/run_me.py
echo "Done"
echo "Balsam submit-launch"
#balsam submit-launch -A plgpynta -q plgrid-short -t 1200 -n 1 --job-mode=mpi --sched-flags="-n 24"
echo "Done"
# python3 $PWD/restart_me.py
echo "executing first balsam launcher"
yes | balsam launcher --job-mode=serial --limit-nodes=1 --num-transition-threads=1 --consume-all &
echo "Done"
sleep 60
echo "executing second balsam launcher"
yes | balsam launcher --job-mode=mpi --offset-nodes=1 --num-transition-threads=1 --consume-all &
echo "Done"
wait
source balsamdeactivate
