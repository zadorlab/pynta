#!/bin/bash
#SBATCH -J rmgcat_to_sella_test_run
#SBATCH --partition=week-long-cpu
#SBATCH --nodes=6
#SBATCH --ntasks=287
#SBATCH -e %x.err
#SBATCH -o %x.out

source balsamactivate ~/myWorkflow_bebop
python3 $PWD/run_me.py
# python3 $PWD/restart_me.py
export SLURM_HOSTS=$(scontrol show hostname)
balsam launcher --job-mode=serial --wf-filter _ --limit-nodes=1 --num-transition-threads=1 &
sleep 60
balsam launcher --job-mode=mpi --wf-filter QE_Sock --offset-nodes=1 --num-transition-threads=1 &
wait
source balsamdeactivate

