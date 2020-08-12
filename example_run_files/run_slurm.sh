#!/bin/bash
#SBATCH --nodes=2
source balsamactivate ~/myWorkflow_bebop
#python3 $PWD/run_me.py 
export SLURM_HOSTS=$(scontrol show hostname)
balsam launcher --job-mode=serial --wf-filter reactions.yaml --limit-nodes=1 --num-transition-threads=1 &
sleep 15
balsam launcher --job-mode=mpi --wf-filter QE_Sock --offset-nodes=1 --num-transition-threads=1 &
wait
source balsamdeactivate
