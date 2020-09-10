#!/bin/bash
#SBATCH -J test_workflow
#SBATCH --partition=day-long-cpu
#SBATCH --nodes=2
#SBATCH --ntasks=96
#SBATCH -e %x.err
#SBATCH -o %x.out

# load your quantum chemistry calculation package, e.g.
module load espresso/6.6_nostress
# activate balsam environment, e.g.
source balsamactivate ~/myWorkflow_bebop
# run python executable script
python3 $PWD/run_me.py
# required environment variable if using balsam branch serial-mode-perf
export SLURM_HOSTS=$(scontrol show hostname)
# launch serial jobs (required)
balsam launcher --job-mode=serial --wf-filter _ --limit-nodes=1 --num-transition-threads=1 &
# give some time to prevent time out before the sockets are ready for the quantum chemistry application, e.g. pw.x for Quantum Espresso
sleep 45
# launch mpi jobs (required)
balsam launcher --job-mode=mpi --wf-filter QE_Sock --offset-nodes=1 --num-transition-threads=1 &
# wait until finished
wait
# deactivate balsam environment
source balsamdeactivate
