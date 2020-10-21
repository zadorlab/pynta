#!/bin/bash
#COBALT --jobname rmgcat_to_sella_test
#COBALT -A catalysis_aesp
#COBALT -n 141
#COBALT -t 2:00:00

source balsamactivate /projects/catalysis_aesp/mgierad/myWorkflow
#conda activate rmgcat_to_sella
module unload cray-libsci
module load cray-hdf5-parallel
python3 $PWD/run_me.py
balsam launcher --job-mode=serial --wf-filter _ --limit-nodes=1 --num-transition-threads=1 &
sleep 45
balsam launcher --job-mode=mpi --wf-filter QE_Sock --offset-nodes=1 --num-transition-threads=1 &
wait
#conda deactivate
source balsamdeactivate
