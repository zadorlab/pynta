#!/bin/bash
#COBALT --jobname rmgcat_to_sella
#COBALT -A catalysis_aesp
#COBALT -n 8
#COBALT -t 0:58:00
#COBALT -q debug-cache-quad

source balsamactivate /projects/catalysis_aesp/brossdh/myWorkflow
module unload cray-libsci
module load cray-hdf5-parallel
python $PWD/run_me.py 
balsam launcher --job-mode=serial --wf-filter reactions.yaml --limit-nodes=1 --num-transition-threads=1 &
sleep 45
balsam launcher --job-mode=mpi --wf-filter QE_Sock --offset-nodes=1 --num-transition-threads=1 &
wait
