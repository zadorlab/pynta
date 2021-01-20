#!/bin/bash
#COBALT --jobname pynta_run
#COBALT -A <project_id>
#COBALT -n 128
#COBALT -t 3:00:00

source balsamactivate ~/Workflow
module unload cray-libsci
module load cray-hdf5-parallel
python3 $PWD/run_me.py
balsam launcher --job-mode=serial --wf-filter _ --limit-nodes=1 --num-transition-threads=1 &
sleep 45
balsam launcher --job-mode=mpi --wf-filter QE_Sock --offset-nodes=1 --num-transition-threads=1 &
wait
source balsamdeactivate
