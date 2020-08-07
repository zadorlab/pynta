#!/bin/bash

# module load quantum-espresso/6.4.1-bdw
#balsam init ~/myWorkflow
source balsamactivate /projects/catalysis_aesp/brossdh/myWorkflow
# export BALSAM_DB_PATH='~/myWorkflow_bebop:$BALSAM_DB_PATH'
# balsam job --name launch_job --workflow test --application Python3 --args $PWD/run_me.py --yes
python3 $PWD/run_me.py
balsam submit-launch -A catalysis_aesp -q debug-cache-quad  -t 60 -n 2 --job-mode=mpi
wait
#workflow 'QE_Socket'
# balsam submit-launch -A CONDO -q day-long-cpu -t 1200 -n 1 --job-mode=serial
#balsam submit-launch -A whatever -q day-long-cpu -t 1200 -n 1 --job-mode=serial --sched-flags="-n 48"
#--job-mode=mpi
#watch balsam ls   #  follow status in realtime from co
