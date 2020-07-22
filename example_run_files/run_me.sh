#!/bin/sh

module load quantum-espresso/6.4.1-bdw
source balsamactivate ~/myWorkflow_bebop
#balsam job --name launch_job --workflow test --application Python --args $PWD/run_me.py --yes
python $PWD/run_me.py
balsam submit-launch -A CONDO -q csed -t 5 -n 1 --job-mode=mpi
#watch balsam ls   #  follow status in realtime from co
