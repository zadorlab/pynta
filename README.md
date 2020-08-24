# rmgcat_to_sella

<!-- ![workflow_idea](./workflow_idea.png) -->
The work-flow designed to automate search for transition states and reaction paths on various surfaces.

`rmgcat_to_sella` reads the .yaml files from the previous RMGCat calculations, puts reactants and products on the surface, calculates TSs and returns reaction energies as well as barriers heights represented as a free energy reaction path diagrams. The general idea of the code can be summarized in the following figure.
<center><img src='./workflow_idea.png' style="width:400px"></center>
<br>

## 1. How to install

1.1 Clone the project in your preferable location.

```
git clone https://gitlab-ex.sandia.gov/mgierad/rmgcat_to_sella.git
```
Usually, `master` branch should be fine. If somehow it is not working, make sure to switch to the latest stable version by checking the tags and looking for `stable`.

1.2. Go to `rmgcat_to_sella`
```
cd rmgcat_to_sella
```

1.3. Create a virtual environment:
```
virtualenv venv
```

1.4. Activate your virtual environment:
```
source venv/bin/activate
```

1.5. Install all dependencies and make sure they are fine. List of all dependencies will be here \[later\]:

1.6a Install `rmgcat_to_sella`:
```
python setup.py install
```

1.6b (optional) If you do not have admin privileges (e.g. you use it on a supercomputer), do the following instead of 1.6a:
```
python setup.py install --user
```

**You should be ready to go now!**

Once finished using the workflow:
```
cd rmgcat_to_sella
deactivate
```

## 2. How to run
### 2.1 Using Balsam

You will need **4** files to run the workflow:
- `run_me.py` a python script that executes the workflow
- `run_me.sh` a bash script that submits job to the SLURM scheduler
- `inputR2S.py` a python script holding all user-modifiable parameters 
- `reactions.yaml` a yaml file with all reactions

An example `run_me.py` file:

```python
#!/usr/bin/env python3
import time
import inputR2S
from rmgcat_to_sella.main import WorkFlow

workflow = WorkFlow()
workflow.gen_job_files()
workflow.execute()
```
An example `run_me.sh` file:

```bash
#!/bin/bash
#SBATCH -J job_name        # name of the job e.g job_name = rmgcat_to_sella_workflow
#SBATCH --partition=queue  # queue name e.g. queue = day-long-cpu
#SBATCH --nodes=x          # number of nodes e.g. x = 2
#SBATCH --ntasks=y         # number of CPUs e.g. 2 x 48 = y = 96
#SBATCH -e %x.err          # error file name
#SBATCH -o %x.out          # out file name

# load your quantum chemistry calculation package, e.g.
module load espresso/6.6_nostress
# activate balsam environment, e.g.
source balsamactivate ~/myWorkflow_bebop
# run python executable script
python3 $PWD/run_me.py
# required environment variable if using balsam branch serial-mode-perf
export SLURM_HOSTS=$(scontrol show hostname)
# launch serial jobs (required)
balsam launcher --job-mode=serial --wf-filter reactions.yaml --limit-nodes=1 --num-transition-threads=1 &
# give some time to prevent time out before the sockets are ready for the quantum chemistry application, e.g. pw.x for Quantum Espresso
sleep 100
# launch mpi jobs (required)
balsam launcher --job-mode=mpi --wf-filter QE_Sock --offset-nodes=1 --num-transition-threads=1 &
# wait until finished
wait
# deactivate balsam environment
source balsamdeactivate
```

An example `reactions.yaml` file:
```yaml
  - index: 0
    reaction: OHX + X <=> OX + HX
    reaction_family: Surface_Abstraction
    reactant: |
        multiplicity -187
        1 *1 O u0 p0 c0 {2,S} {4,S}
        2 *2 H u0 p0 c0 {1,S}
        3 *3 X u0 p0 c0 
        4    X u0 p0 c0 {1,S}
    product: |
        multiplicity -187
        1 *1 O u0 p0 c0 {4,S}
        2 *2 H u0 p0 c0 {3,S}
        3 *3 X u0 p0 c0 {2,S}
        4    X u0 p0 c0 {1,S}
    - index: 1
    reaction: H2OX + X <=> OHX + HX
    reaction_family: Surface_Abstraction
    reactant: |
        multiplicity -187
        1 *1 O u0 p0 c0 {2,S} {3,S} {4,S}
        2 *2 H u0 p0 c0 {1,S}
        3    H u0 p0 c0 {1,S}
        4    X u0 p0 c0 {1,S}
        5 *3 X u0 p0 c0 
    product: |
        multiplicity -187
        1 *1 O u0 p0 c0 {2,S} {4,S}
        2 *2 H u0 p0 c0 {1,S}
        3    H u0 p0 c0 {5,S}
        4    X u0 p0 c0 {1,S}
        5 *3 X u0 p0 c0 {3,S}
```


An example input files are also located at `./rmgcat_to_sella/example_run_files/`.

If you do not have a `.yaml` file with the reaction list but still want to use the work-flow, let me know. Also, stay tuned, as a version of `rmgcat_to_sella` that can work without `.yaml` file is currently under development


### 2.2 Using SLURM only
**Warning `dev` branch uses SLURM scheduler to deal with the job dependencies. Be aware that it might be a bit buggy and do not fully support all the features implemented in the `master` branch.**

An example script (using `dev` branch - SLURM):
```python
#!/usr/bin/env python3
#SBATCH -J job_name        # name of the job e.g job_name = rmgcat_to_sella_workflow
#SBATCH --partition=queue  # queue name e.g. queue = day-long-cpu
#SBATCH --nodes=x          # number of nodes e.g. x = 2
#SBATCH --ntasks=y         # number of CPUs e.g. 2 x 48 = y = 96
#SBATCH -e %x.err          # error file name
#SBATCH -o %x.out          # out file name

import os
import sys

# get environmental variable
submitDir = os.environ['SLURM_SUBMIT_DIR']
# change directory to $SLURM_SUBMIT_DIR
os.chdir(submitDir)
# add current working directory to the path
sys.path.append(os.getcwd())
# import input file with - can be done only after sys.path.append(os.getcwd())
import inputR2S
# import executable class of rmgcat_to_sella
from rmgcat_to_sella.main import WorkFlow
# instantiate the WorkFlow class
workflow = WorkFlow()
# generate input files
workflow.gen_job_files()
# execute the work-flow
workflow.execute()
```

## 3. Documentation

Documentation is currently under development.

## Acknowledgments

If you are using `rmgcat_to_sella` or you wish to use it, let me know!

This work was supported by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, Chemical Sciences, Geosciences and Biosciences Division, as part of the Computational Chemistry Sciences Program.


