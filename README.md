# rmgcat_to_sella

A work-flow that reads the .yaml files from the previous RMGCat calculations, puts reactants and products on the surface, calculates TSs and gives reaction energies as well as barriers heights represented as a free energy reaction path diagrams.

## How to install

Make sure all dependencies are correctly installed. List of all dependencies will be here:

Clone the project in your preferable location. Make sure you switch to the latest stable version by checking the tags and looking for `stable`. Be aware that if you decided to use a `dev` branch, it might be a bit buggy.

On your machine go to `rmgcat_to_sella`. Run the following:
```
python setup.py install
```
If you do not have admin priveleges (e.g. you use it on a supercomputer), do the following instead:
```
python setup.py install --user
```
## How to run
An example script:
```python
#!/usr/bin/env python3
# your SBATCH options goes here
#SBATCH ...
import os
import sys
submitDir = os.environ['SLURM_SUBMIT_DIR']
os.chdir(submitDir)
sys.path.append(os.getcwd())
import inputR2S
from rmgcat_to_sella.main import WorkFlow

workflow = WorkFlow()
workflow.gen_job_files()
workflow.execute()
```
An example input files is located at `/example_run_files/inputR2S.py`. You will also need `reactions.yaml` which can be found in the same location. 

If do not hava a `.yaml` file with the reaction list but still want to use the work-flow, let me know. Also, stay tuned, as a version of rmgcat_to_sella that can work without `.yaml` file is currently under development

If you are using rmgcat_to_sella or you wish to use it, let me know!

## Documentation

Documentation is currently under development.

## Acknowledgments

This work was supported by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, Chemical Sciences, Geosciences and Biosciences Division, as part of the Computational Chemistry Sciences Program.


