# rmgcat_to_sella

Still a quite buggy work-flow that reads the .yaml files from the previous RMGCat calculations, puts reactants and products on the surface, calculates TSs and gives reaction energies as well as barriers heights.

An example script
```python
#!/usr/bin/env python3
# your SBATCH options goes here
#SBATCH ...

import os
import sys
submitDir = os.environ['SLURM_SUBMIT_DIR']
os.chdir(submitDir)
sys.path.append(os.getcwd())
import time
import inputR2S
from rmgcat_to_sella.main import WorkFlow

workflow = WorkFlow()
workflow.gen_job_files()
workflow.execute()
```
An example input files is located at `/example_run_files/inputR2S.py`

If you are using rmgcat_to_sella or you wish to use it, let me know!

## Documentation

Documentation is currently under development.

## Acknowledgments

This work was supported by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, Chemical Sciences, Geosciences and Biosciences Division, as part of the Computational Chemistry Sciences Program.


