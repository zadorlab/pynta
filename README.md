# Pynta

<!-- ![workflow_idea](./workflow_idea.png) -->

The work-flow designed to automate search for transition states and reaction paths on various surfaces.

`pynta` reads the .yaml files from the previous RMGCat calculations, puts reactants and products on the surface, calculates TSs and returns reaction energies as well as barriers heights represented as a free energy reaction path diagrams. The general idea of the code can be summarized in the following figure.

<center><img src='./workflow_idea.png' style="width:400px"></center>
<br>

# 1. How to install

**Instalation can be a bit tricky because `pynta` depends on a several different packages that are under constant development (`balsam`, `xtb-python`)**

The following instruction assumes that you have several softwares installed on your system, such as:

- `Slurm`
- `OpenMPI`
- `OpenBlas`
- `GCC`

## 1.1 Install all prerequisites

1.1.1 Chooose the location you want to build `pynta`. It is suggested to create a virtual environment `pynta` for the `pynta` builds:

```bash
cd <path/whr/to/build/pynta>
python3 -m venv pynta
```

1.1.2 Activate your virtual environment:

```bash
source pynta/bin/activate
```

1.1.3 (_Optional_) Install required python packages (If you skip this process, `pynta` installer will do this later.)

```bash
pip3 install matplotlib<3.2 spglib==1.14.1.post0 networkx<2.4 ase==3.19.0 scipy==1.3.1 numpy==1.18.1 PyYAML==5.3.1 sella==1.0.3

```

1.1.4 Download [`PostgreSQL`](https://www.enterprisedb.com/download-postgresql-binaries) precompiled binaries that suits your system and add `path_to_PostgreSQL/pgsql/bin` to your `PATH` by modifying `~/.bashrc`

```bash
echo 'export PATH=path_to_PostgreSQL/pgsql/bin:$PATH' >> ~/.bashrc
```

1.1.5 Install [`mpi4py`](https://github.com/mpi4py/mpi4py.git):

```bash
git clone https://github.com/mpi4py/mpi4py.git
cd mpi4py
python3 setup.py install --user
cd ../
```

Make sure it works by running

```bash
srun -n 2 python3 -c 'from mpi4py import MPI; print(MPI.COMM_WORLD.Get_rank())'
```

Which should give

```
0
1
```

1.1.6 Install [`balsam`](https://github.com/balsam-alcf/balsam.git) using [`serial-mode-perf`](https://github.com/balsam-alcf/balsam/tree/serial-mode-perf) branch.

```bash
git clone https://github.com/balsam-alcf/balsam.git -b serial-mode-perf
cd balsam
python3 setup.py install --user
cd ../
```

Make sure it works by running tests posted on the `balsam` GitHub page.

1.1.7 Install [`xTB-python`](https://github.com/grimme-lab/xtb-python) following instruction provided there. Make sure to correctly link all required libraries. For example:

- using `OpenBlas` and `GNU` based compilers:

```bash
git clone https://github.com/grimme-lab/xtb-python.git
cd xtb-python
git submodule update --init
LDFLAGS="-L/opt/custom/OpenBLAS/0.3.7/lib" meson setup build --prefix=$PWD --libdir=xtb/xtb --buildtype release --optimization 2 -Dla_backend=openblas
ninja -C build install
pip install --user -e .
```

- using `MKL` and Intel Compilers:

```bash
git clone https://github.com/grimme-lab/xtb-python.git
cd xtb-python
git submodule update --init
# (Theta specific)
# conda instal cffi
# module swap PrgEnv-intel PrgEnv-cray; module swap PrgEnv-cray PrgEnv-intel
CC=icc CXX=icpc FC=ifort meson setup build --prefix=$PWD --libdir=xtb -Dla_backed=mkl -Dpy=3 --buildtype release --optimization 2
ninja -C build install
pip install --user -e .
```

Make sure it works by running:

```python
from ase.build import molecule
from xtb.ase.calculator import XTB
atoms = molecule('H2O')
atoms.calc = XTB(method="GFN2-xTB")
total_ener = atoms.get_potential_energy()
print(total_ener)
```

The expected output should be something around `-137.9677758730299`

**Warning - You might be getting SEGFAULT error -**

`Segmentation Fault (Core dumped)`

**while executing any** `xTB-python` **job, especially for a relatively large molecules. The easiest solution is to unlimit the system stack to avoid stack overflows. In** `bash` **try:**

```
ulimit -s unlimited
```

If `xTB-python` still fails, try to install [`xtb`](https://github.com/grimme-lab/xtb) and test `xTB` itself for any errors.

```bash
git clone https://github.com/grimme-lab/xtb.git
cd xtb
mkdir build
pushd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_FC_COMPILER=ifort ..
make
ctest
popd
echo 'export LD_LIBRARY_PATH=path/to_xtb/xtb/build:$LD_LIBRARY_PATH' >> ~/.bashrc
echo 'export PATH=$HOME/.local/bin:\$PATH' >> ~/.bashrc
```

Then rebuild `xTB-python` on your system ignoring `git submodule update --init` and linking you current `xTB` installation.

## 1.2 Install `pynta`

1.2.1 Clone the project in your preferable location.

```
git clone https://gitlab-ex.sandia.gov/mgierad/pynta.git
```

Usually, `master` branch should be fine. If somehow it is not working, make sure to switch to the latest stable version by checking the tags.

1.2.2 Go to `pynta`

```
cd pynta
```

1.2.3a Install `pynta`:

```
python setup.py install
```

1.2.3b (optional) If you do not have admin privileges (e.g. you use it on a supercomputer), do the following instead of 1.6a:

```
python setup.py install --user
```

**You should be ready to go now!**

Once finished using the workflow:

```
cd pynta
deactivate
```

## 2. How to run

### 2.1 Using Balsam

Before you run any `pynta` calculations, make sure your `balsam` DB is initialized and activated, e.g.

```bash
balsam init ~/myWorkflow
source balsamactivate ~/myWorkflow
```

You will need **4** files to run the workflow:

- `run_me.py` a python script that executes the workflow
- `run_me.sh` a bash script that submits jobs to the `balsam` database
- `inputR2S.py` a python script holding all user-modifiable parameters of the `pynta`
- `reactions.yaml` a yaml file with all reactions to be studied

An example `run_me.py` file:

```python
#!/usr/bin/env python3
from pynta.main import WorkFlow

# instantiate a WorkFlow() class
workflow = WorkFlow()
# create all input files
workflow.gen_job_files()
# execute the workflow
workflow.execute_all()
```

An example `run_me.sh` file:

```bash
#!/bin/bash
#SBATCH -J job_name        # name of the job e.g job_name = pynta_workflow
#SBATCH --partition=queue  # queue name e.g. queue = day-long-cpu
#SBATCH --nodes=x          # number of nodes e.g. x = 2
#SBATCH --ntasks=y         # number of CPUs e.g. 2 x 48 = y = 96
#SBATCH -e %x.err          # error file name
#SBATCH -o %x.out          # out file name

# load your quantum chemistry calculation package.
# Alternatively, provide a path to the preferred executable in 'inputR2S.py'
module load espresso
# activate balsam environment, e.g.
source balsamactivate ~/myWorkflow
# run python executable script
python3 $PWD/run_me.py
# required environment variable if using balsam branch serial-mode-perf and SLURM
export SLURM_HOSTS=$(scontrol show hostname)
# launch serial jobs
balsam launcher --job-mode=serial --wf-filter _ --limit-nodes=1 --num-transition-threads=1 &
# give some time to prevent time out before the sockets are ready
# for the quantum chemistry application, e.g. pw.x for Quantum Espresso
sleep 45
# launch mpi jobs
balsam launcher --job-mode=mpi --wf-filter QE_Sock --offset-nodes=x-1 --num-transition-threads=1 &
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

An example `inputR2S.py` file:

```python
from pathlib import Path
'''
####################################################
                    Basic Input
####################################################
'''
####################################################
# do you want to run surface optimization
optimize_slab = True
####################################################
# specify facet orientation, repeats of the slab+ads
# and repeats of the slab_opt unit cell
surface_types_and_repeats = {'fcc111': [(3, 3, 1), (1, 1, 4)],
                             'fcc100': [(3, 4, 1), (1, 1, 4)]}
####################################################
# surface atoms
symbol = 'Cu'
####################################################
# lattice constant
a = 3.6
####################################################
# vacuum in the z direction (Angstrem)
vacuum = 8.0
####################################################
# Quantum Espresso pseudopotantials and exe settings
# for DFT calculations
pseudo_dir = '/home/mgierad/espresso/pseudo'

pseudopotentials = "dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',"\
    + "H='H.pbe-kjpaw_psl.1.0.0.UPF',"\
    + "O='O.pbe-n-kjpaw_psl.1.0.0.UPF'," \
    + "C='C.pbe-n-kjpaw_psl.1.0.0.UPF')"

executable = '/home/mgierad/00_codes/build/q-e-qe-6.4.1/build/bin/pw.x'
####################################################
# Baslam settings
balsam_exe_settings = {'num_nodes': 1,  # nodes per each balsam job
                       'ranks_per_node': 48,  # cores per node
                       'threads_per_rank': 1
                       }
calc_keywords = {'kpts': (3, 3, 1),
                 'occupations': 'smearing',
                 'smearing': 'marzari-vanderbilt',
                 'degauss': 0.01,  # Rydberg
                 'ecutwfc': 40,  # Rydberg
                 'nosym': True,  # Allow symmetry breaking during optimization
                 'conv_thr': 1e-11,
                 'mixing_mode': 'local-TF'
                 }
####################################################
# Set up a working directory (this is default)
creation_dir = Path.cwd().as_posix()
####################################################
# filename of the .yaml file with reactions
yamlfile = 'reactions.yaml'
####################################################
# specify repeats of the surface in (x, y, z) direction
# repeats_surface = {'fcc111': (1, 1, 4), 'fcc100': (3, 4, 5)}
####################################################
# specify repeats of the surface in (x, y, z) direction
# repeats = {'fcc111': (3, 3, 1), 'fcc100': (3, 4, 1)}
####################################################
# specify the scaling factor to scale the bond distance
# between two atoms taking part in the reaction
scfactor = 1.4
####################################################
# specify the scaling factor to scale the target distance
# i.e. the average bond distance between adsorbate and
# the nearest surface metal atom
scfactor_surface = 1.0
####################################################
# species list
species_dict = {'rxn1': ['O', 'H'], 'rxn2': ['O', 'H']}
####################################################
# do you want to apply the scfactor_surface to the species 1?
scaled1 = False
####################################################
# do you want to apply scfactor_surface to the species 2?
scaled2 = False
####################################################

```

An example input files are also located at `./rmgcat_to_sella/example_run_files/`.

If you do not have a `.yaml` file with the reaction list but still want to use the work-flow, let me know. Also, stay tuned, as a version of `pynta` that can work without `.yaml` file is currently under development

### 2.2 Using SLURM only

**Warning `dev` branch uses SLURM scheduler to deal with the job dependencies. Be aware that it might be a bit buggy and do not fully support all the features implemented in the `master` branch.**

An example script (using `dev` branch - SLURM):

```python
#!/usr/bin/env python3
#SBATCH -J job_name        # name of the job e.g job_name = pynta_workflow
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
# import executable class of pynta
from pynta.main import WorkFlow
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

If you are using `pynta` or you wish to use it, let me know!

This work was supported by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, Chemical Sciences, Geosciences and Biosciences Division, as part of the Computational Chemistry Sciences Program.
