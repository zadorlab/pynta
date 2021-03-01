![Build status](https://github.com/zadorlab/pynta/workflows/Build%20status/badge.svg)
![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/zadorlab/pynta?label=version)
![GitHub last commit](https://img.shields.io/github/last-commit/zadorlab/pynta?label=last%20modified)
![GitHub commits since tagged version](https://img.shields.io/github/commits-since/zadorlab/pynta/v1.1.0?color=9cf)
![GitHub](https://img.shields.io/badge/License-GPLv3-orange)
![GitHub top language](https://img.shields.io/github/languages/top/zadorlab/pynta?color=brightgreen)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/zadorlab/pynta.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/zadorlab/pynta/context:python)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/zadorlab/pynta.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/zadorlab/pynta/alerts/)
![GitHub repo size](https://img.shields.io/github/repo-size/zadorlab/pynta)
![GitHub all releases](https://img.shields.io/github/downloads/zadorlab/pynta/total?color=yellow)

# Pynta

<!-- ![workflow_idea](./workflow_idea.png) -->

`pynta` is an automated workflow for reaction path exploration on metallic surfaces.

The `pynta` code is designed to automatically characterize chemical reactions
relevant to heterogeneous catalysis. In particular, it spawns and processes a
large number of ab initio quantum chemistry calculations to study gas-phase
reactions on crystal facets. It is designed to run on petascale and upcoming
exascale machines.

The code systematically places adsorbates on crystal facets, by enumerating the
various unique crystal sites and considering the symmetry of the adsorbates and
the surface. The structures are systematically perturbed, to try to investigate
all possible adsorption geometries. Following this, the code performs
modifications to these structures and searches for saddle points on the
potential energy landscape in order to characterize the reactions these
adsorbates can undergo. After a series of such calculations the code arrives at
well-characterized reaction pathways, which can be in turn used to calculate
rate coefficients and can be implemented in microkinetic mechanisms/models.

`pynta` is designed to work with the workflow code
[`balsam`](https://balsam.readthedocs.io/en/latest/), which
enables a seamless running of embarrassingly parallel computations on
supercomputers. `pynta` includes several so-called apps that are run through
`balsam`, and appear as one monolithic job submission in the queue system.
Another level of parallelism is gained from the ab initio programs, which may
or may not include GPU acceleration. We use the Atomic Simulation Environment
(ASE) that enables coupling to a large variety of ab initio quantum
chemistry codes.

<center><img src='./workflow_idea.png' style="width:400px"></center>
<br>

# 1. Installation

The following instruction assumes that you have several softwares installed on your system, such as:

- A queue system e.g. `Slurm`, `PBS` or `Cobalt`
- An MPI e.g. `OpenMPI`
- A math library e.g. `OpenBlas/LAPACK` or `Intel's MKL`
- A compilier suite e.g. `GCC` or `Intel`

## 1.1 Install all prerequisites

### 1.1.1 Set the location you want to use for `pynta` build. It is suggested to create a virtual environment e.g 'pynta':

```bash
cd <path/whr/to/build/pynta>
python3 -m venv pynta
```

### 1.1.2 Activate your virtual environment:

```bash
source pynta/bin/activate
```

### 1.1.3 (_Optional_) Install required python packages (if you skip this process, `pynta` installer will do this later.)

```bash
pip3 install matplotlib<3.2 spglib==1.14.1.post0 networkx<2.4 ase==3.19.0 scipy==1.3.1 numpy==1.18.1 PyYAML==5.3.1 sella==1.0.3

```

### 1.1.4 Download [PostgreSQL](https://www.enterprisedb.com/download-postgresql-binaries) precompiled binaries that suits your system and add `path_to_PostgreSQL/pgsql/bin` to your `PATH` by modifying `~/.bashrc`

```bash
echo 'export PATH=path_to_PostgreSQL/pgsql/bin:$PATH' >> ~/.bashrc
```

1.1.5 Install [mpi4py](https://github.com/mpi4py/mpi4py.git):

```bash
git clone https://github.com/mpi4py/mpi4py.git
cd mpi4py
python3 setup.py install
cd ../
```

Make sure it works by running

```python
>>> mpirun -n 2 python3 -c 'from mpi4py import MPI; print(MPI.COMM_WORLD.Get_rank())'
0
1
```

1.1.6 Install [balsam](https://github.com/balsam-alcf/balsam.git) using [serial-mode-perf](https://github.com/balsam-alcf/balsam/tree/serial-mode-perf) branch.

```bash
git clone https://github.com/balsam-alcf/balsam.git -b serial-mode-perf
cd balsam
python3 setup.py install
cd ../
```

Make sure it works by running tests posted on the [balsam](https://github.com/balsam-alcf/balsam.git) GitHub page.

### 1.1.7 Install [xTB-python](https://github.com/grimme-lab/xtb-python) following instruction provided there. Make sure to correctly link all required libraries. For example:

- using `OpenBlas` and `GNU` based compilers:

```bash
git clone https://github.com/grimme-lab/xtb-python.git
cd xtb-python
git submodule update --init
LDFLAGS="-L/opt/custom/OpenBLAS/0.3.7/lib" meson setup build --prefix=$PWD --libdir=xtb/xtb --buildtype release --optimization 2 -Dla_backend=openblas
ninja -C build install
pip install -e .
```

- using `MKL` and Intel Compilers:

```bash
git clone https://github.com/grimme-lab/xtb-python.git
cd xtb-python
git submodule update --init
# (ALCF's Theta specific)
# conda instal cffi
# module swap PrgEnv-intel PrgEnv-cray; module swap PrgEnv-cray PrgEnv-intel
CC=icc CXX=icpc FC=ifort meson setup build --prefix=$PWD --libdir=xtb -Dla_backed=mkl -Dpy=3 --buildtype release --optimization 2
ninja -C build install
python3 setup.py install
ln -s ./xtb/xtb/_libxtb*.so ./xtb
```

Make sure it works by running:

```python
>>> from xtb.ase.calculator import XTB
>>> from ase.build import molecule
>>> import math
>>> atoms = molecule('H2O')
>>> atoms.calc = XTB(method="GFN2-xTB")
>>> total_ener = atoms.get_potential_energy()
>>> is_close = math.isclose(total_ener,-137.9677758730299)
>>> print (total_ener)  # True -> that's good; False -> something is wrong
```

**Warning - You might be getting SEGFAULT error -**

`Segmentation Fault (Core dumped)`

**while executing any** `xTB-python` **job, especially for a relatively large molecules. The easiest solution is to unlimit the system stack to avoid stack overflows. In** `bash` **try:**

```
ulimit -s unlimited
```

If `xTB-python` still fails, try to install [xtb](https://github.com/grimme-lab/xtb) and test `xTB` itself for any errors.

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

Then, rebuild `xTB-python` on your system ignoring `git submodule update --init` and linking you current `xTB` installation.

## 1.2 Install `pynta`

### 1.2.1 Clone the project in your preferable location.

```
git clone https://gitlab-ex.sandia.gov/mgierad/pynta.git
```

Usually, `master` branch should be fine. If somehow it is not working, make sure to switch to the latest stable release version by checking the tags.

### 1.2.2 Go to `pynta` directory

```
cd pynta
```

### 1.2.3a Install `pynta`:

```
python setup.py install
```

### 1.2.3b (optional) If you do not have admin privileges (e.g. you will be using `pynta` on a supercomputer), do the following instead of 1.2.3a:

```
python setup.py install --user
```

**You should be good to go with `pynta`**

Once finished using the workflow, type:

```
cd pynta
deactivate
```

# 2. How to run

Before you run any `pynta` calculations, make sure your `balsam` DB is initialized and activated, e.g.

```bash
balsam init ~/myWorkflow
source balsamactivate ~/myWorkflow
```

Your prompt should change to something like:

```bash
~[BalsamDB: myWorkflow] <username>@<host>:
```

You will need **4** files to run the workflow:

- `run_me.py` a python script that executes the workflow or alternatively,
  `restart_me.py` to restart unfinished calculations
- `run_me.sh` a bash script that submits jobs to the `balsam` database
- `input.json` a json file holding all user-modifiable parameters of the `pynta`
- `reactions.yaml` a yaml file with all reactions to be studied

An example `run_me.py` file:

```python
#!/usr/bin/env python3
from pynta.main import WorkFlow


def run():
    # instantiate a WorkFlow() class
    workflow = WorkFlow()
    # create all input files
    workflow.gen_job_files()
    # execute the workflow
    workflow.execute_all()


if __name__ == '__main__':
    run()
```

An example `run_me.py` file:

```python
#!/usr/bin/env python3
from pynta.main import WorkFlow


def restart():
    return WorkFlow().restart()


if __name__ == '__main__':
    restart()
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

# load your quantum chemistry calculation package or provide a path to the
# executable in 'inputR2S.py'
module load espresso

# activate balsam environment, e.g.
source balsamactivate ~/myWorkflow

# run python executable script
python3 $PWD/run_me.py

# required environment variable if using balsam branch serial-mode-perf and SLURM
export SLURM_HOSTS=$(scontrol show hostname)

# launch serial jobs
balsam launcher --job-mode=serial --wf-filter _ --limit-nodes=1 --num-transition-threads=1 &

# wait a bit to prevent timeing out you jobs before the sockets are ready
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

```json
{
  "quantum_chemistry": "espresso",
  "optimize_slab": true,
  "surface_types_and_repeats": {
    "fcc111": [
      [3, 3, 1],
      [1, 1, 4]
    ]
  },
  "metal_atom": "Cu",
  "a": 3.6,
  "vacuum": 8.0,
  "pseudo_dir": "/home/mgierad/espresso/pseudo",
  "pseudopotentials": "dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF',H='H.pbe-kjpaw_psl.1.0.0.UPF',O='O.pbe-n-kjpaw_psl.1.0.0.UPF',C='C.pbe-n-kjpaw_psl.1.0.0.UPF',N='N.pbe-n-kjpaw_psl.1.0.0.UPF')",
  "executable": "/home/mgierad/00_codes/build/q-e-qe-6.4.1/build/bin/pw.x",
  "node_packing_count": 48,
  "balsam_exe_settings": {
    "num_nodes": 1,
    "ranks_per_node": 48,
    "threads_per_rank": 1
  },
  "calc_keywords": {
    "kpts": [3, 3, 1],
    "occupations": "smearing",
    "smearing": "marzari-vanderbilt",
    "degauss": 0.01,
    "ecutwfc": 40,
    "nosym": true,
    "conv_thr": 1e-11,
    "mixing_mode": "local-TF"
  },
  "yamlfile": "reactions.yaml",
  "scfactor": 1.4,
  "scfactor_surface": 1.0,
  "scaled1": false,
  "scaled2": false
}
```

For meaning of theses keywords, please check [documentation]()

An example input files are also located at `./pynta/example_run_files/`

# 3. Documentation

A full documentation is available [here](https://zadorlab.github.io/pynta/).

# Acknowledgments

If you are using `pynta` or you wish to use it, let me know!

# Copyright Notice

For five (5) years from 1/19/2021 the United States Government is granted for
itself and others acting on its behalf a paid-up, nonexclusive, irrevocable
worldwide license in this data to reproduce, prepare derivative works, and
perform publicly and display publicly, by or on behalf of the Government. There
is provision for the possible extension of the term of this license. Subsequent
to that period or any extension granted, the United States Government is
granted for itself and others acting on its behalf a paid-up, nonexclusive,
irrevocable worldwide license in this data to reproduce, prepare derivative
works, distribute copies to the public, perform publicly and display publicly,
and to permit others to do so. The specific term of the license can be
identified by inquiry made to National Technology and Engineering Solutions of
Sandia, LLC or DOE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF
ENERGY, NOR NATIONAL TECHNOLOGY AND ENGINEERING SOLUTIONS OF SANDIA, LLC, NOR
ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LEGAL RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS
USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

Any licensee of this software has the obligation and responsibility to abide by
the applicable export control laws, regulations, and general prohibitions
relating to the export of technical data. Failure to obtain an export control
license or other authority from the Government may result in criminal liability
under U.S. laws.

# License Notice

Copyright 2021 National Technology & Engineering Solutions of Sandia,
LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the
U.S. Government retains certain rights in this software.
