.. role:: raw-html-m2r(raw)
   :format: html


Pynta
=====

``pynta`` is an automated workflow for reaction path exploration on metallic surfaces.

The ``pynta`` code is designed to automatically characterize chemical reactions
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

Pynta is designed to work with the workflow code
`balsam <https://balsam.readthedocs.io/en/latest/>`_ , which enables a seamless
running of embarrassingly parallel computations on supercomputers. ``pynta``
includes several so-called apps that are run through ``balsam``, and appear as
one monolithic job submission in the queue system. Another level of parallelism
is gained from the ab initio programs, which may or may not include GPU
acceleration. We use the Atomic Simulation Environment (ASE) that enables
coupling to a large variety of ab initio quantum chemistry codes.


.. image:: ../workflow_idea.png
   :width: 800
   :align: center

Installation
============


The following instruction assumes that you have several softwares installed on your system, such as:

* A queue system e.g. ``Slurm``\ , ``PBS`` or ``Cobalt``
* An MPI e.g. ``OpenMPI``
* A math library e.g. ``OpenBlas/LAPACK`` or ``Intel's MKL``
* A compilier suite e.g. ``GCC`` or ``Intel``

Install all prerequisites
-------------------------

Set the location you want to use for ``pynta`` build. It is suggested to create a virtual environment e.g 'pynta':
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   $ cd <path/whr/to/build/pynta>
   $ python3 -m venv pynta

Activate your virtual environment:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   $ source pynta/bin/activate

(\ *Optional*\ ) Install required python packages (if you skip this process, ``pynta`` installer will do this later.)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   $ pip3 install matplotlib<3.2 spglib==1.14.1.post0 networkx<2.4 ase==3.19.0 scipy==1.3.1 numpy==1.18.1 PyYAML==5.3.1 sella==1.0.3

Download `PostgreSQL <https://www.enterprisedb.com/download-postgresql-binaries>`_ precompiled binaries that suits your system and add ``path_to_PostgreSQL/pgsql/bin`` to your ``PATH`` by modifying ``~/.bashrc``:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   $ echo "export PATH=path_to_PostgreSQL/pgsql/bin:$PATH' >> ~/.bashrc"


Install `mpi4py <https://github.com/mpi4py/mpi4py.git>`_\ :
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   $ git clone https://github.com/mpi4py/mpi4py.git
   $ cd mpi4py
   $ python3 setup.py install
   $ cd ../

Make sure it works by running

.. code-block:: bash

   $ srun -n 2 python3 -c 'from mpi4py import MPI; print(MPI.COMM_WORLD.Get_rank())'
   0
   1

Install `balsam <https://github.com/balsam-alcf/balsam.git>`_ using `serial-mode-perf <https://github.com/balsam-alcf/balsam/tree/serial-mode-perf>`_ branch.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   git clone https://github.com/balsam-alcf/balsam.git -b serial-mode-perf
   cd balsam
   python3 setup.py install
   cd ../

Make sure it works by running tests posted on the `balsam <https://github.com/balsam-alcf/balsam.git>`_ GitHub page.

Install `xTB-python <https://github.com/grimme-lab/xtb-python>`_ following instruction provided there. Make sure to correctly link all required libraries. For example:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


* using ``OpenBlas`` and ``GNU`` based compilers:

.. code-block:: bash

   git clone https://github.com/grimme-lab/xtb-python.git
   cd xtb-python
   git submodule update --init
   LDFLAGS="-L/opt/custom/OpenBLAS/0.3.7/lib" meson setup build --prefix=$PWD --libdir=xtb/xtb --buildtype release --optimization 2 -Dla_backend=openblas
   ninja -C build install
   pip install -e .


* using ``MKL`` and Intel Compilers:

.. code-block:: bash

   git clone https://github.com/grimme-lab/xtb-python.git
   cd xtb-python
   git submodule update --init
   # (Theta specific)
   # conda instal cffi
   # module swap PrgEnv-intel PrgEnv-cray; module swap PrgEnv-cray PrgEnv-intel
   CC=icc CXX=icpc FC=ifort meson setup build --prefix=$PWD --libdir=xtb -Dla_backed=mkl -Dpy=3 --buildtype release --optimization 2
   ninja -C build install
   pip install -e .

Make sure it works by running:

.. code-block:: python

   >>> from ase.build import molecule
   >>> from xtb.ase.calculator import XTB
   >>> atoms = molecule('H2O')
   >>> atoms.calc = XTB(method="GFN2-xTB")
   >>> total_ener = atoms.get_potential_energy()
   >>> total_ener
   -137.9677758730299

.. warning:: You might be getting SEGFAULT error - ``Segmentation Fault (Core dumped)`` while executing any ``xTB-python`` job, especially for a relatively large molecules. The easiest solution is to unlimit the system stack to avoid stack overflows. In ``bash`` try

.. code-block::

   ulimit -s unlimited

If ``xTB-python`` still fails, try to install `xtb <https://github.com/grimme-lab/xtb>`_ and test ``xTB`` itself for any errors.

.. code-block:: bash

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

Then, rebuild ``xTB-python`` on your system ignoring ``git submodule update --init`` and linking you current ``xTB`` installation.

Install ``pynta``
-----------------

Clone the project in your preferable location.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

   git clone https://gitlab-ex.sandia.gov/mgierad/pynta.git

Usually, ``master`` branch should be fine. If somehow it is not working, make sure to switch to the latest stable version by checking the tags.

Go to ``pynta`` directory
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

   cd pynta

Install ``pynta``\ :
^^^^^^^^^^^^^^^^^^^^

.. code-block::

   python setup.py install

(\ *Optional*\ )  If you do not have admin privileges (e.g. you use it on a supercomputer), do the following instead of 1.6a:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

   python setup.py install --user

.. note::

   You should be ready to use ``pynta`` at this point!

Once finished using the workflow:

.. code-block::

   cd pynta
   deactivate

How to run
==========

Using Balsam
------------

Before you run any ``pynta`` calculations, make sure your ``balsam`` DB is initialized and activated, e.g.

.. code-block:: bash

   $ balsam init ~/myWorkflow
   $ source balsamactivate ~/myWorkflow

You will need **4** files to run the workflow:


* ``run_me.py`` a python script that executes the workflow or alternatively, ``restart_me.py`` to restart unfinished calculations
* ``run_me.sh`` a bash script that submits jobs to the ``balsam`` database
* ``input.json`` a json file holding all user-modifiable parameters of the ``pynta``
* ``reactions.yaml`` a yaml file with all reactions to be studied

An example ``run_me.py`` file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

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


An example ``restart_me.py`` file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   ##!/usr/bin/env python3
   from pynta.main import WorkFlow


   def restart():
      return WorkFlow().restart()


   if __name__ == '__main__':
      restart()


An example ``run_me.sh`` file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

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

An example ``reactions.yaml`` file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: yaml

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

An example ``inputR2S.py`` file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

   {
    "quantum_chemistry": "espresso",
    "optimize_slab": true,
    "surface_types_and_repeats": {
        "fcc111": [
            [
                3,
                3,
                1
            ],
            [
                1,
                1,
                4
            ]
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
        "kpts": [
            3,
            3,
            1
        ],
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
    "scaled2": false}

An example input files are also located at ``./pynta/example_run_files/``.

Keyword parameters
==================

.. list-table::
   :widths: 25, 5, 25, 25
   :header-rows: 1

   * - Parameter
     - Description
     - Type
     - Avaiable options/examples
   * - `quantum_chemistry`
     - specify which quantum chemistry calculator to use
     - str
     - ``'espresso'`` , ``'nwchem'`` 
   * - `yamlfile`
     - a name of the  yml file with reactions
     - str
     - ``reactions.yaml``
   * - `optimize_slab`
     - do you want to run surface optimization
     - bool
     - ``True``, ``False``
   * - `surface_types_and_repeats`
     - specify facet orientation, repeats of the slab+ads and repeats of the slab_opt unit cell
     - Dict[str:List[Tuple[int]]]
     - ``{'fcc111': [(3, 3, 1), (1, 1, 4)]}``
   * - `metal_atom`
     - surface atoms
     - str
     - ``'Cu'``
   * - `a`
     - lattice constant
     - float
     - ``3.6``
   * - `vacuum`
     - vacuum in the z direction (Angstrem)
     - float
     - ``8.0``
   * - `pseudo_dir`
     - for Quantum Espresso only - path to pseudopotetions directory
     - str
     - ``'/home/user/espresso/pseudo'``
   * - `pseudopotentials`
     - for Quantum Espresso only - which pseudopotentials to use (str preferable)
     - str -> Dict[str, str]
     - ``"dict(Cu='Cu.pbe-spn-kjpaw_psl.1.0.0.UPF', H='H.pbe-kjpaw_psl.1.0.0.UPF')"``
   * - `executable`
     - a path to executable of quantum chemistry package
     - str
     - ``'/home/user/build/qe/bin/pw.x'``
   * - `node_packing_count`
     - number of cores per node - balsam specific keyword
     - int
     - ``48``
   * - `balsam_exe_settings`
     - how to use avaiable computational resources
     - Dict[str, any]
     - ``balsam_exe_settings = {'num_nodes': 1, 'ranks_per_node': node_packing_count, 'threads_per_rank': 1``
   * - `calc_keywords`
     - quantum chemistry packages options - see manual of your chosen package to put a correct keywords here
     - Dict[str, any]
     - ``{'kpts': (3, 3, 1), 'occupations': 'smearing', 'smearing': 'marzari-vanderbilt', 'degauss': 0.01, 'ecutwfc': 40, 'nosym': True, 'conv_thr': 1e-11, 'mixing_mode': 'local-TF'}``
   * - `scfactor`
     - specify a scaling factor to scale the bond distance between two atoms taking part in the reaction
     - float
     - ``1.4``
   * - `scfactor_surface`
     - specify the scaling factor to scale the target distance i.e. the average bond distance between adsorbate and the nearest surface metal atom
     - float
     - ``1.0``
   * - `scaled1`
     - apply the `scfactor_surface` to the species 1?
     - bool
     - ``True`` or ``False``
   * - `scaled2`
     - apply the `scfactor_surface` to the species 2?
     - bool
     - ``True`` or ``False``