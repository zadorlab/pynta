# runprep.py  --  FireWorks lattice-constant optimization driver
#
# Mirrors the runpynta.py pattern: build the class, call one method, and let it
# build + add + launch the FireWorks workflow internally.
#   pyn = Pynta(...); pyn.generate_slab()
#   prep = Prep(...); prep.lattice_optimization()
#
# prep.lattice_optimization() builds N single-point energy Fireworks (one per
# scaled geometry) + a collector, adds the workflow to the launchpad, calls
# self.launch() (queue or local per how Prep was constructed), and -- with the
# default wait=True -- blocks until the collector writes lattice_constant.json,
# then returns the optimized lattice constant.
#
# NOTE on the Espresso API (ASE >= 3.23 / EspressoProfile):
# software_kwargs below intentionally stays a plain, JSON-serializable dict
# using the old-style "command" and "pseudo_dir" keys, because FireWorks
# serializes the spec into MongoDB (an EspressoProfile object cannot be
# stored there). The conversion to the new API,
#     EspressoProfile(command="srun ... pw.x", pseudo_dir=...),
# happens at calculator-construction time inside each rocket, via
# pynta.utils.adapt_espresso_kwargs / make_calculator. Do NOT put an
# EspressoProfile object in this dict, and do NOT include shell redirection
# ("< PREFIX.pwi > PREFIX.pwo") in the command -- the profile handles i/o
# itself (the adapter strips a redirect if one is present, for backwards
# compatibility with old submit scripts).

import json
import os
import re
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(".")))
from preprocessing import Prep


def mypause(flag):
    print("\033[31;49;1m {0:20s}\033[0m".format(flag))


def slurm_ntasks(default=1):
    """Resolve the MPI rank count from the SLURM environment, robustly."""
    for var in ("SLURM_NTASKS", "SLURM_NPROCS"):
        v = os.environ.get(var, "")
        if v.isdigit() and int(v) > 0:
            return int(v)
    tpn = os.environ.get("SLURM_NTASKS_PER_NODE", "")
    nodes = os.environ.get("SLURM_JOB_NUM_NODES") or os.environ.get("SLURM_NNODES") or ""
    m = re.match(r"(\d+)", tpn)
    if m and nodes.isdigit():
        return int(m.group(1)) * int(nodes)
    c = os.environ.get("SLURM_CPUS_ON_NODE", "")
    if c.isdigit() and int(c) > 0:
        return int(c)
    return default


# --- FireWorks config + inputs (point these at your real files) -------------
LAUNCHPAD_PATH = "/home/shikim/local_launchpad.yaml"
FWORKER_PATH   = "/home/shikim/my_fworker.yaml"
QADAPTER_PATH  = "/home/shikim/qe_qadapter.yaml"
CUSTOM_SLAB    = "/home/shikim/pynta/pynta/preprocessing/custom_surfaces/single_atom_PtAu.xyz"

# Set queue=True to fan the scan points out as separate SLURM jobs via the
# qadapter (true parallel). Set queue=False to run them locally in this process
# (sequential rapidfire). Match runpynta.py, which uses queue=True.
USE_QUEUE = True


def run():

    # The launcher each energy Firework uses to run pw.x. With queue=True the
    # qadapter provides the allocation, so resolving ntasks here is only a
    # sanity print; the rocket inherits SLURM_NTASKS from its own job.
    # Bare launch command, no "< PREFIX.pwi > PREFIX.pwo" (see API note above).
    ntasks = slurm_ntasks()
    command = f"srun -n {ntasks} --mpi=pmi2 --cpu-bind=cores pw.x"
    print(f"QE launch command (per Firework): {command}  (ntasks here={ntasks})")

    # Old-style flat dict on purpose: JSON-serializable for the FireWorks spec.
    # adapt_espresso_kwargs converts it to
    #   profile=EspressoProfile(command=..., pseudo_dir=...) + input_data={...}
    # inside each rocket on ASE >= 3.23 (and leaves it as-is on older ASE).
    espresso_kwargs = {
        "kpts": (5, 5, 1),
        "ecutwfc": 70,
        "degauss": 0.02,
        "input_dft": "BEEF-VDW",
        "mixing_mode": "plain",
        "occupations": "smearing",
        "smearing": "marzari-vanderbilt",
        "pseudo_dir": "/home/shikim/espresso/pseudo",
        "pseudopotentials": {
            "Pt": "Pt.pbe-n-kjpaw_psl.1.0.0.UPF",
            "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
            "Au": "Au.pbe-n-kjpaw_psl.1.0.0.UPF",
        },
        "command": command,
    }

    prep = Prep(
        metal="Pt",
        surface_type="fcc111",
        a0=3.96,
        software="Espresso",
        vacuum=10,
        slab=CUSTOM_SLAB,
        path=os.getcwd(),
        label="PtAu_lattice",
        pbc=(True, True, False),
        launchpad_path=LAUNCHPAD_PATH,
        fworker_path=FWORKER_PATH,
        queue=USE_QUEUE,
        queue_adapter_path=QADAPTER_PATH if USE_QUEUE else None,
        njobs_queue=20,
        software_kwargs=espresso_kwargs,
        lattice_opt_software_kwargs=espresso_kwargs,
    )

    # Build + add + launch the workflow, then block for the result.
    a_opt = prep.lattice_optimization(da=0.1, scan_step=0.02, centered=False)
    print(f"Optimized lattice constant: a = {a_opt:.6f} Angstrom")
    print("Lattice optimization finished")


if __name__ == "__main__":
    mypause("Before Run")
    run()
    mypause("After  Run")