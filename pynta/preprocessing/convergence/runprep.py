# preprocessing_input.py
import os
import re
from preprocessing import Prep
# from ase.io import read   # only needed if you set prep.slab = read(...) manually


def mypause(flag):
    print("\033[31;49;1m {0:20s}\033[0m".format(flag))


def slurm_ntasks(default=1):
    """Resolve the MPI rank count from the SLURM environment, robustly.

    Falls back through SLURM_NTASKS -> SLURM_NPROCS -> (per-node x nodes) ->
    CPUs-on-node. Returns `default` only if nothing is set (i.e. not in a job).
    """
    for var in ("SLURM_NTASKS", "SLURM_NPROCS"):
        v = os.environ.get(var, "")
        if v.isdigit() and int(v) > 0:
            return int(v)
    tpn = os.environ.get("SLURM_NTASKS_PER_NODE", "")          # may be "36" or "36(x2)"
    nodes = os.environ.get("SLURM_JOB_NUM_NODES") or os.environ.get("SLURM_NNODES") or ""
    m = re.match(r"(\d+)", tpn)
    if m and nodes.isdigit():
        return int(m.group(1)) * int(nodes)
    c = os.environ.get("SLURM_CPUS_ON_NODE", "")
    if c.isdigit() and int(c) > 0:
        return int(c)
    return default


def run():

    # QE run command for ASE >=3.23: launcher + binary ONLY.
    # Do NOT add `-input`/`-in` or `< PREFIX.pwi > PREFIX.pwo`; ASE appends those.
    ntasks = slurm_ntasks()
    if ntasks <= 1:
        print(f"WARNING: resolved ntasks={ntasks}; pw.x will run SERIAL. "
              f"Are you inside the 16-task allocation (Job B), not a 2-core shell?")

    # Launch with srun: it inherits SLURM_NTASKS from the *current allocation* and
    # pins each rank to an allocated core. No `-n` needed -> uses the allocation's
    # task count (16 in Job B). The cores come from the allocation, not this line.
    #
    # IMPORTANT: srun must hand pw.x the right PMI or you get N independent serial
    # copies again ("running on 1 processors"). Pick the plugin to match pw.x's MPI:
    #   openmpi  -> --mpi=pmix      mpich/intel -> --mpi=pmi2
    # Check:  ldd $(which pw.x) | grep -i -E 'mpi|pmi'   and   srun --mpi=list
    # pw.x links conda MPICH (libmpi.so.12), and this cluster's srun offers pmi2.
    # MPICH speaks PMI-2, so --mpi=pmi2 forms ONE N-rank job; plain srun / pmix do
    # not bootstrap it -> N serial singletons ("running on 1 processors").
    # Must run inside the 16-core allocation (Job B); inherits SLURM_NTASKS (=16).
    command = "srun --mpi=pmi2 --cpu-bind=cores pw.x"
    # Equivalent (single node): the matching conda mpirun ->
    #   command = f"mpirun -np {ntasks} pw.x"
    print(f"QE launch command: {command}  (allocation tasks: {ntasks})")

    CUSTOM_SLAB = "/home/shikim/pynta/pynta/preprocessing/custom_surfaces/single_atom_PtAu.xyz"

    prep = Prep(
        metal="Pt",
        surface_type="fcc111",
        a0=3.96,
        software="Espresso",
        fmax=0.05,
        vacuum=10,
        frozen_layers=3,          # used by optimize_slab() to freeze bottom layers
        slab=CUSTOM_SLAB,         # path OR an ase.Atoms; read in automatically
        software_kwargs={         # production / convergence settings
            "kpts": (3, 3, 1),
            "ecutwfc": 40,
            "ecutrho": 160,
            "occupations": "smearing",
            "pseudo_dir": "/home/shikim/espresso/pseudo",
            "smearing": "marzari-vanderbilt",
            "degauss": 0.01,
            "input_dft": "BEEF-VDW",
            "nosym": True,
            "conv_thr": 1e-6,
            "pseudopotentials": {
                "Pt": "Pt.pbe-n-kjpaw_psl.1.0.0.UPF",
                "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
                "Au": "Au.pbe-n-kjpaw_psl.1.0.0.UPF",
            },
            "command": command,
        },
        lattice_opt_software_kwargs={   # used by calculate_lattice_parameters()
            "kpts": (5, 5, 1),          # slab: sample in-plane only, not vacuum (z)
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
        },
    )

    # =====================================================================
    # MODE A -- EXPLICIT: pick exactly the function you want on this slab.
    # =====================================================================

    # Lattice constant of the custom surface (in-plane scan, vacuum z fixed):
    a_opt = prep.calculate_lattice_parameters(da=0.1, scan_step=0.02, inplane_only=None,
                                     centered=False)
    print(f"Optimized lattice constant: a = {a_opt:.6f} Angstrom")

    # Other explicit functions available on the same slab (uncomment as needed):
    # prep.optimize_slab(fmax=0.05)                       # relax atomic positions
    # prep.run_ecut_convergence(range(30, 71, 10))        # ecut scan @ kpts from software_kwargs
    # prep.run_kpoint_convergence([(2,2,1),(3,3,1),(4,4,1)])   # k-point scan @ fixed ecut
    # prep.run_convergence(range(30,71,10), [(2,2,1),(3,3,1)]) # full 2D grid

    # =====================================================================
    # MODE B -- KEYWORD: pass keywords and let run_selected() dispatch.
    # (Don't combine with Mode A in one run.)
    #   prep = Prep(..., optimize_lattice=True)                  -> lattice opt
    #   prep = Prep(..., ecut_values=range(30,71,10))            -> ecut scan
    #   prep = Prep(..., kmesh_values=[(2,2,1),(3,3,1)])         -> k-point scan
    #   prep = Prep(..., ecut_values=..., kmesh_values=...)      -> both 1D scans
    #   prep = Prep(..., ecut_values=..., kmesh_values=...,
    #               convergence_grid=True)                       -> add the 2D grid
    #   prep.run_selected()
    # =====================================================================

    print("Preprocessing finished")


if __name__ == "__main__":

    mypause("Before Run")
    run()
    mypause("After  Run")