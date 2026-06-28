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
    # Resolve ranks from SLURM (not a literal $VAR: ASE runs pw.x without a shell,
    # so a literal $SLURM_NTASKS would never expand).
    ntasks = slurm_ntasks()
    if ntasks <= 1:
        print(f"WARNING: resolved ntasks={ntasks}; pw.x will run SERIAL. "
              f"Check that the job requested >1 task and that python is NOT "
              f"launched via `srun` (nested srun is limited to 1 task).")

    # IMPORTANT: launch with the mpirun that BELONGS TO THE SAME ENV as pw.x.
    # `srun pw.x` with a conda-built QE spawns N *independent serial* pw.x copies
    # (each reports "running on 1 processors") instead of one N-rank job. Using the
    # env's own mpirun forms a single MPI communicator -> one parallel pw.x.
    command = f"mpirun -np {ntasks} pw.x"
    # Single node assumed. Alternatives if mpirun isn't right for your build:
    #   - srun with explicit PMI:   command = f"srun --mpi=pmi2 -n {ntasks} pw.x"
    #                          or:  command = f"srun --mpi=pmix -n {ntasks} pw.x"
    #     (pick what `srun --mpi=list` shows AND what `ldd $(which pw.x)` links).
    print(f"QE launch command: {command}")

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
    a_opt = prep.calculate_lattice_parameters(da=0.1, scan_step=0.01)
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