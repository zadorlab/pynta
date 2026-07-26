"""Curated pynta run log -- a human-readable timeline of what a run is actually doing, kept
SEPARATE from FireWorks' queue/rocket chatter ("job submission successful", "N jobs in the queue").

Log meaningful run events through this module's ``log`` (the dedicated "pynta" logger): stage
starts/finishes, per-calculation start/finish, and warnings you'd want to notice (an optimization
killed on its time limit, a central with no ZPE, ...). Because it is its OWN logger, FireWorks'
loggers never leak into it.

Opt-in: set ``PYNTA_LOG_FILE=/path/to/run.log`` (e.g. in your qadapter's ``pre_rocket`` or shell rc).

PER-PROCESS files (important on a cluster): every firework runs in its own process, often on a
different node, so if they all appended to ONE shared file the writes would interleave -- and on a
network filesystem (NFS/Lustre) append is not atomic, producing null-byte-corrupted output. Instead
each process writes its OWN file under ``<PYNTA_LOG_FILE>.d/<host>-<pid>.log`` (sole writer -> safe on
any filesystem). Each file is already in chronological order, so merge the whole run's timeline on
demand with:

    sort -k1,2 /path/to/run.log.d/*.log        # or: -m, the per-file streams are pre-sorted

The ``<base>.d`` directory keeps the per-process files tidily out of the way. If the var is unset,
nothing is written (events still propagate to the console/FW stdout as before).
"""
import os
import socket
import logging

log = logging.getLogger("pynta")


def configure_from_env():
    """Attach a per-process FileHandler under $PYNTA_LOG_FILE.d/ (idempotent within a process)."""
    base = os.environ.get("PYNTA_LOG_FILE")
    if not base:
        return log
    base = os.path.abspath(base)
    logdir = base + ".d"
    try:
        os.makedirs(logdir, exist_ok=True)  # concurrent creation is fine with exist_ok
    except OSError:
        pass
    # sole-writer file for THIS process: no cross-process/cross-node contention -> no corruption
    logpath = os.path.join(logdir, "%s-%d.log" % (socket.gethostname(), os.getpid()))
    for h in log.handlers:  # don't double-attach if this process re-imports/re-configures
        if isinstance(h, logging.FileHandler) and getattr(h, "baseFilename", None) == logpath:
            return log
    fh = logging.FileHandler(logpath)  # append mode, but this pid is the only writer
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log.addHandler(fh)
    if log.level == logging.NOTSET:
        log.setLevel(logging.INFO)
    # propagate=True (default): events also reach the console/FW stdout, so nothing is lost
    return log


configure_from_env()
