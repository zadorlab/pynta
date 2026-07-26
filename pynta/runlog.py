"""Curated pynta run log -- a human-readable timeline of what a run is actually doing, kept
SEPARATE from FireWorks' queue/rocket chatter ("job submission successful", "N jobs in the queue").

Log meaningful run events through this module's ``log`` (the dedicated "pynta" logger): stage
starts/finishes, per-calculation start/finish, and warnings you'd want to notice (an optimization
killed on its time limit, a central with no ZPE, ...). Because it is its OWN logger, FireWorks'
loggers never leak into it.

Opt-in file: set the env var ``PYNTA_LOG_FILE=/path/to/run.log`` (e.g. in your qadapter's
``pre_rocket`` or your shell rc) and every firework appends its pynta events to that one file, so the
timeline survives across fireworks instead of being scattered through per-job SLURM stdout. If the
var is unset, nothing changes (events still propagate to the console/FW stdout as before).
"""
import os
import logging

log = logging.getLogger("pynta")


def configure_from_env():
    """Attach a FileHandler for $PYNTA_LOG_FILE to the pynta logger (idempotent per process)."""
    logpath = os.environ.get("PYNTA_LOG_FILE")
    if not logpath:
        return log
    logpath = os.path.abspath(logpath)
    for h in log.handlers:  # don't double-attach if this process re-imports/re-configures
        if isinstance(h, logging.FileHandler) and getattr(h, "baseFilename", None) == logpath:
            return log
    fh = logging.FileHandler(logpath)  # append mode; multiple fireworks share the file
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log.addHandler(fh)
    if log.level == logging.NOTSET:
        log.setLevel(logging.INFO)
    # propagate=True (default): events also reach the console/FW stdout, so nothing is lost
    return log


configure_from_env()
