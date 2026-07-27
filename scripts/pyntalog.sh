#!/usr/bin/env bash
# pyntalog -- merge the per-process pynta run logs into one chronological stream.
#
# When PYNTA_LOG_FILE is set, each firework writes its own file
#   <PYNTA_LOG_FILE>.d/<host>-<pid>.log
# (per-process, so concurrent writers on a network FS can't corrupt one shared file). Each file is
# already time-ordered and every line starts with "YYYY-MM-DD HH:MM:SS,mmm", so sorting by the first
# two whitespace fields reconstructs the whole run's timeline.
#
# Usage:
#   pyntalog.sh                 # merged, filtered to the meaningful stage/warning events, in a pager
#   pyntalog.sh -a              # merged, ALL lines (includes per-opt start/finish -- verbose)
#   pyntalog.sh -g 'PATTERN'    # merged, filtered by your own egrep PATTERN
#   pyntalog.sh -o FILE         # write the merged (filtered) log to FILE instead of paging
#   pyntalog.sh -d DIR          # use DIR as the log dir (default: $PYNTA_LOG_FILE.d)
# Flags combine, e.g.  pyntalog.sh -a -o ~/run.merged.log

usage() { sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; }

# default view: stage boundaries + warnings, without the high-volume per-opt start/finish lines
filter='covdep iter|config energies|covdep select|central penalty|KILLED|WARNING'
dir="${PYNTA_LOG_FILE:+${PYNTA_LOG_FILE}.d}"
out=""

while getopts "ad:g:o:h" opt; do
  case "$opt" in
    a) filter="" ;;
    d) dir="$OPTARG" ;;
    g) filter="$OPTARG" ;;
    o) out="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) echo "try: $(basename "$0") -h" >&2; exit 2 ;;
  esac
done

if [ -z "$dir" ]; then
  echo "pyntalog: no log directory. Set PYNTA_LOG_FILE, or pass -d <dir>." >&2
  exit 2
fi
if [ ! -d "$dir" ]; then
  echo "pyntalog: '$dir' is not a directory (is PYNTA_LOG_FILE set to the same path the run used?)." >&2
  exit 2
fi

# -exec cat {} + streams every per-process file without hitting the shell arg-list limit
merged() { find "$dir" -name '*.log' -exec cat {} + | sort -k1,2; }
if [ -n "$filter" ]; then
  stream() { merged | grep -E "$filter"; }
else
  stream() { merged; }
fi

if [ -n "$out" ]; then
  stream > "$out"
  echo "wrote $(wc -l < "$out") lines to $out"
elif [ -t 1 ]; then
  stream | less
else
  stream
fi
