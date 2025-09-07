#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ------------------------
# Defaults for SBATCH
# ------------------------
ACCOUNT="hpc2n2024-161"
JOB_NAME="${PWD##*/}"
NTASKS=1
CPUS_PER_TASK=1
CONSTRAINT="skylake"
TIME="0-06:00:00"
LOG_DIR="logs"
OUT_PATTERN="${LOG_DIR}/%x-%j.out"
ERR_PATTERN="${LOG_DIR}/%x-%j.err"
DEPENDENCY=""

usage() {
  cat <<EOF >&2
Usage: $0 [options] script_path [script_args...]
Options:
  -A ACCOUNT        SLURM account (default: $ACCOUNT)
  -J JOB_NAME       job name (default: current dir: $JOB_NAME)
  -n NTASKS         number of tasks (default: $NTASKS)
  -c CPUS           cpus-per-task (default: $CPUS_PER_TASK)
  -C CONSTRAINT     SLURM constraint (default: $CONSTRAINT)
  -t TIME           time limit, D-HH:MM:SS (default: $TIME)
  -o OUT_PATTERN    stdout file pattern (default: $OUT_PATTERN)
  -e ERR_PATTERN    stderr file pattern (default: $ERR_PATTERN)
  -d DEP            dependency string for sbatch, e.g. "afterok:123:456"
  -h                show this help and exit

Note: prints the submitted JobID to stdout.
EOF
  exit 1
}

while getopts "A:J:n:c:C:t:o:e:d:h" opt; do
  case "$opt" in
    A) ACCOUNT="$OPTARG" ;;
    J) JOB_NAME="$OPTARG" ;;
    n) NTASKS="$OPTARG" ;;
    c) CPUS_PER_TASK="$OPTARG" ;;
    C) CONSTRAINT="$OPTARG" ;;
    t) TIME="$OPTARG" ;;
    o) OUT_PATTERN="$OPTARG" ;;
    e) ERR_PATTERN="$OPTARG" ;;
    d) DEPENDENCY="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done
shift $((OPTIND-1))

if [ "$#" -lt 1 ]; then usage; fi
SCRIPT="$1"; shift
SCRIPT_ARGS=("$@")

mkdir -p "$LOG_DIR"

SBATCH_ARGS=(
  -A "$ACCOUNT"
  -J "$JOB_NAME"
  -n "$NTASKS"
  -c "$CPUS_PER_TASK"
  --constraint="$CONSTRAINT"
  -t "$TIME"
  --output="$OUT_PATTERN"
  --error="$ERR_PATTERN"
  --parsable
)
if [[ -n "$DEPENDENCY" ]]; then
  SBATCH_ARGS+=( --dependency="$DEPENDENCY" )
fi

# Submit and print just the JobID
JOBID=$(sbatch "${SBATCH_ARGS[@]}" "$SCRIPT" "${SCRIPT_ARGS[@]}")
echo "$JOBID"
