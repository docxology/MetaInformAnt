#!/usr/bin/env bash
set -euo pipefail

# Hands-off end-to-end amalgkit runner with progress logging and setup
#
# ALL configuration stays in the YAML file provided via --config.
# This script:
#   - Runs repository setup (uv, venv, project install, amalgkit install) with a default NCBI_EMAIL
#   - Enables streaming logs unless disabled
#   - Optionally pre-checks genome presence (based on YAML's genome.dest_dir)
#   - Runs the pipeline and prints pointers to logs/manifest
#
# Usage:
#   scripts/run_amalgkit.sh \
#     --config config/amalgkit_pbarbatus.yaml \
#     [--stream] [--no-stream] [--check]

CONFIG=""
STREAM=1
CHECK_FLAG=0

usage() {
  sed -n '1,60p' "$0" | sed -n '1,40p' >&2
  exit 2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="${2:-}"; shift 2;;
    --stream) STREAM=1; shift;;
    --no-stream) STREAM=0; shift;;
    --check) CHECK_FLAG=1; shift;;
    -h|--help) usage;;
    *) echo "Unknown option: $1" >&2; usage;;
  esac
done

if [[ -z "$CONFIG" ]]; then
  echo "--config is required" >&2
  usage
fi

# 0) Setup environment (uv, venv, project, amalgkit, default email if none)
#    Pass --skip-tests to avoid running the full repo tests here.
bash scripts/setup_uv.sh --with-amalgkit --skip-tests || true

# Use the repo's virtualenv Python if available; fall back to system python3
VENV_PY=".venv/bin/python"
if [[ ! -x "$VENV_PY" ]]; then
  VENV_PY="python3"
fi

# Enable streaming logs for long-running steps unless disabled
if [[ "$STREAM" -eq 1 ]]; then
  export MI_STREAM_AMALGKIT_LOGS=1
fi

# Use the repo's Python sources
export PYTHONPATH="${PYTHONPATH:-}:$(pwd)/src"

# Ensure the venv bin is on PATH so subprocesses can find 'amalgkit'
if [[ -d ".venv/bin" ]]; then
  export PATH="$(pwd)/.venv/bin:$PATH"
fi

# Quick preflight: show amalgkit location and version/help
if command -v amalgkit >/dev/null 2>&1; then
  echo "amalgkit CLI: $(command -v amalgkit)"
  amalgkit -h | head -n 3 || true
else
  echo "amalgkit CLI not found on PATH after setup; attempting to proceed via Python module runner"
fi

# Ensure NCBI_EMAIL is exported in this shell, using value recorded by setup if present
if [[ -z "${NCBI_EMAIL:-}" && -f output/setup/ncbi_email.txt ]]; then
  export NCBI_EMAIL="$(cat output/setup/ncbi_email.txt)"
fi

# Determine genome dest_dir and effective log_dir from YAML for quick hints
GENOME_DIR=$($VENV_PY - <<PY
from metainformant.rna.workflow import load_workflow_config
cfg = load_workflow_config("$CONFIG")
print(str((cfg.genome or {}).get("dest_dir", "")))
PY
)
if [[ -n "$GENOME_DIR" ]]; then
  echo "Genome destination (from YAML): $GENOME_DIR"
  if [[ -d "$GENOME_DIR" ]]; then
    echo "Genome directory exists; runner will still verify and proceed."
  fi
fi

# Compute log dir from config
LOG_DIR_CFG=$($VENV_PY - <<PY
from metainformant.rna.workflow import load_workflow_config
cfg = load_workflow_config("$CONFIG")
print(str((cfg.log_dir or (cfg.work_dir/"logs")).resolve()))
PY
)

# Run the full workflow (hands-off)
CMD=( "$VENV_PY" -m metainformant rna run-config --config "$CONFIG" )
if [[ "$CHECK_FLAG" -eq 1 ]]; then CMD+=( --check ); fi

echo "Running: ${CMD[*]}"
"${CMD[@]}"

# Print pointers for monitoring
WORK_DIR=$($VENV_PY - <<PY
from metainformant.rna.workflow import load_workflow_config
from pathlib import Path
cfg = load_workflow_config("$CONFIG")
print(str(cfg.work_dir.resolve()))
PY
)

LOG_DIR="$LOG_DIR_CFG"

cat <<INFO

--- Monitoring ---
Manifest:   $WORK_DIR/amalgkit.manifest.jsonl
Logs:        $LOG_DIR
Tail logs:   tail -f $LOG_DIR/*.stdout.log
INFO

if [[ -f "$WORK_DIR/amalgkit.manifest.jsonl" ]]; then
  echo "\n--- MANIFEST HEAD ---"
  sed -n '1,80p' "$WORK_DIR/amalgkit.manifest.jsonl" || true
fi

# Check for metadata table existence and show a quick head
META_TSV=$($VENV_PY - <<PY
from metainformant.rna.workflow import load_workflow_config
cfg = load_workflow_config("$CONFIG")
print(str((cfg.work_dir/"metadata"/"metadata.tsv").resolve()))
PY
)

echo "\nExpected metadata table: $META_TSV"
if [[ -f "$META_TSV" ]]; then
  echo "Metadata table found. Showing first lines:"
  sed -n '1,20p' "$META_TSV" || true
else
  echo "Metadata table not found yet. If this run just started, allow time for 'amalgkit metadata' to complete."
  echo "Check logs in: $LOG_DIR (look for *metadata* logs)."
fi
