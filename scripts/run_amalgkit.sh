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
#    Try to install external CLI deps as well for an end-to-end run.
bash scripts/setup_uv.sh --with-amalgkit --with-deps --skip-tests || true

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

echo "[RUN] $(date -u +%Y-%m-%dT%H:%M:%SZ) ${CMD[*]}"
start_ts=$(date +%s)
# Allow the workflow to fail without aborting this orchestrator, so we can run verification
set +e
"${CMD[@]}"
rc=$?
set -e
end_ts=$(date +%s)
echo "[DONE] $(date -u +%Y-%m-%dT%H:%M:%SZ) code=$rc duration=$((end_ts-start_ts))s"

# After initial run, ensure all filtered/selected SRR entries have FASTQs present.
# We will retry missing ones with robust fallbacks (amalgkit with pfd=no, and as a last resort prefetch+fasterq-dump).

# Determine effective work_dir, getfastq out_dir, and threads from config
WORK_DIR=$($VENV_PY - <<PY
from metainformant.rna.workflow import load_workflow_config
cfg = load_workflow_config("$CONFIG")
print(str(cfg.work_dir.resolve()))
PY
)

THREADS=$($VENV_PY - <<PY
from metainformant.rna.workflow import load_workflow_config
cfg = load_workflow_config("$CONFIG")
print(int(cfg.threads))
PY
)

FASTQ_DIR_CFG=$($VENV_PY - <<PY
from metainformant.rna.workflow import load_workflow_config
from pathlib import Path
cfg = load_workflow_config("$CONFIG")
ps = cfg.per_step.get("getfastq", {}) if isinstance(cfg.per_step, dict) else {}
outd = ps.get("out_dir") or str(cfg.work_dir/"fastq")
print(str(Path(outd).expanduser().resolve()))
PY
)

META_FILTERED="$WORK_DIR/metadata/metadata.filtered.tissue.tsv"
META_DEFAULT="$WORK_DIR/metadata/metadata.tsv"
META=""
if [[ -f "$META_FILTERED" ]]; then META="$META_FILTERED"; elif [[ -f "$META_DEFAULT" ]]; then META="$META_DEFAULT"; fi

if [[ -n "$META" ]]; then
  echo "[VERIFY] Ensuring FASTQs present for all SRR in: $META"
  # Find the 'run' column index (SRR IDs)
  RUN_COL=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++) if($i=="run"){print i; exit}}' "$META") || RUN_COL=""
  if [[ -n "$RUN_COL" ]]; then
    echo "[VERIFY] Enumerating SRR IDs from metadata..."
    missing_any=0
    pass=1
    max_pass=2
    while [[ $pass -le $max_pass ]]; do
      echo "[PASS $pass/$max_pass] Checking and fetching missing FASTQs into: $FASTQ_DIR_CFG/getfastq/<SRR>"
      missing_this_pass=0
      # Count how many SRRs currently have FASTQs
      HAVE=$(find "$FASTQ_DIR_CFG/getfastq" -type f -name "*.fastq*" 2>/dev/null | sed -E 's@.*/(SRR[0-9]+)/.*@\1@' | sort -u | wc -l | tr -d ' ')
      TOTAL=$(tail -n +2 "$META" | awk -v c="$RUN_COL" -F"\t" 'NF>=c {print $c}' | sort -u | wc -l | tr -d ' ')
      echo "[VERIFY] Current FASTQ coverage: $HAVE/$TOTAL"
      # Launch accelerated parallel downloader to fill gaps
      JOBS=$(( THREADS > 8 ? 8 : THREADS ))
      [[ "$JOBS" -lt 2 ]] && JOBS=2
      TPJ=$(( THREADS / JOBS ))
      [[ "$TPJ" -lt 2 ]] && TPJ=2
      if [[ -x scripts/force_fasterq_parallel.sh ]]; then
        echo "[ACCEL] Starting parallel downloader: jobs=$JOBS threads-per-job=$TPJ"
        bash scripts/force_fasterq_parallel.sh --config "$CONFIG" --jobs "$JOBS" --threads-per-job "$TPJ" | sed -n '1,60p' || true
      else
        echo "[WARN] scripts/force_fasterq_parallel.sh not found or not executable; skipping parallel acceleration"
      fi
      # Recompute coverage to determine if another pass is needed
      HAVE2=$(find "$FASTQ_DIR_CFG/getfastq" -type f -name "*.fastq*" 2>/dev/null | sed -E 's@.*/(SRR[0-9]+)/.*@\1@' | sort -u | wc -l | tr -d ' ')
      if [[ "$HAVE2" -gt "$HAVE" ]]; then
        missing_any=1
        missing_this_pass=1
      fi
      if [[ $missing_this_pass -eq 0 ]]; then
        break
      fi
      pass=$((pass+1))
    done
    if [[ $missing_any -eq 1 ]]; then
      echo "[VERIFY] Completed retries. Remaining missing (if any) can be inspected under: $FASTQ_DIR_CFG"
      echo "[CONTINUE] Launching a second pipeline pass to proceed into quant/merge with newly downloaded FASTQs"
      set +e
      "${CMD[@]}"
      rc2=$?
      set -e
      echo "[CONTINUE DONE] code=$rc2"
    else
      echo "[VERIFY] All SRR FASTQs present."
      echo "[CONTINUE] All reads present; running pipeline to finalize downstream steps"
      set +e
      "${CMD[@]}"
      rc2=$?
      set -e
      echo "[CONTINUE DONE] code=$rc2"
    fi
  else
    echo "[WARN] Could not locate 'run' column in $META; skipping FASTQ verification"
  fi
else
  echo "[INFO] Metadata table not found; skipping FASTQ verification loop"
fi

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
