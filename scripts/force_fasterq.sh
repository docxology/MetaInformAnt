#!/usr/bin/env bash
set -euo pipefail

# Force-download missing SRR FASTQs via sra-tools fasterq-dump.
# - Reads config to locate work_dir, per-step getfastq out_dir, and threads
# - Enumerates SRR IDs from metadata (prefers filtered.tissue.tsv)
# - For each SRR missing FASTQs, runs fasterq-dump with a temp directory, then compresses outputs
# - Logs to output/amalgkit/pbarbatus/logs/getfastq.force_sratools.log
#
# Usage:
#   scripts/force_fasterq.sh --config config/amalgkit_pbarbatus.yaml [--threads N] [--parallel M]

CONFIG=""
THREADS_OVERRIDE=""
PARALLEL_JOBS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="${2:-}"; shift 2;;
    --threads) THREADS_OVERRIDE="${2:-}"; shift 2;;
    --parallel) PARALLEL_JOBS="${2:-}"; shift 2;;
    -h|--help)
      echo "Usage: $0 --config <yaml> [--threads N]"; exit 2;;
    *) echo "Unknown option: $1" >&2; exit 2;;
  esac
done

if [[ -z "$CONFIG" ]]; then
  echo "--config is required" >&2
  exit 2
fi

if [[ -x ".venv/bin/python" ]]; then
  PY=".venv/bin/python"
else
  PY="python3"
fi

# Resolve paths from config
WORK_DIR=$($PY - <<PY
from metainformant.rna.workflow import load_workflow_config
cfg = load_workflow_config("$CONFIG")
print(str(cfg.work_dir.resolve()))
PY
)

FASTQ_DIR=$($PY - <<PY
from metainformant.rna.workflow import load_workflow_config
from pathlib import Path
cfg = load_workflow_config("$CONFIG")
ps = cfg.per_step.get("getfastq", {}) if isinstance(cfg.per_step, dict) else {}
outd = ps.get("out_dir") or str(cfg.work_dir/"fastq")
print(str(Path(outd).expanduser().resolve()))
PY
)

CFG_THREADS=$($PY - <<PY
from metainformant.rna.workflow import load_workflow_config
cfg = load_workflow_config("$CONFIG")
print(int(cfg.threads))
PY
)

if [[ -n "$THREADS_OVERRIDE" ]]; then
  THREADS="$THREADS_OVERRIDE"
else
  THREADS="$CFG_THREADS"
fi

LOG="output/amalgkit/pbarbatus/logs/getfastq.force_sratools.log"
mkdir -p "$(dirname "$LOG")" "$FASTQ_DIR/getfastq" "output/amalgkit/_tmp/fasterq"

# Select metadata file
META=""
if [[ -f "$WORK_DIR/metadata/metadata.filtered.tissue.tsv" ]]; then
  META="$WORK_DIR/metadata/metadata.filtered.tissue.tsv"
elif [[ -f "$WORK_DIR/metadata/metadata.tsv" ]]; then
  META="$WORK_DIR/metadata/metadata.tsv"
fi

{
  echo "[START $(date -u +%Y-%m-%dT%H:%M:%SZ)] force sra-tools loop (threads=$THREADS)"
  if [[ -z "$META" ]]; then
    echo "[ERROR] No metadata table under $WORK_DIR/metadata"
    echo "[END $(date -u +%Y-%m-%dT%H:%M:%SZ)]"
    exit 0
  fi

  # Extract SRR list from metadata's 'run' column
  RUN_COL=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++) if($i=="run"){print i; exit}}' "$META")
  if [[ -z "$RUN_COL" ]]; then
    echo "[ERROR] No run column in $META"
    echo "[END $(date -u +%Y-%m-%dT%H:%M:%SZ)]"
    exit 0
  fi

  # Build SRR list
  tail -n +2 "$META" | awk -v c="$RUN_COL" -F"\t" 'NF>=c {print $c}' | sort -u > "$WORK_DIR/metadata/srr_list.txt"
  TOTAL=$(wc -l < "$WORK_DIR/metadata/srr_list.txt" | tr -d ' ')
  PROG_DIR="$WORK_DIR/progress_force"
  mkdir -p "$PROG_DIR"

  # Define worker function (compatible with macOS bash without wait -n)
  process_one_srr() {
    local SRR="$1"
    local SRR_DIR="$FASTQ_DIR/getfastq/$SRR"
    local TMP_DIR="output/amalgkit/_tmp/fasterq/$SRR"
    local f1gz="$SRR_DIR/${SRR}_1.fastq.gz"; local f2gz="$SRR_DIR/${SRR}_2.fastq.gz"; local fsing="$SRR_DIR/${SRR}.fastq.gz"
    touch "$PROG_DIR/${SRR}.start" 2>/dev/null || true
    if [[ -f "$f1gz" || -f "$f2gz" || -f "$fsing" ]]; then
      echo "[SKIP] $SRR already has FASTQ"
      ln -sf "${SRR}.start" "$PROG_DIR/${SRR}.done" 2>/dev/null || true
      return 0
    fi
    mkdir -p "$SRR_DIR" "$TMP_DIR"
    echo "[FETCH] $SRR via prefetch+fasterq-dump"
    # Ensure sra-tools present
    if ! command -v fasterq-dump >/dev/null 2>&1; then
      echo "[ERROR] fasterq-dump not on PATH"
      return 0
    fi
    # Prefetch locally to avoid network streaming errors
    if command -v prefetch >/dev/null 2>&1; then
      local pattempt
      for pattempt in 1 2; do
        stdbuf -oL -eL prefetch --progress --output-directory "$SRR_DIR" "$SRR" 2>&1 \
          | sed -u "s/^/[SRR $SRR] prefetch /" \
          >> "$LOG" && break || true
        echo "[RETRY] $SRR prefetch attempt $pattempt failed; retrying" >> "$LOG"
        sleep 3
      done
    fi
    local attempt
    for attempt in 1 2; do
      # Prefer local .sra if present
      local INPUT_ACC
      if [[ -f "$SRR_DIR/$SRR.sra" ]]; then INPUT_ACC="$SRR_DIR/$SRR.sra"; else INPUT_ACC="$SRR"; fi
      stdbuf -oL -eL fasterq-dump -e "$THREADS" -p --split-files -O "$SRR_DIR" -t "$TMP_DIR" "$INPUT_ACC" 2>&1 \
        | sed -u "s/^/[SRR $SRR] /" \
        >> "$LOG" && break || true
      echo "[RETRY] $SRR attempt $attempt failed; retrying"
      sleep 3
    done
    if ls "$SRR_DIR"/${SRR}_*.fastq >/dev/null 2>&1 || [[ -f "$SRR_DIR/${SRR}.fastq" ]]; then
      if command -v pigz >/dev/null 2>&1; then
        pigz -f "$SRR_DIR"/${SRR}_*.fastq 2>/dev/null | sed -u "s/^/[SRR $SRR] /" >> "$LOG" || true
        [[ -f "$SRR_DIR/${SRR}.fastq" ]] && pigz -f "$SRR_DIR/${SRR}.fastq" 2>/dev/null | sed -u "s/^/[SRR $SRR] /" >> "$LOG" || true
      else
        gzip -f "$SRR_DIR"/${SRR}_*.fastq 2>/dev/null | sed -u "s/^/[SRR $SRR] /" >> "$LOG" || true
        [[ -f "$SRR_DIR/${SRR}.fastq" ]] && gzip -f "$SRR_DIR/${SRR}.fastq" 2>/dev/null | sed -u "s/^/[SRR $SRR] /" >> "$LOG" || true
      fi
      echo "[DONE] $SRR"
      # Remove local .sra to save space
      rm -f "$SRR_DIR/$SRR.sra" 2>/dev/null || true
      ln -sf "${SRR}.start" "$PROG_DIR/${SRR}.done" 2>/dev/null || true
      # Progress summary
      local DONE_COUNT
      DONE_COUNT=$(ls "$PROG_DIR"/*.done 2>/dev/null | wc -l | tr -d ' ' || echo 0)
      local PCT
      if [[ "$TOTAL" -gt 0 ]]; then
        PCT=$(( 100 * DONE_COUNT / TOTAL ))
      else
        PCT=0
      fi
      echo "[PROGRESS] $DONE_COUNT/$TOTAL ($PCT%) FASTQs ready"
    else
      echo "[WARN] No FASTQ files produced for $SRR"
    fi
    rm -rf "$TMP_DIR" 2>/dev/null || true
  }

  export -f process_one_srr || true

  # Throttled parallel execution
  active=0
  while read -r SRR; do
    [[ -z "$SRR" ]] && continue
    process_one_srr "$SRR" &
    active=$((active+1))
    # Throttle when reaching PARALLEL_JOBS
    while :; do
      running=$(jobs -rp | wc -l | tr -d ' ')
      if [[ "$running" -lt "$PARALLEL_JOBS" ]]; then break; fi
      sleep 2
    done
  done < "$WORK_DIR/metadata/srr_list.txt"
  wait || true
  echo "[END $(date -u +%Y-%m-%dT%H:%M:%SZ)]"
} >> "$LOG" 2>&1


