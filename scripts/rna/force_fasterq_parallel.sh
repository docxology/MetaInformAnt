#!/usr/bin/env bash
set -euo pipefail

# Parallel SRR downloader using prefetch + fasterq-dump per SRR.
# Usage:
#   scripts/force_fasterq_parallel.sh --config config/amalgkit_pbarbatus.yaml \
#     [--jobs 4] [--threads-per-job 4]

CONFIG=""; JOBS=4; TPJ=4
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="${2:-}"; shift 2;;
    --jobs) JOBS="${2:-}"; shift 2;;
    --threads-per-job) TPJ="${2:-}"; shift 2;;
    -h|--help) echo "Usage: $0 --config <yaml> [--jobs N] [--threads-per-job M]"; exit 2;;
    *) echo "Unknown option: $1" >&2; exit 2;;
  esac
done

[[ -z "$CONFIG" ]] && { echo "--config is required" >&2; exit 2; }

if [[ -x .venv/bin/python ]]; then PY=.venv/bin/python; else PY=python3; fi

WORK_DIR=$($PY - <<PY
from metainformant.rna.workflow import load_workflow_config
cfg=load_workflow_config("$CONFIG")
print(str(cfg.work_dir.resolve()))
PY
)

FASTQ_DIR=$($PY - <<PY
from metainformant.rna.workflow import load_workflow_config
from pathlib import Path
cfg=load_workflow_config("$CONFIG")
ps=cfg.per_step.get("getfastq", {}) if isinstance(cfg.per_step, dict) else {}
outd = ps.get("out_dir") or str(cfg.work_dir/"fastq")
print(str(Path(outd).expanduser().resolve()))
PY
)

THREADS=$($PY - <<PY
from metainformant.rna.workflow import load_workflow_config
cfg=load_workflow_config("$CONFIG")
print(int(cfg.threads))
PY
)

LOG="output/amalgkit/pbarbatus/logs/getfastq.force_parallel.log"
TMP_ROOT="output/amalgkit/_tmp/fasterq"
mkdir -p "$(dirname "$LOG")" "$FASTQ_DIR/getfastq" "$TMP_ROOT" "$WORK_DIR/metadata"

META="$WORK_DIR/metadata/metadata.filtered.tissue.tsv"
[[ -f "$META" ]] || META="$WORK_DIR/metadata/metadata.tsv"
[[ -f "$META" ]] || { echo "[ERROR] No metadata table under $WORK_DIR/metadata" | tee -a "$LOG"; exit 0; }

RUN_COL=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++) if($i=="run"){print i; exit}}' "$META")
[[ -n "$RUN_COL" ]] || { echo "[ERROR] No run column in $META" | tee -a "$LOG"; exit 0; }

SRR_LIST_FILE="$WORK_DIR/metadata/srr_list.txt"
tail -n +2 "$META" | awk -v c="$RUN_COL" -F"\t" 'NF>=c {print $c}' | sort -u > "$SRR_LIST_FILE"
TOTAL=$(wc -l < "$SRR_LIST_FILE" | tr -d ' ')

echo "[START $(date -u +%Y-%m-%dT%H:%M:%SZ)] parallel getfastq jobs=$JOBS tpj=$TPJ total=$TOTAL" | tee -a "$LOG"

active=0
while read -r SRR; do
  [[ -z "$SRR" ]] && continue
  bash scripts/process_one_srr.sh "$SRR" "$FASTQ_DIR" "$TPJ" "$TMP_ROOT" "$LOG" &
  while :; do
    running=$(jobs -rp | wc -l | tr -d ' ')
    if [[ "$running" -lt "$JOBS" ]]; then break; fi
    done_count=$(find "$FASTQ_DIR/getfastq" -type f -name "*.fastq.gz" 2>/dev/null | sed -E 's@.*/(SRR[0-9]+)/.*@\1@' | sort -u | wc -l | tr -d ' ')
    pct=0; [[ "$TOTAL" -gt 0 ]] && pct=$((100*done_count/TOTAL))
    echo "[PROGRESS $(date -u +%H:%M:%SZ)] $done_count/$TOTAL ($pct%)" | tee -a "$LOG"
    sleep 2
  done
done < "$SRR_LIST_FILE"

wait || true

done_count=$(find "$FASTQ_DIR/getfastq" -type f -name "*.fastq.gz" 2>/dev/null | sed -E 's@.*/(SRR[0-9]+)/.*@\1@' | sort -u | wc -l | tr -d ' ')
pct=0; [[ "$TOTAL" -gt 0 ]] && pct=$((100*done_count/TOTAL))
echo "[END $(date -u +%Y-%m-%dT%H:%M:%SZ)] completed=$done_count/$TOTAL ($pct%)" | tee -a "$LOG"
