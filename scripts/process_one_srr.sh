#!/usr/bin/env bash
set -euo pipefail

# Usage: scripts/process_one_srr.sh <SRR> <FASTQ_DIR> <THREADS_PER_JOB> <TMP_ROOT> <LOG_FILE>
# Produces FASTQs under <FASTQ_DIR>/getfastq/<SRR>/ and logs progress to <LOG_FILE>

SRR="${1:?SRR required}"
FASTQ_DIR="${2:?FASTQ_DIR required}"
TPJ="${3:?THREADS_PER_JOB required}"
TMP_ROOT="${4:?TMP_ROOT required}"
LOG_FILE="${5:?LOG_FILE required}"

SRR_DIR="$FASTQ_DIR/getfastq/$SRR"
TMP_DIR="$TMP_ROOT/$SRR"
mkdir -p "$SRR_DIR" "$TMP_DIR" "$(dirname "$LOG_FILE")"

log() { echo "[SRR $SRR] $*" >> "$LOG_FILE"; }

f1gz="$SRR_DIR/${SRR}_1.fastq.gz"; f2gz="$SRR_DIR/${SRR}_2.fastq.gz"; fsing="$SRR_DIR/${SRR}.fastq.gz"
if [[ -f "$f1gz" || -f "$f2gz" || -f "$fsing" ]]; then
  log "skip: already has FASTQ"
  exit 0
fi

# 1) Try AWS ODP S3 .sra direct (often faster than prefetch)
AWS_URL="https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRR}/${SRR}"
log "AWS ODP check: $AWS_URL"
if command -v curl >/dev/null 2>&1 && curl -sfI "$AWS_URL" >/dev/null 2>&1; then
  log "AWS ODP download start"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x 16 -s 16 -k 1M -d "$SRR_DIR" -o "${SRR}.sra" "$AWS_URL" >> "$LOG_FILE" 2>&1 || true
  else
    curl -fSL -o "$SRR_DIR/${SRR}.sra" "$AWS_URL" >> "$LOG_FILE" 2>&1 || true
  fi
  if [[ -f "$SRR_DIR/${SRR}.sra" ]]; then
    INPUT="$SRR_DIR/${SRR}.sra"
    log "AWS ODP download done"
    # Proceed to fasterq-dump below
  fi
fi

# 2) Try ENA direct FASTQ links (often much faster than SRA)
log "ENA filereport query"
ENA_TSV="$TMP_DIR/ena.tsv"
mkdir -p "$TMP_DIR"
if command -v curl >/dev/null 2>&1; then
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp,fastq_http,fastq_galaxy,submitted_ftp,submitted_http&format=tsv&download=false" -o "$ENA_TSV" 2>>"$LOG_FILE" || true
  if [[ -s "$ENA_TSV" ]]; then
    FASTQ_HTTP=$(awk -F"\t" 'NR==2{print $3}' "$ENA_TSV" 2>/dev/null || true)
    FASTQ_FTP=$(awk -F"\t" 'NR==2{print $2}' "$ENA_TSV" 2>/dev/null || true)
    DL_LINKS=()
    if [[ -n "${FASTQ_HTTP:-}" ]]; then
      IFS=';' read -r -a LINKS <<< "$FASTQ_HTTP"
      for L in "${LINKS[@]}"; do
        [[ -n "$L" ]] && DL_LINKS+=("$L")
      done
    elif [[ -n "${FASTQ_FTP:-}" ]]; then
      IFS=';' read -r -a LINKS <<< "$FASTQ_FTP"
      for L in "${LINKS[@]}"; do
        [[ -n "$L" ]] && DL_LINKS+=("https://$L")
      done
    fi
    if [[ ${#DL_LINKS[@]} -gt 0 ]]; then
      log "ENA download start (${#DL_LINKS[@]} files)"
      if command -v aria2c >/dev/null 2>&1; then
        # Multi-connection per file
        for URL in "${DL_LINKS[@]}"; do
          FN=$(basename "$URL")
          aria2c -x 16 -s 16 -k 1M -d "$SRR_DIR" -o "$FN" "$URL" >> "$LOG_FILE" 2>&1 || true
        done
      else
        # Fallback: parallel curls
        pids=()
        for URL in "${DL_LINKS[@]}"; do
          ( cd "$SRR_DIR" && curl -fSL -O "$URL" ) >> "$LOG_FILE" 2>&1 &
          pids+=("$!")
        done
        for p in "${pids[@]}"; do wait "$p" || true; done
      fi
      if [[ -f "$f1gz" || -f "$f2gz" || -f "$fsing" ]]; then
        log "ENA download done"
        rm -rf "$TMP_DIR" 2>/dev/null || true
        log "done"
        exit 0
      else
        log "ENA not available for $SRR or download failed; falling back to SRA"
      fi
    fi
  fi
fi

# 3) SRA path: prefetch then fasterq-dump (if not already downloaded via AWS)
log "prefetch start"
if command -v prefetch >/dev/null 2>&1; then
  prefetch --progress --output-directory "$SRR_DIR" "$SRR" >> "$LOG_FILE" 2>&1 || true
fi

INPUT="$SRR"
if [[ -f "$SRR_DIR/$SRR.sra" ]]; then INPUT="$SRR_DIR/$SRR.sra"; fi

log "fasterq-dump start (threads=$TPJ)"
attempt=1
until [[ $attempt -gt 2 ]]; do
  if fasterq-dump -e "$TPJ" -p --split-files -O "$SRR_DIR" -t "$TMP_DIR" "$INPUT" >> "$LOG_FILE" 2>&1; then
    break
  fi
  log "retry: fasterq-dump attempt $attempt failed"
  attempt=$((attempt+1))
  sleep 2
done

produced=0
if ls "$SRR_DIR"/${SRR}_*.fastq >/dev/null 2>&1 || [[ -f "$SRR_DIR/${SRR}.fastq" ]]; then
  produced=1
fi

if [[ "$produced" -eq 1 ]]; then
  log "compressing"
  if command -v pigz >/dev/null 2>&1; then
    pigz -f "$SRR_DIR"/${SRR}_*.fastq 2>/dev/null || true
    [[ -f "$SRR_DIR/${SRR}.fastq" ]] && pigz -f "$SRR_DIR/${SRR}.fastq" 2>/dev/null || true
  else
    gzip -f "$SRR_DIR"/${SRR}_*.fastq 2>/dev/null || true
    [[ -f "$SRR_DIR/${SRR}.fastq" ]] && gzip -f "$SRR_DIR/${SRR}.fastq" 2>/dev/null || true
  fi
  rm -rf "$TMP_DIR" 2>/dev/null || true
  rm -f "$SRR_DIR/$SRR.sra" 2>/dev/null || true
  log "done"
  exit 0
else
  log "warn: no FASTQ output"
  rm -rf "$TMP_DIR" 2>/dev/null || true
  exit 0
fi


