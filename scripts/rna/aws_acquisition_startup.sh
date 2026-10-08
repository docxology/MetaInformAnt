#!/bin/bash
# Rendered by the completion controller with immutable input object bindings.
set -Eeuo pipefail
BUCKET=@@BUCKET@@
REGION=@@REGION@@
COHORT=@@COHORT@@
SPECIES=@@SPECIES@@
SOURCE_KEY=@@SOURCE_KEY@@
SOURCE_SHA=@@SOURCE_SHA@@
INPUT_KEY=@@INPUT_KEY@@
INPUT_SHA=@@INPUT_SHA@@
JOB_PREFIX=@@JOB_PREFIX@@
LIMIT_SECONDS=@@LIMIT_SECONDS@@
RAW_BYTES=@@RAW_BYTES@@
mkdir -p /mnt/completion /opt/amalgkit /mnt/snapshot
export TMPDIR=/mnt/completion/tmp
mkdir -p "$TMPDIR"
exec > >(tee -a /mnt/completion/startup.log) 2>&1
export AWS_DEFAULT_REGION="$REGION"
RUN_ID=unknown

finish() {
  code=$?
  trap - EXIT
  set +e
  printf '{"status":"finished","worker_code":%s,"instance_id":"%s"}\n' "$code" "$RUN_ID" > /mnt/completion/status.json
  aws s3 cp /mnt/completion/startup.log "s3://$BUCKET/$JOB_PREFIX/startup.log" --region "$REGION"
  log_code=$?
  if [ "$log_code" -ne 0 ]; then printf 'Final log upload failed: %s\n' "$log_code" >&2; fi
  aws s3 cp /mnt/completion/status.json "s3://$BUCKET/$JOB_PREFIX/status.json" --region "$REGION"
  status_code=$?
  if [ "$status_code" -ne 0 ]; then printf 'Final status upload failed: %s\n' "$status_code" >&2; fi
  shutdown -h now
  shutdown_code=$?
  if [ "$shutdown_code" -ne 0 ]; then
    printf 'OS shutdown failed: %s; requesting EC2 termination\n' "$shutdown_code" >&2
    aws ec2 terminate-instances --region "$REGION" --instance-ids "$RUN_ID"
  fi
  exit "$code"
}
trap finish EXIT
(sleep "$((LIMIT_SECONDS + 600))"; shutdown -h now) &
TOKEN=$(curl -fsS -m 15 -X PUT http://169.254.169.254/latest/api/token -H 'X-aws-ec2-metadata-token-ttl-seconds: 21600')
RUN_ID=$(curl -fsS -m 15 -H "X-aws-ec2-metadata-token: $TOKEN" http://169.254.169.254/latest/meta-data/instance-id)
aws --version
dnf install -y --setopt=install_weak_deps=False git pigz libgomp zlib bzip2 xz ncurses
if ! command -v uv; then
  curl -fLsS https://astral.sh/uv/install.sh | env UV_INSTALL_DIR=/usr/local/bin sh
fi
uv python install 3.12
curl -fL --retry 5 https://github.com/pachterlab/kallisto/releases/download/v0.52.0/kallisto_linux-v0.52.0.tar.gz -o "$TMPDIR/kallisto.tar.gz"
tar -xzf "$TMPDIR/kallisto.tar.gz" -C "$TMPDIR"
install -m 0755 "$TMPDIR/kallisto/kallisto" /usr/local/bin/kallisto
kallisto version
curl -fL --retry 5 https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.4.1/sratoolkit.3.4.1-ubuntu64.tar.gz -o "$TMPDIR/sra.tar.gz"
tar -xzf "$TMPDIR/sra.tar.gz" -C /opt
export PATH="/opt/sratoolkit.3.4.1-ubuntu64/bin:/usr/local/bin:$PATH"
fasterq-dump --version
aws s3 cp "s3://$BUCKET/$SOURCE_KEY" /mnt/completion/source.tar --region "$REGION"
printf '%s  /mnt/completion/source.tar\n' "$SOURCE_SHA" | sha256sum --check
aws s3 cp "s3://$BUCKET/$INPUT_KEY" /mnt/completion/inputs.tar --region "$REGION"
printf '%s  /mnt/completion/inputs.tar\n' "$INPUT_SHA" | sha256sum --check
tar -xf /mnt/completion/source.tar -C /opt/amalgkit
tar -xf /mnt/completion/inputs.tar -C /mnt/snapshot
export AMALGKIT_DATA_ROOT="/mnt/amalgkit/$SPECIES"
mkdir -p "$AMALGKIT_DATA_ROOT"
cp -a /mnt/snapshot/data/. "$AMALGKIT_DATA_ROOT/"
mkdir -p "$AMALGKIT_DATA_ROOT/$SPECIES/genome/index"
cp "$AMALGKIT_DATA_ROOT/$SPECIES/work/index/"*.idx "$AMALGKIT_DATA_ROOT/$SPECIES/genome/index/"
export AMALGKIT_DURABLE_BUCKET="$BUCKET" AMALGKIT_DURABLE_COHORT="$COHORT"
export AMALGKIT_CLOUD_MAX_RAW_BYTES="$RAW_BYTES"
export AMALGKIT_RECLAIM_RAW_AFTER_QUANT=yes
CPU_THREADS=$(getconf _NPROCESSORS_ONLN)
export AMALGKIT_PIPELINE_QUANT_SLOTS=@@QUANT_SLOTS@@ AMALGKIT_PIPELINE_FASTQ_SLOTS=@@FASTQ_SLOTS@@
export AMALGKIT_PIPELINE_ENA_RETRY_DELAY_SECONDS=45
export AMALGKIT_MIN_EXTERNAL_FREE_GB=64 AMALGKIT_MIN_SYSTEM_FREE_GB=32
export AMALGKIT_PIPELINE_DOWNLOAD_TIMEOUT_SECONDS="$LIMIT_SECONDS" AMALGKIT_PIPELINE_FASTQ_TIMEOUT_SECONDS="$LIMIT_SECONDS"
export AMALGKIT_PIPELINE_QUANT_TIMEOUT_SECONDS="$LIMIT_SECONDS"
upload_progress() {
  while sleep 60; do
    # Diagnostic telemetry is best-effort; quant receipts use verified writes.
    if ! aws s3 cp /mnt/completion/startup.log "s3://$BUCKET/$JOB_PREFIX/startup.log" \
      --region "$REGION" --cli-connect-timeout 10 --cli-read-timeout 30; then
      printf 'Rolling diagnostic log upload failed\n' >&2
    fi
  done
}
upload_progress &
cd /opt/amalgkit
export PYTHONPATH="$PWD/src"
timeout --signal=TERM --kill-after=120 "$LIMIT_SECONDS" \
  uv run --frozen --no-dev --extra rna --extra aws --python 3.12 python \
  scripts/rna/acquisition_worker.py \
  --manifest /mnt/snapshot/manifest.jsonl --data-root "$AMALGKIT_DATA_ROOT" \
  --config-dir /mnt/snapshot/config/amalgkit \
  --workers @@WORKERS@@ --threads @@THREADS@@ --quant-slots @@QUANT_SLOTS@@ --fastq-slots @@FASTQ_SLOTS@@ \
  --max-in-flight @@MAX_IN_FLIGHT@@ --fastq-threads @@FASTQ_THREADS@@ --compression-threads @@COMPRESSION_THREADS@@ --compression-level @@COMPRESSION_LEVEL@@ --validation-slots @@VALIDATION_SLOTS@@
