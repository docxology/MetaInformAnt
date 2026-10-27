#!/bin/bash
# Host aliases come from `gcloud compute config-ssh`. Set GOOGLE_CLOUD_SDK_BIN to put a
# non-standard gcloud install on PATH (needed only for IAP-tunnelled ssh).
if [ -n "${GOOGLE_CLOUD_SDK_BIN:-}" ]; then
  export PATH="$GOOGLE_CLOUD_SDK_BIN:$PATH"
fi
# accept-new trusts a host on first contact and refuses a changed key afterwards.
echo "Starting continuous cloud down-sync..."
while true; do
  rsync -az --exclude 'fastq' --exclude '*.fastq' --exclude '*.sra' --exclude '*.fastq.gz' --exclude '.downloads' -e "ssh -o StrictHostKeyChecking=accept-new" metainformant-pipeline.us-central1-a.cryptoptera:/opt/MetaInformAnt/projects/hymenoptera_amalgkit/data/ projects/hymenoptera_amalgkit/data/ || echo "Rsync failed, retrying in 60s..."
  sleep 60
done
