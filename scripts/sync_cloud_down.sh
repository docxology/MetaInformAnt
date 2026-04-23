#!/bin/bash
export PATH="/home/trim/google-cloud-sdk/bin:$PATH"
echo "Starting continuous cloud down-sync..."
while true; do
  rsync -az --exclude 'fastq' --exclude '*.fastq' --exclude '*.sra' --exclude '*.fastq.gz' --exclude '.downloads' -e "ssh -o StrictHostKeyChecking=no" metainformant-pipeline.us-central1-a.cryptoptera:/opt/MetaInformAnt/projects/hymenoptera_amalgkit/data/ projects/hymenoptera_amalgkit/data/ || echo "Rsync failed, retrying in 60s..."
  sleep 60
done
