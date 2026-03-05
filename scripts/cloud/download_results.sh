#!/usr/bin/env bash
# download_results.sh — Download amalgkit quant results from GCP VM
#
# Usage: bash scripts/cloud/download_results.sh [--output DIR] [--species SPECIES]
#
# Downloads quant directories from the Docker container on GCP to local machine.
# Only downloads the small quant output files (~2 MB/sample), NOT raw FASTQs.
#
# Estimated sizes:
#   - 1,700 samples → ~3.5 GB total quant output
#   - 5,000 samples → ~10 GB total quant output

set -euo pipefail

PROJECT="${GCP_PROJECT:-cryptoptera}"
ZONE="${GCP_ZONE:-us-central1-a}"
VM_NAME="${GCP_VM:-metainformant-pipeline}"
CONTAINER="${CONTAINER_NAME:-metainformant-pipeline}"
LOCAL_OUTPUT="output/amalgkit"
SPECIES=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --output)
      LOCAL_OUTPUT="$2"
      shift 2
      ;;
    --species)
      SPECIES="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

echo "╔═══════════════════════════════════════════════════╗"
echo "║  Downloading amalgkit results from GCP            ║"
echo "╚═══════════════════════════════════════════════════╝"
echo ""
echo "  Project:   $PROJECT"
echo "  VM:        $VM_NAME"
echo "  Container: $CONTAINER"
echo "  Local dir: $LOCAL_OUTPUT"
echo ""

# Step 1: Copy quant results from Docker container to VM filesystem
echo "📦 Step 1/3: Copying quant output from container to VM /tmp/quant_export..."
gcloud compute ssh "$VM_NAME" \
    --project "$PROJECT" \
    --zone "$ZONE" \
    --quiet \
    --command="
        sudo rm -rf /tmp/quant_export
        sudo mkdir -p /tmp/quant_export

        # Copy all quant directories from container
        for sp_dir in \$(sudo docker exec $CONTAINER ls output/amalgkit/ 2>/dev/null); do
            if [ -n "$SPECIES" ] && [ "\$sp_dir" != "$SPECIES" ]; then
                continue
            fi
            if sudo docker exec $CONTAINER test -d output/amalgkit/\$sp_dir/work/quant; then
                count=\$(sudo docker exec $CONTAINER ls output/amalgkit/\$sp_dir/work/quant 2>/dev/null | wc -l)
                if [ \"\$count\" -gt 0 ]; then
                    sudo mkdir -p /tmp/quant_export/\$sp_dir/work/quant
                    sudo docker cp $CONTAINER:/app/output/amalgkit/\$sp_dir/work/quant /tmp/quant_export/\$sp_dir/work/
                    echo \"  ✓ \$sp_dir: \$count samples\"
                fi
            fi
        done

        # Also copy progress DB and any merged results
        sudo docker cp $CONTAINER:/app/output/amalgkit/pipeline_progress.db /tmp/quant_export/ 2>/dev/null || true
        for sp_dir in \$(sudo docker exec $CONTAINER ls output/amalgkit/ 2>/dev/null); do
            if [ -n "$SPECIES" ] && [ "\$sp_dir" != "$SPECIES" ]; then
                continue
            fi
            if sudo docker exec $CONTAINER test -d output/amalgkit/\$sp_dir/merged; then
                sudo mkdir -p /tmp/quant_export/\$sp_dir/
                sudo docker cp $CONTAINER:/app/output/amalgkit/\$sp_dir/merged /tmp/quant_export/\$sp_dir/ 2>/dev/null || true
                echo \"  ✓ \$sp_dir/merged\"
            fi
        done

        # Fix permissions for SCP
        sudo chmod -R 755 /tmp/quant_export
        sudo chown -R \$(whoami) /tmp/quant_export

        echo ''
        echo \"Total export size:\"
        du -sh /tmp/quant_export
    "

echo ""
echo "📥 Step 2/3: Downloading from VM to local machine..."
mkdir -p "$LOCAL_OUTPUT"

gcloud compute scp --recurse \
    --project "$PROJECT" \
    --zone "$ZONE" \
    "${VM_NAME}:/tmp/quant_export/*" \
    "$LOCAL_OUTPUT/"

echo ""
echo "🧹 Step 3/3: Cleaning up VM temp files..."
gcloud compute ssh "$VM_NAME" \
    --project "$PROJECT" \
    --zone "$ZONE" \
    --quiet \
    --command="sudo rm -rf /tmp/quant_export"

echo ""
echo "╔═══════════════════════════════════════════════════╗"
echo "║  ✅ Download complete!                            ║"
echo "╚═══════════════════════════════════════════════════╝"
echo ""
echo "Results saved to: $LOCAL_OUTPUT/"
echo ""
echo "Contents:"
du -sh "$LOCAL_OUTPUT"/*/work/quant 2>/dev/null | sort -rh || echo "  (no quant dirs found)"
echo ""
echo "Total:"
du -sh "$LOCAL_OUTPUT"
