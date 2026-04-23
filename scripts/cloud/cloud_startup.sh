#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────
# MetaInformAnt — GCP VM Startup Script
#
# This script runs automatically when the VM boots.  It:
#   1. Installs Docker
#   2. Clones the MetaInformAnt repository
#   3. Builds the pipeline Docker image
#   4. Starts the pipeline with params from instance metadata
#   5. Optionally syncs results to GCS
#
# Logs: sudo journalctl -u google-startup-scripts.service -f
# ─────────────────────────────────────────────────────────────────────────
set -euo pipefail

LOGFILE="/var/log/metainformant-startup.log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "═══════════════════════════════════════════════════════════════════"
echo "  MetaInformAnt Pipeline — VM Startup"
echo "  $(date -Iseconds)"
echo "═══════════════════════════════════════════════════════════════════"

# ── Read instance metadata ──────────────────────────────────────────────
META_URL="http://metadata.google.internal/computeMetadata/v1/instance/attributes"
META_HEADER="Metadata-Flavor: Google"

get_meta() {
    curl -sf -H "$META_HEADER" "$META_URL/$1" 2>/dev/null || echo "$2"
}

REPO_URL=$(get_meta "pipeline-repo-url" "https://github.com/docxology/MetaInformAnt.git")
REPO_BRANCH=$(get_meta "pipeline-repo-branch" "main")
MAX_GB=$(get_meta "pipeline-max-gb" "20.0")
WORKERS=$(get_meta "pipeline-workers" "80")
THREADS=$(get_meta "pipeline-threads" "96")
GCS_BUCKET=$(get_meta "pipeline-gcs-bucket" "")
DOCKER_IMAGE=$(get_meta "pipeline-docker-image" "metainformant-pipeline")

echo "Config: max_gb=$MAX_GB workers=$WORKERS threads=$THREADS"
echo "Repo:   $REPO_URL @ $REPO_BRANCH"
echo "GCS:    ${GCS_BUCKET:-none}"

# ── Install Docker ──────────────────────────────────────────────────────
echo ""
echo "▸ Installing Docker..."
if ! command -v docker &>/dev/null; then
    apt-get update -qq
    apt-get install -y -qq docker.io
    systemctl enable --now docker
    echo "  ✓ Docker installed"
else
    echo "  ✓ Docker already installed"
fi

# ── Mount NVMe Local SSDs (If Attached) ───────────────────────────────
echo ""
echo "▸ Checking for Ephemeral NVMe Local SSDs..."
apt-get install -y -qq mdadm
NVME_DEVS=$(find /dev/disk/by-id/ -name 'google-local-nvme-ssd-*' | sort | xargs -r readlink -f | uniq || true)

WORK_DIR="/opt/MetaInformAnt"

if [ -n "$NVME_DEVS" ]; then
    DEVICES=($NVME_DEVS)
    NUM_DEVS=${#DEVICES[@]}
    echo "  Found $NUM_DEVS NVMe interfaces: ${DEVICES[*]}"
    
    # Check if md0 is already active to prevent re-formatting on ghost re-evaluations
    if ! grep -q md0 /proc/mdstat; then
        mdadm --create /dev/md0 --level=0 --raid-devices=$NUM_DEVS "${DEVICES[@]}" --force --run
        mkfs.ext4 -m 0 -F -q /dev/md0
    fi
    
    mkdir -p "$WORK_DIR"
    mount -o discard,defaults /dev/md0 "$WORK_DIR"
    chmod a+w "$WORK_DIR"
    echo "  ✓ Mounted NVMe RAID0 array to $WORK_DIR ($NUM_DEVS x 375GB)"
else
    echo "  ⚠ No NVMe devices found. Booting with Standard Persistent Disk."
    mkdir -p "$WORK_DIR"
fi

# ── Clone repository ────────────────────────────────────────────────────
echo ""
echo "▸ Cloning repository..."
if [ -d "$WORK_DIR/.git" ]; then
    cd "$WORK_DIR"
    git fetch origin "$REPO_BRANCH" && git checkout "$REPO_BRANCH" && git pull
    echo "  ✓ Repository updated"
else
    git clone --branch "$REPO_BRANCH" --depth 1 "$REPO_URL" "$WORK_DIR"
    echo "  ✓ Repository cloned"
fi
cd "$WORK_DIR"

# ── Build Docker image ──────────────────────────────────────────────────
# ── Build Docker image ──────────────────────────────────────────────────
echo ""
echo "▸ Building Docker image..."
docker build -t "$DOCKER_IMAGE" . 2>&1 | tail -5
echo "  ✓ Image built: $DOCKER_IMAGE"

# ── Create output directory ─────────────────────────────────────────────
mkdir -p "$WORK_DIR/projects"

# ── Start pipeline ──────────────────────────────────────────────────────
echo ""
echo "▸ Starting pipeline..."
docker run -d \
    --name metainformant-pipeline \
    --restart unless-stopped \
    -v "$WORK_DIR/projects:/app/projects" \
    -v "$WORK_DIR/config:/app/config" \
    -e "PIPELINE_MAX_GB=$MAX_GB" \
    -e "PIPELINE_WORKERS=$WORKERS" \
    -e "PIPELINE_THREADS=$THREADS" \
    --cpus="$THREADS" \
    --memory="90g" \
    "$DOCKER_IMAGE" \
    python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_apis_mellifera.yaml

echo "  ✓ Pipeline container started"
echo "  Container: $(docker ps --filter name=metainformant-pipeline --format '{{.ID}} {{.Status}}')"

# ── Set up GCS sync (if configured) ────────────────────────────────────
if [ -n "$GCS_BUCKET" ]; then
    echo ""
    echo "▸ Setting up periodic GCS sync..."

    # Sync every 30 minutes
    cat > /etc/cron.d/metainformant-sync << CRON
*/30 * * * * root gsutil -m rsync -r /opt/MetaInformAnt/projects/ gs://${GCS_BUCKET}/projects/ >> /var/log/metainformant-sync.log 2>&1
CRON
    echo "  ✓ GCS sync configured (every 30 min → gs://$GCS_BUCKET/amalgkit/)"
fi

# ── Done ────────────────────────────────────────────────────────────────
echo ""
echo "═══════════════════════════════════════════════════════════════════"
echo "  ✓ Pipeline deployment complete!"
echo "  Monitor: docker logs -f metainformant-pipeline"
echo "  Status:  docker ps"
echo "═══════════════════════════════════════════════════════════════════"
