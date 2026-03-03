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

# ── Clone repository ────────────────────────────────────────────────────
echo ""
echo "▸ Cloning repository..."
WORK_DIR="/opt/MetaInformAnt"
if [ -d "$WORK_DIR" ]; then
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
mkdir -p "$WORK_DIR/output/amalgkit"

# ── Start pipeline ──────────────────────────────────────────────────────
echo ""
echo "▸ Starting pipeline..."
docker run -d \
    --name metainformant-pipeline \
    --restart unless-stopped \
    -v "$WORK_DIR/output:/app/output" \
    -v "$WORK_DIR/config:/app/config" \
    -e "PIPELINE_MAX_GB=$MAX_GB" \
    -e "PIPELINE_WORKERS=$WORKERS" \
    -e "PIPELINE_THREADS=$THREADS" \
    --cpus="$THREADS" \
    --memory="90g" \
    "$DOCKER_IMAGE"

echo "  ✓ Pipeline container started"
echo "  Container: $(docker ps --filter name=metainformant-pipeline --format '{{.ID}} {{.Status}}')"

# ── Set up GCS sync (if configured) ────────────────────────────────────
if [ -n "$GCS_BUCKET" ]; then
    echo ""
    echo "▸ Setting up periodic GCS sync..."

    # Sync every 30 minutes
    cat > /etc/cron.d/metainformant-sync << CRON
*/30 * * * * root gsutil -m rsync -r /opt/MetaInformAnt/output/amalgkit/ gs://${GCS_BUCKET}/amalgkit/ >> /var/log/metainformant-sync.log 2>&1
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
