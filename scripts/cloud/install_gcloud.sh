#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────
# Install Google Cloud CLI on Linux
#
# Usage:  bash scripts/cloud/install_gcloud.sh
# Docs:   https://cloud.google.com/sdk/docs/install
# ─────────────────────────────────────────────────────────────────────────
set -euo pipefail

echo "▸ Installing Google Cloud CLI..."

# Check if already installed
if command -v gcloud &>/dev/null; then
    echo "  ✓ gcloud already installed: $(gcloud --version 2>/dev/null | head -1)"
    exit 0
fi

# Install via apt (Debian/Ubuntu/Parrot)
if command -v apt-get &>/dev/null; then
    echo "  Using apt repository..."

    # Add Google Cloud GPG key
    curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
        | sudo gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg 2>/dev/null

    # Add repository
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
        | sudo tee /etc/apt/sources.list.d/google-cloud-sdk.list > /dev/null

    sudo apt-get update -qq
    sudo apt-get install -y -qq google-cloud-cli

# Fallback: standalone installer
else
    echo "  Using standalone installer..."
    cd /tmp
    curl -fsSLO https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-linux-x86_64.tar.gz
    tar xzf google-cloud-cli-linux-x86_64.tar.gz
    ./google-cloud-sdk/install.sh --quiet --path-update true
    rm -f google-cloud-cli-linux-x86_64.tar.gz
    export PATH="$PATH:/tmp/google-cloud-sdk/bin"
fi

echo ""
echo "  ✓ gcloud installed: $(gcloud --version 2>/dev/null | head -1)"
echo ""
echo "▸ Next steps:"
echo "  1. Authenticate:  gcloud auth login"
echo "  2. Set project:   gcloud config set project YOUR_PROJECT_ID"
echo "  3. Deploy:        python scripts/cloud/deploy_gcp.py deploy --project YOUR_PROJECT_ID"
