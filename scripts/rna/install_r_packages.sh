#!/bin/bash
# Install R packages with proper temp directory configuration
# Usage: ./scripts/rna/install_r_packages.sh [package1] [package2] ...

set -e

# Get repository root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

# Set temp directory to repository .tmp/R
export TMPDIR="$REPO_ROOT/.tmp/R"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
mkdir -p "$TMPDIR"

echo "Using temp directory: $TMPDIR"
echo "Available space: $(df -h "$TMPDIR" | tail -1 | awk '{print $4}')"

# Install packages
if [ $# -eq 0 ]; then
    echo "Installing ggplot2..."
    Rscript -e "install.packages('ggplot2', repos='https://cloud.r-project.org', lib='~/R/library')"
else
    for package in "$@"; do
        echo "Installing $package..."
        Rscript -e "install.packages('$package', repos='https://cloud.r-project.org', lib='~/R/library')"
    done
fi

echo "✅ Package installation complete"
