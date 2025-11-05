#!/bin/bash
# Genome Setup Scripts - Run Sequentially
# This script runs all genome setup steps in order

set -e  # Exit on error

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

echo "================================================================================"
echo "GENOME SETUP PIPELINE"
echo "================================================================================"
echo "Repository root: $REPO_ROOT"
echo ""

# Step 1: Verify current status
echo "================================================================================"
echo "Step 1: Verifying genome downloads and indexes"
echo "================================================================================"
python3 scripts/rna/verify_genomes_and_indexes.py
echo ""

# Step 2: Download missing genomes
echo "================================================================================"
echo "Step 2: Downloading missing genomes"
echo "================================================================================"
python3 scripts/rna/download_missing_genomes.py
echo ""

# Step 3: Prepare transcriptomes
echo "================================================================================"
echo "Step 3: Preparing transcriptomes"
echo "================================================================================"
python3 scripts/rna/prepare_transcriptomes.py
echo ""

# Step 4: Build kallisto indexes
echo "================================================================================"
echo "Step 4: Building kallisto indexes"
echo "================================================================================"
python3 scripts/rna/build_kallisto_indexes.py
echo ""

# Step 5: Final verification
echo "================================================================================"
echo "Step 5: Final verification"
echo "================================================================================"
python3 scripts/rna/verify_genomes_and_indexes.py
echo ""

echo "================================================================================"
echo "GENOME SETUP PIPELINE COMPLETE"
echo "================================================================================"

