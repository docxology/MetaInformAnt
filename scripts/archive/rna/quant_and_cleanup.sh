#!/bin/bash
# Quantify completed samples and delete their FASTQ files to free disk space

set -e

CONFIG_FILE="config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml"
SPECIES="pogonomyrmex_barbatus"
FASTQ_DIR="output/amalgkit/${SPECIES}/fastq/getfastq"
QUANT_DIR="output/amalgkit/${SPECIES}/quant"

echo "═══════════════════════════════════════════════════════════════════════════════"
echo "🔬 Quantify and Cleanup Completed Samples"
echo "═══════════════════════════════════════════════════════════════════════════════"
echo ""

# Check if config exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ Config file not found: $CONFIG_FILE"
    exit 1
fi

# Find samples with FASTQ files
echo "📋 Finding samples with completed FASTQ files..."
COMPLETED_SAMPLES=()
for sample_dir in "$FASTQ_DIR"/SRR*/; do
    if [ -d "$sample_dir" ]; then
        sample=$(basename "$sample_dir")
        has_fastq=$(find "$sample_dir" -name "*.fastq*" ! -name "*.sra" 2>/dev/null | wc -l)
        if [ "$has_fastq" -gt 0 ]; then
            COMPLETED_SAMPLES+=("$sample")
        fi
    fi
done

if [ ${#COMPLETED_SAMPLES[@]} -eq 0 ]; then
    echo "⚠️  No samples with completed FASTQ files found"
    echo "   (Samples may still be downloading or converting)"
    exit 0
fi

echo "   Found ${#COMPLETED_SAMPLES[@]} samples with FASTQ files:"
for sample in "${COMPLETED_SAMPLES[@]}"; do
    echo "   - $sample"
done
echo ""

# Check which are already quantified
echo "🔍 Checking which samples are already quantified..."
NEEDS_QUANT=()
for sample in "${COMPLETED_SAMPLES[@]}"; do
    # Check if abundance file exists
    abundance_file="$QUANT_DIR/${sample}/abundance.tsv"
    if [ ! -f "$abundance_file" ]; then
        NEEDS_QUANT+=("$sample")
    else
        echo "   ✓ $sample already quantified"
    fi
done

if [ ${#NEEDS_QUANT[@]} -eq 0 ]; then
    echo ""
    echo "✅ All completed samples are already quantified"
    echo ""
    echo "🗑️  Proceeding to delete FASTQ files..."
    for sample in "${COMPLETED_SAMPLES[@]}"; do
        sample_dir="$FASTQ_DIR/$sample"
        if [ -d "$sample_dir" ]; then
            echo "   Deleting $sample..."
            rm -rf "$sample_dir"
            echo "   ✓ Deleted $sample"
        fi
    done
    echo ""
    echo "✅ Cleanup complete!"
    exit 0
fi

echo ""
echo "📊 Samples needing quantification: ${#NEEDS_QUANT[@]}"
echo ""

# Option 1: Use Python cleanup function (recommended)
echo "═══════════════════════════════════════════════════════════════════════════════"
echo "Option 1: Use built-in cleanup function (RECOMMENDED)"
echo "═══════════════════════════════════════════════════════════════════════════════"
echo ""
echo "This will:"
echo "  1. Quantify all downloaded but unquantified samples"
echo "  2. Delete FASTQ files after successful quantification"
echo ""
echo "Run:"
echo "  python3 -c \"
from metainformant.rna.orchestration import cleanup_unquantified_samples
from pathlib import Path

config_path = Path('$CONFIG_FILE')
quantified, failed = cleanup_unquantified_samples(config_path)
print(f'✅ Quantified: {quantified}')
print(f'❌ Failed: {failed}')
\""
echo ""

# Option 2: Run workflow with quant step
echo "═══════════════════════════════════════════════════════════════════════════════"
echo "Option 2: Run workflow with quant step (automatic sequential mode)"
echo "═══════════════════════════════════════════════════════════════════════════════"
echo ""
echo "This will automatically:"
echo "  - Process remaining downloads"
echo "  - Quantify completed samples"
echo "  - Delete FASTQ files after quantification"
echo ""
echo "Run:"
echo "  python3 scripts/rna/run_workflow.py --config $CONFIG_FILE --steps getfastq quant"
echo ""

# Option 3: Manual per-sample processing
echo "═══════════════════════════════════════════════════════════════════════════════"
echo "Option 3: Manual per-sample processing"
echo "═══════════════════════════════════════════════════════════════════════════════"
echo ""
echo "For each sample, run:"
echo "  1. Quantify: amalgkit quant --id <SRR_ID> --out_dir <quant_dir>"
echo "  2. Delete: rm -rf $FASTQ_DIR/<SRR_ID>"
echo ""

