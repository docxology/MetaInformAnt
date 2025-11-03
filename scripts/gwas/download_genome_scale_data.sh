#!/bin/bash
# Download real Apis mellifera samples for genome-scale GWAS

set -e

WORK_DIR="/home/q/Documents/GitHub/MetaInformAnt"
SRA_DIR="$WORK_DIR/data/raw/sra"
THREADS=8

# Target samples (Scout and Recruit behavioral phenotypes)
SAMPLES=("SRR2096937" "SRR2096938" "SRR2096939")

echo "╔════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                            ║"
echo "║           DOWNLOADING REAL APIS MELLIFERA GENOME-SCALE DATA               ║"
echo "║                                                                            ║"
echo "╚════════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "BioProject: PRJNA292680 (Scout/Recruit Behavioral Study)"
echo "Samples: ${SAMPLES[@]}"
echo "Expected size: ~10-15 GB"
echo "Estimated time: 1-2 hours"
echo ""

mkdir -p "$SRA_DIR"
cd "$SRA_DIR"

for sample in "${SAMPLES[@]}"; do
    echo "════════════════════════════════════════════════════════════════"
    echo "Downloading: $sample"
    echo "════════════════════════════════════════════════════════════════"
    
    # Check if already downloaded
    if [ -f "${sample}_1.fastq" ] && [ -f "${sample}_2.fastq" ]; then
        echo "✓ $sample already downloaded, skipping"
        continue
    fi
    
    # Download with fasterq-dump
    echo "Starting download..."
    fasterq-dump $sample -e $THREADS --split-files -p || {
        echo "✗ Failed to download $sample"
        continue
    }
    
    # Check result
    if [ -f "${sample}_1.fastq" ] && [ -f "${sample}_2.fastq" ]; then
        SIZE1=$(du -h "${sample}_1.fastq" | cut -f1)
        SIZE2=$(du -h "${sample}_2.fastq" | cut -f1)
        echo "✓ Downloaded successfully:"
        echo "  Read 1: $SIZE1"
        echo "  Read 2: $SIZE2"
    else
        echo "✗ Download incomplete"
    fi
    echo ""
done

echo "════════════════════════════════════════════════════════════════"
echo "Download Summary"
echo "════════════════════════════════════════════════════════════════"
ls -lh "$SRA_DIR"/*.fastq 2>/dev/null || echo "No FASTQ files found"
echo ""
echo "Total size:"
du -sh "$SRA_DIR"

echo ""
echo "✓ Download complete! Next steps:"
echo "  1. Install BWA, SAMtools, bcftools (see INSTALL.md)"
echo "  2. Run alignment: bash scripts/align_samples.sh"
echo "  3. Call variants: bash scripts/call_variants.sh"

