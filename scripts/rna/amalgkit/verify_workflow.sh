#!/bin/bash
# Centralized amalgkit workflow verification script
# Works for ANY species directory
# 
# USAGE:
#   scripts/rna/amalgkit/verify_workflow.sh <species_name>
#   scripts/rna/amalgkit/verify_workflow.sh pbarbatus
#   scripts/rna/amalgkit/verify_workflow.sh sinvicta
#
# Or run from species directory:
#   cd output/amalgkit/pbarbatus
#   ../../../scripts/rna/amalgkit/verify_workflow.sh .

set -e  # Exit on error

# Determine species directory
if [[ -z "$1" ]]; then
    # No argument - check if we're in a species directory
    if [[ -d "work/quant" ]]; then
        SPECIES_DIR="$(pwd)"
    else
        echo "❌ ERROR: No species specified and not in a species directory"
        echo ""
        echo "USAGE:"
        echo "  scripts/rna/amalgkit/verify_workflow.sh <species_name>"
        echo "  scripts/rna/amalgkit/verify_workflow.sh pbarbatus"
        echo ""
        echo "Or run from species directory:"
        echo "  cd output/amalgkit/pbarbatus"
        echo "  ../../../scripts/rna/amalgkit/verify_workflow.sh ."
        exit 1
    fi
elif [[ "$1" == "." ]]; then
    # Current directory
    SPECIES_DIR="$(pwd)"
elif [[ -d "$1" ]]; then
    # Full path provided
    SPECIES_DIR="$1"
elif [[ -d "output/amalgkit/$1" ]]; then
    # Species name provided (from repo root)
    SPECIES_DIR="output/amalgkit/$1"
else
    echo "❌ ERROR: Species directory not found: $1"
    exit 1
fi

# Convert to absolute path
SPECIES_DIR="$(cd "$SPECIES_DIR" && pwd)"

echo "================================================================================"
echo "🔍 AMALGKIT WORKFLOW VERIFICATION"
echo "================================================================================"
echo "Species directory: $SPECIES_DIR"
echo "Date: $(date)"
echo ""

# Check directory structure
if [[ ! -d "$SPECIES_DIR/work" ]]; then
    echo "❌ ERROR: Not a valid amalgkit species directory!"
    echo "Expected structure:"
    echo "  $SPECIES_DIR/"
    echo "  └── work/"
    echo "      ├── metadata/"
    echo "      ├── genome/"
    echo "      ├── quant/"
    echo "      ├── merge/"
    echo "      ├── curate/"
    echo "      └── sanity/"
    exit 1
fi

# Get species name
SPECIES_NAME=$(basename "$SPECIES_DIR")
echo "✅ Species: $SPECIES_NAME"
echo "✅ Valid amalgkit directory structure"
echo ""

# Change to species directory for amalgkit commands
cd "$SPECIES_DIR"

# ============================================================================
# SANITY CHECK
# ============================================================================

echo "================================================================================"
echo "📊 STEP 1: SANITY CHECK"
echo "================================================================================"

if [[ ! -d "work/quant" ]]; then
    echo "⚠️  No quantification directory found - workflow may not be complete"
    echo ""
else
    echo "Running: amalgkit sanity --out_dir work --all"
    echo ""
    
    amalgkit sanity --out_dir work --all 2>&1 | tail -20
    
    SANITY_EXIT=$?
    
    echo ""
    echo "Sanity exit code: $SANITY_EXIT"
    
    if [[ $SANITY_EXIT -eq 0 ]]; then
        echo "✅ Sanity check PASSED"
    else
        echo "❌ Sanity check FAILED"
        exit 1
    fi
    
    # Check outputs
    echo ""
    echo "Sanity outputs:"
    if [[ -f "work/sanity/SRA_IDs_without_fastq.txt" ]]; then
        FASTQ_COUNT=$(wc -l < work/sanity/SRA_IDs_without_fastq.txt)
        echo "  ✅ SRA_IDs_without_fastq.txt: $FASTQ_COUNT samples (expected if FASTQs deleted)"
    fi
    
    if [[ -f "work/sanity/SRA_IDs_without_quant.txt" ]]; then
        QUANT_COUNT=$(wc -l < work/sanity/SRA_IDs_without_quant.txt)
        if [[ $QUANT_COUNT -gt 0 ]]; then
            echo "  ⚠️  SRA_IDs_without_quant.txt: $QUANT_COUNT samples missing quantification!"
            echo "      This indicates incomplete workflow"
        fi
    else
        echo "  ✅ No SRA_IDs_without_quant.txt (all samples validated!)"
    fi
fi

echo ""

# ============================================================================
# CURATE CHECK
# ============================================================================

echo "================================================================================"
echo "📊 STEP 2: CURATE CHECK"
echo "================================================================================"

echo "Checking curate outputs..."
echo ""

CURATE_FILES=$(find work/curate -type f 2>/dev/null | wc -l)
CURATE_PDFS=$(find work/curate -name "*.pdf" 2>/dev/null | wc -l)
CURATE_TSVS=$(find work/curate -name "*.tsv" 2>/dev/null | wc -l)
CURATE_RDATA=$(find work/curate -name "*.RData" 2>/dev/null | wc -l)

echo "Current curate outputs:"
echo "  Total files: $CURATE_FILES"
echo "  PDF visualizations: $CURATE_PDFS"
echo "  TSV data tables: $CURATE_TSVS"
echo "  RData files: $CURATE_RDATA"
echo ""

if [[ $CURATE_FILES -ge 17 ]] && [[ $CURATE_PDFS -ge 6 ]]; then
    echo "✅ Curate outputs complete"
elif [[ -d "work/merge" ]]; then
    echo "⚠️  Curate outputs incomplete - regenerating..."
    echo ""
    echo "Running: amalgkit curate --out_dir work --batch_effect_alg no"
    echo ""
    
    amalgkit curate --out_dir work --batch_effect_alg no 2>&1 | tail -30
    
    CURATE_EXIT=$?
    echo ""
    echo "Curate exit code: $CURATE_EXIT"
    
    if [[ $CURATE_EXIT -eq 0 ]]; then
        echo "✅ Curate step PASSED"
        # Recount files
        CURATE_FILES=$(find work/curate -type f 2>/dev/null | wc -l)
        CURATE_PDFS=$(find work/curate -name "*.pdf" 2>/dev/null | wc -l)
        CURATE_TSVS=$(find work/curate -name "*.tsv" 2>/dev/null | wc -l)
    else
        echo "❌ Curate step FAILED"
    fi
else
    echo "⚠️  Merge step not yet complete - curate cannot run"
fi

echo ""

# ============================================================================
# FINAL SUMMARY
# ============================================================================

echo "================================================================================"
echo "✅ VERIFICATION SUMMARY"
echo "================================================================================"

# Sample count
SAMPLE_COUNT=0
if [[ -d "work/quant" ]]; then
    SAMPLE_COUNT=$(find work/quant -maxdepth 1 -type d ! -name quant 2>/dev/null | wc -l | tr -d ' ')
fi

echo ""
echo "Species: $SPECIES_NAME"
echo "Location: $SPECIES_DIR"
echo ""
echo "Workflow Status:"
if [[ -f "work/metadata/metadata_original.tsv" ]]; then
    METADATA_COUNT=$(tail -n +2 work/metadata/metadata_original.tsv 2>/dev/null | wc -l | tr -d ' ')
    echo "  ✅ Metadata: $METADATA_COUNT samples"
else
    echo "  ⏳ Metadata: Not yet retrieved"
fi

if [[ -d "work/genome" ]]; then
    GENOME_FILES=$(find work/genome -type f 2>/dev/null | wc -l)
    echo "  ✅ Genome: $GENOME_FILES files downloaded"
else
    echo "  ⏳ Genome: Not yet downloaded"
fi

if [[ $SAMPLE_COUNT -gt 0 ]]; then
    echo "  ✅ Quantification: $SAMPLE_COUNT samples"
else
    echo "  ⏳ Quantification: Not yet started"
fi

if [[ -d "work/merge" ]]; then
    echo "  ✅ Merge: Complete"
else
    echo "  ⏳ Merge: Not yet complete"
fi

if [[ $CURATE_PDFS -ge 6 ]]; then
    echo "  ✅ Curate: $CURATE_PDFS PDFs, $CURATE_TSVS tables"
else
    echo "  ⏳ Curate: Not yet complete"
fi

echo ""
echo "Quick Access:"
if [[ $CURATE_PDFS -gt 0 ]]; then
    echo "  View plots: open $SPECIES_DIR/work/curate/*/plots/*.pdf"
fi
if [[ $CURATE_TSVS -gt 0 ]]; then
    echo "  Expression matrix: $SPECIES_DIR/work/curate/*/tables/*.no.tc.tsv"
fi

echo ""
echo "================================================================================"

exit 0

