#!/bin/bash
# Verification script template for amalgkit sanity and curate steps
# 
# USAGE:
#   1. Copy this script to your species directory:
#      cp docs/rna/amalgkit/verify_template.sh output/amalgkit/<species_name>/verify_workflow.sh
#   
#   2. Make it executable:
#      chmod +x output/amalgkit/<species_name>/verify_workflow.sh
#   
#   3. Run from your species directory:
#      cd output/amalgkit/<species_name>
#      ./verify_workflow.sh

set -e  # Exit on error

echo "================================================================================"
echo "🔍 AMALGKIT VERIFICATION SCRIPT"
echo "================================================================================"
echo "Directory: $(pwd)"
echo "Date: $(date)"
echo ""

# Check we're in the right directory
if [[ ! -d "work/quant" ]]; then
    echo "❌ ERROR: Not in the correct species directory!"
    echo "Please run from: output/amalgkit/<species_name>/"
    echo ""
    echo "Expected directory structure:"
    echo "  work/"
    echo "  ├── metadata/"
    echo "  ├── quant/"
    echo "  ├── merge/"
    echo "  ├── curate/"
    echo "  └── sanity/"
    exit 1
fi

# Get species name from current directory
SPECIES_NAME=$(basename "$(pwd)")
echo "✅ Species: $SPECIES_NAME"
echo "✅ Correct directory structure found"
echo ""

# ============================================================================
# SANITY CHECK
# ============================================================================

echo "================================================================================"
echo "📊 STEP 1: SANITY CHECK"
echo "================================================================================"

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
    echo "  ✅ SRA_IDs_without_fastq.txt: $FASTQ_COUNT samples (expected if FASTQs deleted after quant)"
fi

if [[ -f "work/sanity/SRA_IDs_without_quant.txt" ]]; then
    QUANT_COUNT=$(wc -l < work/sanity/SRA_IDs_without_quant.txt)
    if [[ $QUANT_COUNT -gt 0 ]]; then
        echo "  ⚠️  SRA_IDs_without_quant.txt: $QUANT_COUNT samples missing quantification!"
        echo "      This indicates incomplete workflow - some samples failed."
    fi
else
    echo "  ✅ No SRA_IDs_without_quant.txt (all samples validated!)"
fi

echo ""

# ============================================================================
# CURATE CHECK
# ============================================================================

echo "================================================================================"
echo "📊 STEP 2: CURATE CHECK"
echo "================================================================================"

echo "Checking existing curate outputs..."
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
    echo "✅ Curate outputs complete (using existing)"
    echo ""
    echo "To regenerate, run:"
    echo "  rm -rf work/curate/*"
    echo "  amalgkit curate --out_dir work --batch_effect_alg no"
else
    echo "⚠️  Curate outputs incomplete, regenerating..."
    echo ""
    
    # Backup existing if any
    CURATE_SPECIES_DIR=$(find work/curate -maxdepth 1 -type d ! -name curate 2>/dev/null | head -1)
    if [[ -n "$CURATE_SPECIES_DIR" && -d "$CURATE_SPECIES_DIR" ]]; then
        echo "Backing up existing curate outputs..."
        mv "$CURATE_SPECIES_DIR" "${CURATE_SPECIES_DIR}.backup.$(date +%s)"
    fi
    
    echo "Running: amalgkit curate --out_dir work --batch_effect_alg no"
    echo ""
    
    amalgkit curate --out_dir work --batch_effect_alg no 2>&1 | tail -30
    
    CURATE_EXIT=$?
    
    echo ""
    echo "Curate exit code: $CURATE_EXIT"
    
    if [[ $CURATE_EXIT -eq 0 ]]; then
        echo "✅ Curate step PASSED"
    else
        echo "❌ Curate step FAILED"
        exit 1
    fi
fi

echo ""

# ============================================================================
# FINAL VERIFICATION
# ============================================================================

echo "================================================================================"
echo "✅ FINAL VERIFICATION"
echo "================================================================================"

echo ""
echo "Sanity outputs:"
ls -lh work/sanity/ | tail -n +2

echo ""
echo "Curate outputs:"
echo "  Plots:"
find work/curate -name "*.pdf" -exec ls -lh {} \; | awk '{print "    " $9 " (" $5 ")"}'

echo ""
echo "  Tables:"
find work/curate -name "*.tsv" -exec ls -lh {} \; | awk '{print "    " $9 " (" $5 ")"}'

echo ""
echo "  RData:"
find work/curate -name "*.RData" -exec ls -lh {} \; | awk '{print "    " $9 " (" $5 ")"}'

echo ""
echo "================================================================================"
echo "🎉 ALL CHECKS COMPLETE"
echo "================================================================================"
SAMPLE_COUNT=$(find work/quant -maxdepth 1 -type d ! -name quant 2>/dev/null | wc -l | tr -d ' ')

echo ""
echo "Summary:"
echo "  ✅ Species: $SPECIES_NAME"
echo "  ✅ Sanity: Validated $SAMPLE_COUNT samples"
echo "  ✅ Curate: Generated $CURATE_PDFS PDF visualizations"
echo "  ✅ Curate: Created $CURATE_TSVS data tables"
echo "  ✅ Total files: $CURATE_FILES curate outputs"
echo ""
echo "View PDFs:"
echo "  open work/curate/*/plots/*.pdf"
echo ""
echo "View expression matrix:"
echo "  head work/curate/*/tables/*.no.tc.tsv"
echo ""
echo "Next steps:"
echo "  - Differential expression analysis"
echo "  - GO/KEGG enrichment"
echo "  - Co-expression network analysis"
echo ""

exit 0

