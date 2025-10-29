#!/bin/bash
# Monitor amalgkit workflow progress for multiple species

echo "================================================================================"
echo "🔍 AMALGKIT MULTI-SPECIES WORKFLOW MONITOR"
echo "================================================================================"
echo "Date: $(date)"
echo ""

# Check main log
MAIN_LOG="output/amalgkit_multi_species_run.log"
if [[ -f "$MAIN_LOG" ]]; then
    echo "=== Latest Log Entries (last 30 lines) ==="
    tail -30 "$MAIN_LOG"
    echo ""
fi

# Check species directories
echo "================================================================================"
echo "📊 SPECIES STATUS"
echo "================================================================================"

for SPECIES in sinvicta cfloridanus; do
    echo ""
    echo "--- $SPECIES ---"
    
    WORK_DIR="output/amalgkit/$SPECIES/work"
    LOG_DIR="output/amalgkit/$SPECIES/logs"
    
    if [[ -d "$WORK_DIR" ]]; then
        echo "Work directory: ✅ EXISTS"
        
        # Check for metadata
        if [[ -f "$WORK_DIR/metadata/metadata_original.tsv" ]]; then
            SAMPLE_COUNT=$(tail -n +2 "$WORK_DIR/metadata/metadata_original.tsv" 2>/dev/null | wc -l | tr -d ' ')
            echo "  Metadata: ✅ $SAMPLE_COUNT samples"
        fi
        
        # Check for quant outputs
        if [[ -d "$WORK_DIR/quant" ]]; then
            QUANT_COUNT=$(find "$WORK_DIR/quant" -maxdepth 1 -type d ! -name quant 2>/dev/null | wc -l | tr -d ' ')
            echo "  Quantifications: $QUANT_COUNT samples"
        fi
        
        # Check for merge outputs
        if [[ -d "$WORK_DIR/merge" ]]; then
            if [[ -f "$WORK_DIR/merge/metadata.tsv" ]]; then
                echo "  Merge: ✅ COMPLETE"
            else
                echo "  Merge: ⏳ IN PROGRESS"
            fi
        fi
        
        # Check for curate outputs
        if [[ -d "$WORK_DIR/curate" ]]; then
            CURATE_FILES=$(find "$WORK_DIR/curate" -type f 2>/dev/null | wc -l | tr -d ' ')
            if [[ $CURATE_FILES -gt 0 ]]; then
                echo "  Curate: ✅ $CURATE_FILES files"
            fi
        fi
        
    else
        echo "Work directory: ⏳ NOT YET CREATED"
    fi
    
    # Check logs
    if [[ -d "$LOG_DIR" ]]; then
        LATEST_LOG=$(ls -t "$LOG_DIR"/*.log 2>/dev/null | head -1)
        if [[ -n "$LATEST_LOG" ]]; then
            echo "  Latest log: $(basename "$LATEST_LOG")"
        fi
    fi
done

echo ""
echo "================================================================================"
echo "💡 USAGE"
echo "================================================================================"
echo "  Watch progress: watch -n 60 ./scripts/rna/monitor_amalgkit_progress.sh"
echo "  View full log:  tail -f output/amalgkit_multi_species_run.log"
echo "  Check species:  ls -lh output/amalgkit/{sinvicta,cfloridanus}/work/"
echo "================================================================================"
echo ""

