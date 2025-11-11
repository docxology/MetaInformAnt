#!/bin/bash
# Monitor getfastq progress with periodic updates

SPECIES="pogonomyrmex_barbatus"
FASTQ_DIR="output/amalgkit/${SPECIES}/fastq/getfastq"
UPDATE_INTERVAL=30  # seconds

echo "📊 getfastq Progress Monitor"
echo "Press Ctrl+C to stop"
echo ""

while true; do
    clear
    echo "═══════════════════════════════════════════════════════════════════════════════"
    echo "📊 getfastq Progress Monitor - $(date '+%Y-%m-%d %H:%M:%S')"
    echo "═══════════════════════════════════════════════════════════════════════════════"
    echo ""
    
    # Run status check
    bash scripts/rna/check_getfastq_status.sh
    
    echo ""
    echo "⏱️  Next update in ${UPDATE_INTERVAL} seconds (Ctrl+C to stop)"
    
    sleep "$UPDATE_INTERVAL"
done

