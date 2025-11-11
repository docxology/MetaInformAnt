#!/bin/bash
# Monitor getfastq progress without blocking

LOG_DIR="output/amalgkit/pogonomyrmex_barbatus/logs"
FASTQ_DIR="output/amalgkit/pogonomyrmex_barbatus/fastq/getfastq"

echo "📊 getfastq Progress Monitor"
echo "Press Ctrl+C to stop"
echo ""

while true; do
    clear
    echo "═══════════════════════════════════════════════════════════════════════════════"
    echo "📊 getfastq Progress Monitor - $(date '+%Y-%m-%d %H:%M:%S')"
    echo "═══════════════════════════════════════════════════════════════════════════════"
    echo ""
    
    # Count samples with FASTQ files
    if [ -d "$FASTQ_DIR" ]; then
        total_samples=$(find "$FASTQ_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
        samples_with_fastq=0
        samples_with_sra=0
        
        for sample_dir in "$FASTQ_DIR"/*/; do
            if [ -d "$sample_dir" ]; then
                sample=$(basename "$sample_dir")
                if ls "$sample_dir"/*.fastq.gz 2>/dev/null | grep -q .; then
                    ((samples_with_fastq++))
                elif ls "$sample_dir"/*.sra 2>/dev/null | grep -q .; then
                    ((samples_with_sra++))
                fi
            fi
        done
        
        echo "📁 Sample Status:"
        echo "   Total sample directories: $total_samples"
        echo "   ✓ With FASTQ files: $samples_with_fastq"
        echo "   ⏳ With SRA (converting): $samples_with_sra"
        echo "   ⏸️  Pending: $((total_samples - samples_with_fastq - samples_with_sra))"
        echo ""
    fi
    
    # Show recent log activity
    echo "📝 Recent Activity (last 5 lines from most recent log):"
    most_recent=$(ls -t "$LOG_DIR"/*getfastq*.log 2>/dev/null | head -1)
    if [ -n "$most_recent" ]; then
        tail -5 "$most_recent" 2>/dev/null | sed 's/^/   /'
    else
        echo "   No log files found"
    fi
    echo ""
    
    # Show running processes
    echo "🔄 Running Processes:"
    ps aux | grep -E "amalgkit.*getfastq" | grep -v grep | head -3 | awk '{print "   PID", $2, "-", $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}' | sed 's/^/   /'
    echo ""
    
    echo "═══════════════════════════════════════════════════════════════════════════════"
    echo "Refreshing in 5 seconds... (Ctrl+C to stop)"
    sleep 5
done


