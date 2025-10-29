#!/bin/bash
# Restart command for ENA batch processing
# Location: /Users/4d/Documents/GitHub/metainformant/output/amalgkit/pbarbatus/

cd /Users/4d/Documents/GitHub/metainformant

echo "🚀 Starting ENA batch processor..."
echo ""
echo "Configuration:"
echo "  • Concurrent downloads: 5"
echo "  • Concurrent quantifications: 3"
echo "  • Kallisto threads: 3 per sample"
echo ""

# Start in background
nohup python3 output/amalgkit/pbarbatus/batch_ena.py \
  > output/amalgkit/pbarbatus/batch_ena.log 2>&1 &

# Save PID
NEW_PID=$!
echo $NEW_PID > output/amalgkit/pbarbatus/batch_ena_pid.txt

echo "✅ Started (PID: $NEW_PID)"
echo ""
echo "Monitor with:"
echo "  # Check progress"
echo "  ls output/amalgkit/pbarbatus/work/quant/*/abundance.tsv | wc -l"
echo ""
echo "  # View active processes"
echo "  ps aux | grep -E 'wget|kallisto' | grep -v grep"
echo ""
echo "  # Check storage"
echo "  du -sh output/amalgkit/pbarbatus/work/fastq"

