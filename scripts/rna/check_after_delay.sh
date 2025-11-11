#!/bin/bash
# Check workflow status after a delay

DELAY_HOURS=${1:-1}
DELAY_SECONDS=$((DELAY_HOURS * 3600))

echo "⏰ Waiting ${DELAY_HOURS} hour(s) before checking status..."
echo "Started at: $(date)"
echo "Will check at: $(date -d "+${DELAY_HOURS} hour")"
echo ""

sleep $DELAY_SECONDS

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📊 STATUS CHECK AFTER ${DELAY_HOURS} HOUR(S)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Check time: $(date)"
echo ""

cd /media/q/ext6/github/MetaInformAnt || exit 1

echo "=== WORKFLOW STATUS ==="
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status 2>&1 | tail -5
echo ""

echo "=== ACTIVE PROCESSES ==="
ACTIVE=$(ps aux | grep -E "amalgkit|run_workflow|process_samples" | grep -v grep | wc -l)
echo "Active processes: $ACTIVE"
echo ""

echo "=== DOWNLOAD PROGRESS ==="
FASTQ_COUNT=$(find output/amalgkit/pogonomyrmex_barbatus/fastq/getfastq -name "*.fastq*" -type f 2>/dev/null | wc -l)
SRA_COUNT=$(find output/amalgkit/pogonomyrmex_barbatus/fastq/getfastq -name "*.sra" -type f 2>/dev/null | wc -l)
QUANT_COUNT=$(find output/amalgkit/pogonomyrmex_barbatus/quant -name "abundance.tsv" 2>/dev/null | wc -l)
echo "FASTQ files: $FASTQ_COUNT"
echo "SRA files: $SRA_COUNT"
echo "Quantified samples: $QUANT_COUNT"
echo ""

echo "=== DISK USAGE ==="
du -sh output/amalgkit/pogonomyrmex_barbatus/fastq 2>/dev/null | awk '{print "Fastq directory: " $1}'
du -sh output/amalgkit/pogonomyrmex_barbatus/quant 2>/dev/null | awk '{print "Quant directory: " $1}'
du -sh output/amalgkit/pogonomyrmex_barbatus/fastq/temp/sra 2>/dev/null | awk '{print "SRA temp: " $1}'
echo ""

echo "=== RECENT ACTIVITY ==="
RECENT=$(find output/amalgkit/pogonomyrmex_barbatus/logs -name "*.log" -type f -mmin -10 2>/dev/null | wc -l)
echo "Log files updated in last 10 minutes: $RECENT"
echo ""

echo "=== DISK SPACE ==="
df -h /media/q/ext6 | tail -1 | awk '{print "External drive: " $4 " free (" $5 " used)"}'
df -h /tmp | tail -1 | awk '{print "/tmp: " $4 " free (" $5 " used)"}'
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅ Status check complete"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

