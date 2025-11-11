#!/bin/bash
# Comprehensive Workflow Monitoring Commands
# Pogonomyrmex barbatus Workflow

echo "=== Quick Status Check ==="
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status --detailed

echo ""
echo "=== Process Count ==="
echo "Active processes: $(ps aux | grep -E 'run_workflow|amalgkit' | grep -v grep | wc -l)"

echo ""
echo "=== Quantification Progress ==="
echo "Quantified: $(find output/amalgkit/pogonomyrmex_barbatus/quant -name 'abundance.tsv' 2>/dev/null | wc -l) / 83"

echo ""
echo "=== Output Verification ==="
echo "Genome index: $(test -f output/amalgkit/pogonomyrmex_barbatus/work/index/Pogonomyrmex_barbatus_transcripts.idx && echo '✓' || echo '✗')"
echo "Merged matrix: $(test -f output/amalgkit/pogonomyrmex_barbatus/merged/merged_abundance.tsv && echo '✓' || echo '⏳')"
echo "Curate results: $(test -d output/amalgkit/pogonomyrmex_barbatus/curate && echo '✓' || echo '⏳')"
echo "CSCA results: $(test -d output/amalgkit/pogonomyrmex_barbatus/csca && echo '✓' || echo '⏳')"

echo ""
echo "=== Manifest Check ==="
if [ -f output/amalgkit/pogonomyrmex_barbatus/work/amalgkit.manifest.jsonl ]; then
    echo "Steps completed: $(wc -l < output/amalgkit/pogonomyrmex_barbatus/work/amalgkit.manifest.jsonl)"
    echo "Last 5 steps:"
    tail -5 output/amalgkit/pogonomyrmex_barbatus/work/amalgkit.manifest.jsonl | python3 -c "import json, sys; [print(f\"  {json.loads(l)['step']}: return_code={json.loads(l).get('return_code', 'N/A')}\") for l in sys.stdin]"
else
    echo "Manifest not yet created"
fi
