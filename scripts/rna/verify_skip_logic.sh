#!/bin/bash
# Comprehensive verification script for skip logic fix

set -e

cd "$(dirname "$0")/../.."

echo "================================================================================="
echo "COMPREHENSIVE SKIP LOGIC VERIFICATION"
echo "================================================================================="
echo ""

echo "Step 1: Cleanup..."
pkill -9 -f "run_multi_species\|amalgkit\|fastq-dump\|fastp" 2>/dev/null || true
sleep 2
find output/amalgkit -name "metadata.batch*.tsv" -delete 2>/dev/null || true
echo "✅ Cleanup complete"
echo ""

echo "Step 2: Verify Pbarbatus quantified samples..."
pbarbatus_quant=$(find output/amalgkit/pbarbatus/quant -name "abundance.tsv" 2>/dev/null | wc -l | tr -d ' ')
echo "   Found: $pbarbatus_quant/83 quantified samples"
if [ "$pbarbatus_quant" -eq 83 ]; then
  echo "   ✅ Perfect - ready for skip test"
else
  echo "   ⚠️  Expected 83, found $pbarbatus_quant"
  exit 1
fi
echo ""

echo "Step 3: Test skip function directly..."
python3 << 'PYTHON_TEST'
import sys
from pathlib import Path
sys.path.insert(0, 'src')
from metainformant.rna.steps.batched_process import _sample_already_quantified

quant_dir = Path('output/amalgkit/pbarbatus/quant').absolute()
print(f"Testing skip function:")
print(f"  Quant dir: {quant_dir}")

# Test a few samples
test_samples = list(quant_dir.glob('SRR*/abundance.tsv'))[:5]
for sample_file in test_samples:
    run_id = sample_file.parent.name
    result = _sample_already_quantified(run_id, quant_dir)
    status = "✅" if result else "❌"
    print(f"  {status} {run_id}: {'quantified' if result else 'NOT quantified'}")

# Count all
all_samples = [d.name for d in quant_dir.glob('SRR*') if d.is_dir()]
quantified = sum(1 for s in all_samples if _sample_already_quantified(s, quant_dir))
print(f"\n  Total: {quantified}/{len(all_samples)} quantified")
if quantified == len(all_samples) == 83:
    print("  ✅ All samples correctly identified as quantified")
else:
    print("  ⚠️  Mismatch")
    sys.exit(1)
PYTHON_TEST

echo ""
echo "Step 4: Starting workflow..."
LOG_FILE="output/skip_logic_verification_$(date +%Y%m%d_%H%M%S).log"
python3 scripts/rna/run_multi_species_amalgkit.py > "$LOG_FILE" 2>&1 &
WORKFLOW_PID=$!
echo "   Workflow PID: $WORKFLOW_PID"
echo "   Log file: $LOG_FILE"
echo ""

echo "Step 5: Monitoring workflow (waiting 6 minutes for Pbarbatus)..."
for i in {1..12}; do
  sleep 30
  if grep -q "Starting workflow for.*pbarbatus" "$LOG_FILE" 2>/dev/null; then
    echo "   ✅ Pbarbatus reached at $(date)"
    break
  fi
  echo "   ... waiting ($i/12) - $(date +%H:%M:%S)"
done
echo ""

echo "Step 6: Analyzing results..."
echo ""

if grep -q "Starting workflow for.*pbarbatus" "$LOG_FILE" 2>/dev/null; then
  echo "=== PBARBATUS PROCESSING SECTION ==="
  awk '/Starting workflow for.*pbarbatus/,/Starting workflow for|CROSS-SPECIES|FINAL SUMMARY/' "$LOG_FILE" 2>/dev/null | head -80
  
  echo ""
  echo "=== SKIP LOGIC MESSAGES ==="
  if grep -q "🔍 Checking for already-quantified" "$LOG_FILE" 2>/dev/null; then
    echo "✅ Found: 'Checking for already-quantified samples' message"
  else
    echo "⚠️  Missing: 'Checking for already-quantified samples' message"
  fi
  
  if grep -q "⏭️.*Skipping.*83\|Skipping 83.*already-quantified" "$LOG_FILE" 2>/dev/null; then
    echo "✅ Found: Skip message for 83 samples"
    grep "⏭️.*Skipping.*83\|Skipping 83.*already-quantified" "$LOG_FILE" 2>/dev/null | head -3
  else
    echo "⚠️  Missing: Skip message for 83 samples"
  fi
  
  if grep -q "All samples already quantified" "$LOG_FILE" 2>/dev/null; then
    echo "✅ Found: 'All samples already quantified' message"
    grep "All samples already quantified" "$LOG_FILE" 2>/dev/null | head -3
  else
    echo "⚠️  Missing: 'All samples already quantified' message"
  fi
  
  echo ""
  echo "=== DOWNLOAD PROCESS CHECK ==="
  PBARBATUS_PROCS=$(ps aux | grep -c "[a]malgkit.*pbarbatus.*getfastq\|[f]astq-dump.*pbarbatus\|[f]astp.*pbarbatus" || echo "0")
  if [ "$PBARBATUS_PROCS" -eq 0 ]; then
    echo "✅ SUCCESS: No Pbarbatus download processes running"
    echo "   Skip logic is working correctly!"
  else
    echo "⚠️  WARNING: Found $PBARBATUS_PROCS Pbarbatus download processes"
    echo "   Skip logic may not be working correctly"
    ps aux | grep "[a]malgkit.*pbarbatus.*getfastq\|[f]astq-dump.*pbarbatus\|[f]astp.*pbarbatus" | head -3
  fi
  
  echo ""
  echo "=== BATCH FILES CHECK ==="
  BATCH_FILES=$(find output/amalgkit/pbarbatus/work/metadata -name "metadata.batch*.tsv" 2>/dev/null | wc -l | tr -d ' ')
  if [ "$BATCH_FILES" -eq 0 ]; then
    echo "✅ SUCCESS: No batch metadata files created"
    echo "   All samples were skipped before batch creation"
  else
    echo "⚠️  WARNING: Found $BATCH_FILES batch metadata files"
    echo "   Some batches were created (should be 0 if all samples skipped)"
  fi
  
else
  echo "⏳ Pbarbatus not yet reached in workflow"
  echo ""
  echo "Current progress:"
  grep "Starting workflow for" "$LOG_FILE" 2>/dev/null | tail -2
  echo ""
  echo "Last 20 lines:"
  tail -20 "$LOG_FILE" 2>/dev/null
fi

echo ""
echo "================================================================================="
echo "VERIFICATION COMPLETE"
echo "================================================================================="
echo ""
echo "Log file: $LOG_FILE"
echo "Workflow PID: $WORKFLOW_PID"
echo ""
echo "To continue monitoring:"
echo "  tail -f $LOG_FILE"
echo ""

