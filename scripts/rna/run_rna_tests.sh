#!/bin/bash
# RNA Test Suite Execution Script
# Runs all RNA tests with coverage reporting

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

echo "=========================================="
echo "RNA Test Suite Execution"
echo "=========================================="
echo ""

# Check if pytest-cov is available
if python3 -c "import pytest_cov" 2>/dev/null; then
    COV_OPTS="--cov=src/metainformant/rna/workflow --cov=src/metainformant/rna/amalgkit --cov=src/metainformant/rna/steps --cov-report=html --cov-report=term-missing"
    echo "✓ Coverage reporting enabled"
else
    COV_OPTS=""
    echo "⚠ Coverage reporting disabled (pytest-cov not installed)"
    echo "  Install with: pip install pytest-cov"
fi

echo ""
echo "Running RNA tests..."
echo ""

# Create output directory
mkdir -p output

# Run tests with PYTHONPATH set
export PYTHONPATH="$REPO_ROOT/src:$PYTHONPATH"
python3 -m pytest tests/test_rna_*.py -v --tb=short $COV_OPTS 2>&1 | tee output/test_results_rna.txt

EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ All tests passed!"
else
    echo "✗ Some tests failed (exit code: $EXIT_CODE)"
fi
echo "=========================================="
echo ""
echo "Results saved to: output/test_results_rna.txt"

if [ -n "$COV_OPTS" ]; then
    echo "Coverage report: output/coverage_html/index.html"
fi

exit $EXIT_CODE

