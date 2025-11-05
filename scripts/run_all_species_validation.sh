#!/bin/bash
# Comprehensive validation and startup test for all species end-to-end workflows

set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "=================================================================================="
echo "COMPREHENSIVE END-TO-END AMALGKIT WORKFLOW VALIDATION"
echo "=================================================================================="
echo

# Step 1: Check virtual environment
echo "Step 1: Checking virtual environment..."
if [ ! -d ".venv" ]; then
    echo "  ⚠️  Virtual environment not found"
    echo "  Setting up virtual environment..."
    python3 -m venv .venv
    source .venv/bin/activate
    pip install -e .
    pip install git+https://github.com/kfuku52/amalgkit
    echo "  ✅ Virtual environment created and configured"
else
    echo "  ✅ Virtual environment exists"
    source .venv/bin/activate
fi

# Step 2: Validate all species configs
echo
echo "Step 2: Validating all species configurations..."
python3 scripts/validate_all_species_workflow.py

# Step 3: Test workflow planning
echo
echo "Step 3: Testing workflow planning for all species..."
python3 scripts/test_end_to_end_startup.py

# Step 4: Test thread configuration
echo
echo "Step 4: Testing thread configuration..."
export AK_THREADS=12
python3 -c "
import sys
sys.path.insert(0, 'src')
from metainformant.rna.workflow import load_workflow_config
from pathlib import Path

config_dir = Path('config/amalgkit')
configs = sorted(config_dir.glob('amalgkit_*.yaml'))
configs = [c for c in configs if 'template' not in c.stem.lower() and 'test' not in c.stem.lower()]

if configs:
    cfg = load_workflow_config(configs[0])
    if cfg.threads == 12:
        print(f'  ✅ AK_THREADS=12 override works (threads: {cfg.threads})')
    else:
        print(f'  ❌ AK_THREADS override failed (got: {cfg.threads})')
        sys.exit(1)
"

# Step 5: Show startup commands
echo
echo "=================================================================================="
echo "STARTUP COMMANDS"
echo "=================================================================================="
echo
echo "✅ All validation passed!"
echo
echo "To start all species with configurable threads (SRA workflow):"
echo "  export AK_THREADS=12"
echo "  python3 scripts/rna/run_multi_species.py"
echo
echo "To start all species with configurable threads (ENA workflow - recommended):"
echo "  for config in config/amalgkit/amalgkit_*.yaml; do"
echo "    [[ \"\$config\" == *template* ]] && continue"
echo "    [[ \"\$config\" == *test* ]] && continue"
echo "    species=\$(basename \"\$config\" .yaml | sed 's/amalgkit_//')"
echo "    nohup python3 scripts/rna/workflow_ena_integrated.py \\"
echo "      --config \"\$config\" \\"
echo "      --batch-size 12 \\"
echo "      --threads 12 \\"
echo "      > output/workflow_\${species}_\$(date +%Y%m%d_%H%M%S).log 2>&1 &"
echo "  done"
echo
echo "Monitor progress:"
echo "  python3 scripts/rna/monitor_comprehensive.py"
echo
echo "=================================================================================="

