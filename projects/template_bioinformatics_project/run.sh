#!/usr/bin/env bash
# =============================================================================
# run.sh — Master pipeline runner for the template bioinformatics project.
# =============================================================================
# Usage:
#   ./run.sh                          # Run all stages
#   ./run.sh --stage process          # Run Stage 1 only
#   ./run.sh --stage analyze
#   ./run.sh --stage visualize
#   ./run.sh --stage all
#   ./run.sh --config config/custom.yaml
#   ./run.sh --clean                  # Remove all generated outputs
#   ./run.sh --synthetic              # Generate synthetic test data first
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CONFIG="config/default.yaml"
STAGE="all"
CLEAN=false
SYNTHETIC=false

# ── Parse arguments ────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)   CONFIG="$2"; shift 2 ;;
        --stage)    STAGE="$2";  shift 2 ;;
        --clean)    CLEAN=true;  shift   ;;
        --synthetic) SYNTHETIC=true; shift ;;
        -h|--help)
            sed -n '3,14p' "$0" | sed 's/^# \?//'
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ── Clean mode ─────────────────────────────────────────────────────────────────
if $CLEAN; then
    echo "🧹 Cleaning generated outputs..."
    rm -rf data/processed/* results/* logs/*
    # Preserve .gitkeep files
    touch data/processed/.gitkeep results/.gitkeep logs/.gitkeep 2>/dev/null || true
    echo "✅ Clean complete."
    exit 0
fi

# ── Ensure uv is available ─────────────────────────────────────────────────────
if ! command -v uv &>/dev/null; then
    echo "❌ 'uv' is required but not found. Install via: curl -LsSf https://astral.sh/uv/install.sh | sh" >&2
    exit 1
fi

# ── Synthetic data ─────────────────────────────────────────────────────────────
if $SYNTHETIC; then
    echo "🔬 Generating synthetic data..."
    uv run scripts/99_create_synthetic_data.py --config "$CONFIG"
fi

# ── Stage runners ──────────────────────────────────────────────────────────────
run_stage() {
    local name="$1"
    local script="$2"
    echo ""
    echo "==> Stage: $name"
    uv run "$script" --config "$CONFIG"
}

run_process()   { run_stage "Data Processing        (01)" scripts/01_process_data.py; }
run_analyze()   { run_stage "Results Analysis       (02)" scripts/02_analyze_results.py; }
run_visualize() { run_stage "Visualization          (03)" scripts/03_visualize.py; }

# ── Dispatch ───────────────────────────────────────────────────────────────────
echo "🚀 Starting pipeline (stage=${STAGE}, config=${CONFIG})"

case "$STAGE" in
    process)   run_process ;;
    analyze)   run_analyze ;;
    visualize) run_visualize ;;
    all)
        run_process
        run_analyze
        run_visualize
        ;;
    *)
        echo "❌ Unknown stage: $STAGE" >&2
        exit 1
        ;;
esac

echo ""
echo "✨ Pipeline finished successfully."
echo "   📂 Raw data:     data/raw/"
echo "   📂 Processed:    data/processed/"
echo "   📂 Results:      results/"
echo "   📂 Logs:         logs/"
