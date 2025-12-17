#!/bin/bash
# UV-based performance profiling
set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"

# Setup environment
setup_environment

print_status "INFO" "Performance profiling with UV..."
case "${1:-help}" in
    "cpu")
        shift
        print_status "INFO" "CPU profiling: ${*:-python -c 'import metainformant; print(metainformant.__version__)'}"
        uv run python -m cProfile -o output/profiles/cpu_profile.stats "${@:-python -c 'import metainformant; print(metainformant.__version__)'}"
        print_status "OK" "CPU profile saved to output/profiles/cpu_profile.stats"
        ;;
    "memory")
        shift
        print_status "INFO" "Memory profiling: ${*}"
        uv run python -m memory_profiler "${@}"
        ;;
    "benchmark")
        print_status "INFO" "Running performance benchmarks..."
        uv run pytest tests/ -m benchmark --benchmark-json=output/benchmarks/results.json
        ;;
    *)
        echo "Usage: $0 [cpu|memory|benchmark] [command...]"
        echo "Examples:"
        echo "  $0 cpu python -m metainformant.dna.sequences --help"
        echo "  $0 memory tests/test_dna_phylogeny.py"
        echo "  $0 benchmark"
        ;;
esac
