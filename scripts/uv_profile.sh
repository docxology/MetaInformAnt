#!/bin/bash
# UV-based performance profiling
set -euo pipefail
cd "$(dirname "$0")/.."

echo "⚡ Performance profiling with UV..."
case "${1:-help}" in
    "cpu")
        shift
        echo "🔥 CPU profiling: ${*:-python -c 'import metainformant; print(metainformant.__version__)'}"
        uv run python -m cProfile -o output/profiles/cpu_profile.stats "${@:-python -c 'import metainformant; print(metainformant.__version__)'}"
        echo "📊 CPU profile saved to output/profiles/cpu_profile.stats"
        ;;
    "memory")
        shift
        echo "💾 Memory profiling: ${*}"
        uv run python -m memory_profiler "${@}"
        ;;
    "benchmark")
        echo "🏁 Running performance benchmarks..."
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
