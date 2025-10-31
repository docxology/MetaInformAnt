#!/bin/bash
# Optimized UV test runner with speed enhancements for METAINFORMANT
set -euo pipefail
cd "$(dirname "$0")/.."

echo "🚀 Running optimized tests with UV..."

case "${1:-fast}" in
    "ultra-fast")
        echo "⚡ Ultra-fast test mode (core functionality only)"
        uv run pytest tests/test_core_functionality.py -x -q --tb=no --disable-warnings --timeout=5
        ;;
    "fast")
        echo "🏃 Fast test mode (essential tests only)"
        uv run pytest tests/test_core_functionality.py tests/test_core_comprehensive.py tests/test_dna_comprehensive.py -x -q --tb=short --timeout=10
        ;;
    "integration")
        echo "🔗 Integration tests (optimized for speed)"
        uv run pytest tests/test_integration_comprehensive.py -v --tb=short --timeout=30 --disable-warnings
        ;;
    "coverage-fast")
        echo "📊 Fast coverage analysis (core modules only)"
        uv run pytest tests/test_core_functionality.py tests/test_core_*.py \
            --cov=src/metainformant/core --cov=src/metainformant/simulation \
            --cov-report=term-missing --cov-report=html:output/coverage_html_fast \
            --cov-fail-under=0 --timeout=15
        ;;
    "smoke")
        echo "💨 Smoke tests (basic functionality verification)"
        uv run pytest tests/test_core_functionality.py::TestCoreIO::test_json_operations \
            tests/test_core_functionality.py::TestDNAAnalysis::test_gc_content_calculation \
            tests/test_core_functionality.py::TestSimulation::test_sequence_generation \
            -v --tb=short --timeout=5
        ;;
    "network-only")
        echo "🕸️  Network tests only"
        uv run pytest tests/test_networks_*.py -v --tb=short --timeout=20
        ;;
    "ml-only")
        echo "🤖 Machine learning tests only"
        uv run pytest tests/test_ml_*.py -v --tb=short --timeout=20
        ;;
    "help"|*)
        echo "Usage: $0 [mode]"
        echo "Modes:"
        echo "  ultra-fast    - Core functionality only (~5s)"
        echo "  fast          - Essential tests (~15s)"
        echo "  integration   - Integration tests (optimized, ~30s)"
        echo "  coverage-fast - Fast coverage on core modules (~20s)"
        echo "  smoke         - Basic smoke tests (~3s)"
        echo "  network-only  - Network analysis tests only"
        echo "  ml-only       - ML/features tests only"
        echo "  help          - Show this help"
        ;;
esac

echo "✅ Test run completed!"
