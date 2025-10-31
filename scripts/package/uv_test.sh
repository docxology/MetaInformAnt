#!/bin/bash
# UV-based test runner with comprehensive coverage
set -euo pipefail
cd "$(dirname "$0")/.."

echo "🧪 Running tests with UV..."
case "${1:-all}" in
    "fast")
        uv run pytest -x -q --tb=short
        ;;
    "coverage")
        uv run pytest --cov=src/metainformant --cov-report=term-missing --cov-report=html:output/coverage_html
        ;;
    "parallel")
        uv run pytest -n auto --tb=short
        ;;
    "benchmark")
        uv run pytest tests/ -m benchmark --benchmark-json=output/benchmarks/results.json
        ;;
    "network")
        uv run pytest tests/ -m network --tb=short
        ;;
    "external")
        uv run pytest tests/ -m external_tool --tb=short
        ;;
    "integration")
        uv run pytest tests/ -m integration --tb=short
        ;;
    "all"|*)
        uv run pytest --cov=src/metainformant --cov-report=term-missing --cov-report=html:output/coverage_html --cov-report=xml:output/coverage.xml
        ;;
esac
