#!/bin/bash
# UV-based code quality checks
set -euo pipefail
cd "$(dirname "$0")/.."

echo "🔍 Running code quality checks with UV..."
case "${1:-all}" in
    "format")
        echo "🎨 Formatting code..."
        uv run black src/ tests/
        uv run isort src/ tests/
        ;;
    "lint")
        echo "🔎 Linting code..."
        uv run flake8 src/ tests/
        ;;
    "typecheck")
        echo "🔍 Type checking..."
        uv run mypy src/
        ;;
    "all"|*)
        echo "🎨 Formatting code..."
        uv run black src/ tests/
        uv run isort src/ tests/
        echo "🔎 Linting code..."
        uv run flake8 src/ tests/
        echo "🔍 Type checking..."
        uv run mypy src/
        ;;
esac
