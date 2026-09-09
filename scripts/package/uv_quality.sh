#!/bin/bash
# UV-based code quality checks
#
# Mirrors the blocking quality-checks job in .github/workflows/test.yml so
# developers can reproduce the CI gate locally before pushing.
set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"

# Setup environment
setup_environment

print_status "INFO" "Running code quality checks with UV..."
run_format() {
    print_status "INFO" "Formatting code..."
    uv run black src/ tests/ scripts/
    uv run isort src/ tests/ scripts/
    print_status "OK" "Code formatting completed"
}

run_lint() {
    print_status "INFO" "Linting code..."
    uv run flake8 src/ tests/ --max-line-length=120 --extend-ignore=E203,W503
    # Script entrypoints intentionally bootstrap src before package imports;
    # Black/isort cover their formatting while F/E9 remains blocking.
    uv run flake8 scripts/ --select=F,E9 --max-line-length=120 --extend-ignore=E203,W503
    print_status "OK" "Code linting completed"
}

run_typecheck() {
    print_status "INFO" "Type checking..."
    uv run python scripts/quality/check_mypy_budget.py src/metainformant/rna
    print_status "OK" "Type checking completed"
}

case "${1:-all}" in
    "format")
        run_format
        ;;
    "lint")
        run_lint
        ;;
    "typecheck")
        run_typecheck
        ;;
    "all"|*)
        run_format
        run_lint
        run_typecheck
        print_status "OK" "Code quality checks completed"
        ;;
esac
