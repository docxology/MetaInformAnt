#!/bin/bash
# UV-based code quality checks
set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"

# Setup environment
setup_environment

print_status "INFO" "Running code quality checks with UV..."
case "${1:-all}" in
    "format")
        print_status "INFO" "Formatting code..."
        uv run black src/ tests/
        uv run isort src/ tests/
        print_status "OK" "Code formatting completed"
        ;;
    "lint")
        print_status "INFO" "Linting code..."
        uv run flake8 src/ tests/
        print_status "OK" "Code linting completed"
        ;;
    "typecheck")
        print_status "INFO" "Type checking..."
        uv run mypy src/
        print_status "OK" "Type checking completed"
        ;;
    "all"|*)
        print_status "INFO" "Formatting code..."
        uv run black src/ tests/
        uv run isort src/ tests/
        print_status "OK" "Code formatting completed"

        print_status "INFO" "Linting code..."
        uv run flake8 src/ tests/
        print_status "OK" "Code linting completed"

        print_status "INFO" "Type checking..."
        uv run mypy src/
        print_status "OK" "Type checking completed"
        ;;
esac
