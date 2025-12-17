#!/bin/bash
# UV-based documentation builder
set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"

# Setup environment
setup_environment

echo "📚 Building documentation with UV..."
case "${1:-build}" in
    "build")
        uv run sphinx-build -b html docs/ docs/_build/html
        print_status "OK" "Documentation built in docs/_build/html/"
        ;;
    "serve")
        if [ -d "docs/_build/html" ]; then
            print_status "INFO" "Serving documentation at http://localhost:8000"
            uv run python -m http.server 8000 -d docs/_build/html
        else
            print_status "ERROR" "Documentation not built yet. Run: ./scripts/uv_docs.sh build"
            exit 1
        fi
        ;;
    "clean")
        rm -rf docs/_build/
        print_status "OK" "Documentation build directory cleaned"
        ;;
    *)
        echo "Usage: $0 [build|serve|clean]"
        exit 1
        ;;
esac
