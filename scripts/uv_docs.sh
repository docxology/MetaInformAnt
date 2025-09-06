#!/bin/bash
# UV-based documentation builder
set -euo pipefail
cd "$(dirname "$0")/.."

echo "📚 Building documentation with UV..."
case "${1:-build}" in
    "build")
        uv run sphinx-build -b html docs/ docs/_build/html
        echo "📖 Documentation built in docs/_build/html/"
        ;;
    "serve")
        if [ -d "docs/_build/html" ]; then
            echo "🌐 Serving documentation at http://localhost:8000"
            uv run python -m http.server 8000 -d docs/_build/html
        else
            echo "❌ Documentation not built yet. Run: ./scripts/uv_docs.sh build"
            exit 1
        fi
        ;;
    "clean")
        rm -rf docs/_build/
        echo "🧹 Documentation build directory cleaned"
        ;;
    *)
        echo "Usage: $0 [build|serve|clean]"
        exit 1
        ;;
esac
