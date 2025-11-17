#!/bin/bash
# Comprehensive UV development setup script for METAINFORMANT
# Follows .cursorrules and NO_MOCKING_POLICY for real implementation testing

set -euo pipefail

echo "🚀 Setting up METAINFORMANT development environment with UV..."

# Ensure we're in the project root
cd "$(dirname "$0")/.."
PROJECT_ROOT="$(pwd)"

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "❌ UV is not installed. Please install it first:"
    echo "   curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
fi

echo "✅ UV version: $(uv --version)"

# Detect filesystem type and configure UV cache
echo "🔍 Detecting filesystem type..."
FS_DETECT=$(python3 <<'PYTHON'
import subprocess
import sys
from pathlib import Path

def detect_fs_type(path):
    """Detect filesystem type using df -T."""
    try:
        result = subprocess.run(
            ["df", "-T", str(path)],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            lines = result.stdout.strip().split("\n")
            if len(lines) > 1:
                parts = lines[1].split()
                if len(parts) >= 2:
                    return parts[1].lower()
    except Exception:
        pass
    return "unknown"

def supports_symlinks(fs_type):
    """Check if filesystem supports symlinks."""
    no_symlink_fs = {"exfat", "fat32", "fat", "vfat", "msdos"}
    return fs_type not in no_symlink_fs

repo_root = Path.cwd()
fs_type = detect_fs_type(repo_root)
has_symlinks = supports_symlinks(fs_type)

print(f"Filesystem type: {fs_type}")
print(f"Symlink support: {has_symlinks}")

if not has_symlinks:
    cache_dir = "/tmp/uv-cache"
    print(f"FAT filesystem detected - using {cache_dir} for UV cache")
    sys.exit(1)  # Signal to use alternative cache location
else:
    cache_dir = ".uv-cache"
    print(f"Using standard cache location: {cache_dir}")
    sys.exit(0)
PYTHON
)

FS_EXIT_CODE=$?

if [[ $FS_EXIT_CODE -eq 1 ]]; then
  # FAT filesystem detected - use /tmp/uv-cache
  export UV_CACHE_DIR="/tmp/uv-cache"
  mkdir -p "$UV_CACHE_DIR"
  echo "  → UV cache directory set to: $UV_CACHE_DIR"
else
  echo "  → Using standard UV cache directory: .uv-cache"
fi

# Create uv.lock if it doesn't exist and sync dependencies
echo "📦 Syncing dependencies with UV..."
uv sync --extra dev --extra scientific --extra all

echo "🔧 Setting up development tools..."

# Install pre-commit hooks
echo "📋 Setting up pre-commit hooks..."
uv run pre-commit install

# Create output directory structure following .cursorrules
echo "📁 Creating output directory structure..."
mkdir -p output/{coverage_html,benchmarks,profiles,docs,test_results}
mkdir -p output/{dna,rna,protein,math,networks,singlecell,visualization}
mkdir -p output/{simulation,ontology,quality,multiomics,epigenome}

# Set up documentation build environment
echo "📚 Setting up documentation environment..."
mkdir -p docs/_build docs/_static docs/_templates

# Create UV task runner wrapper scripts
echo "⚙️  Creating UV task runner scripts..."

# Test runner script
cat > scripts/uv_test.sh << 'EOF'
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
EOF
chmod +x scripts/uv_test.sh

# Code quality script
cat > scripts/uv_quality.sh << 'EOF'
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
EOF
chmod +x scripts/uv_quality.sh

# Documentation builder script
cat > scripts/uv_docs.sh << 'EOF'
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
EOF
chmod +x scripts/uv_docs.sh

# Performance profiling script
cat > scripts/uv_profile.sh << 'EOF'
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
EOF
chmod +x scripts/uv_profile.sh

echo "✅ Development environment setup complete!"
echo ""
echo "🎯 Available UV workflows:"
echo "  ./scripts/uv_test.sh [fast|coverage|parallel|benchmark|network|external|integration|all]"
echo "  ./scripts/uv_quality.sh [format|lint|typecheck|all]"
echo "  ./scripts/uv_docs.sh [build|serve|clean]"
echo "  ./scripts/uv_profile.sh [cpu|memory|benchmark] [command...]"
echo ""
echo "📋 Quick start commands:"
echo "  uv run python -m metainformant --help"
echo "  ./scripts/uv_test.sh fast"
echo "  ./scripts/uv_quality.sh all"
echo "  ./scripts/uv_docs.sh build"
echo ""
echo "🚀 Happy coding with UV!"
