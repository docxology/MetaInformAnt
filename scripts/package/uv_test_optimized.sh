#!/bin/bash
# Optimized UV test runner with speed enhancements for METAINFORMANT
set -euo pipefail
cd "$(dirname "$0")/.."

# Detect filesystem type and configure UV cache/venv
# Temporarily disable exit on error to capture Python script exit code
set +e
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

if not has_symlinks:
    cache_dir = "/tmp/uv-cache"
    venv_dir = "/tmp/metainformant_venv"
    print(f"FAT:{cache_dir}:{venv_dir}")
    sys.exit(1)
else:
    cache_dir = ".uv-cache"
    venv_dir = ".venv"
    print(f"STANDARD:{cache_dir}:{venv_dir}")
    sys.exit(0)
PYTHON
)

FS_EXIT_CODE=$?
set -e

if [[ $FS_EXIT_CODE -eq 1 ]]; then
  # FAT filesystem detected
  export UV_CACHE_DIR="/tmp/uv-cache"
  VENV_DIR="/tmp/metainformant_venv"
  mkdir -p "$UV_CACHE_DIR"
  # Check if venv exists, if not try .venv
  if [[ ! -d "$VENV_DIR" ]] && [[ -d ".venv" ]]; then
    VENV_DIR=".venv"
  fi
  USE_VENV_PYTHON=1
else
  VENV_DIR=".venv"
  USE_VENV_PYTHON=0
fi

# Determine pytest command
if [[ $USE_VENV_PYTHON -eq 1 ]]; then
  if [[ -d "$VENV_DIR" ]] && [[ -f "$VENV_DIR/bin/pytest" ]]; then
    # Use venv pytest directly on FAT filesystems
    PYTEST_CMD="$VENV_DIR/bin/pytest"
  elif [[ -d "$VENV_DIR" ]] && [[ -f "$VENV_DIR/bin/python3" ]]; then
    # Use venv Python with -m pytest if pytest not directly available
    PYTEST_CMD="$VENV_DIR/bin/python3 -m pytest"
  else
    # FAT filesystem but venv doesn't exist - provide helpful error
    echo "❌ Virtual environment not found at $VENV_DIR" >&2
    echo "   FAT filesystem detected - venv must be created in /tmp" >&2
    echo "   Run: bash scripts/package/setup_uv.sh" >&2
    exit 1
  fi
else
  # Standard filesystem - use uv run
  PYTEST_CMD="uv run pytest"
fi

echo "🚀 Running optimized tests with UV..."
echo "  → Filesystem: $(python3 -c "from pathlib import Path; import subprocess; result = subprocess.run(['df', '-T', str(Path.cwd())], capture_output=True, text=True, timeout=5); print(result.stdout.split('\n')[1].split()[1] if result.returncode == 0 else 'unknown')" 2>/dev/null || echo "unknown")"
echo "  → UV cache: ${UV_CACHE_DIR:-.uv-cache}"
echo "  → Virtual env: $VENV_DIR"

case "${1:-fast}" in
    "ultra-fast")
        echo "⚡ Ultra-fast test mode (core functionality only)"
        $PYTEST_CMD tests/test_core_functionality.py -x -q --tb=no --disable-warnings --timeout=5
        ;;
    "fast")
        echo "🏃 Fast test mode (essential tests only)"
        $PYTEST_CMD tests/test_core_functionality.py tests/test_core_comprehensive.py tests/test_dna_comprehensive.py -x -q --tb=short --timeout=10
        ;;
    "integration")
        echo "🔗 Integration tests (optimized for speed)"
        $PYTEST_CMD tests/test_integration_comprehensive.py -v --tb=short --timeout=30 --disable-warnings
        ;;
    "coverage-fast")
        echo "📊 Fast coverage analysis (core modules only)"
        $PYTEST_CMD tests/test_core_functionality.py tests/test_core_*.py \
            --cov=src/metainformant/core --cov=src/metainformant/simulation \
            --cov-report=term-missing --cov-report=html:output/coverage_html_fast \
            --cov-fail-under=0 --timeout=15
        ;;
    "smoke")
        echo "💨 Smoke tests (basic functionality verification)"
        $PYTEST_CMD tests/test_core_functionality.py::TestCoreIO::test_json_operations \
            tests/test_core_functionality.py::TestDNAAnalysis::test_gc_content_calculation \
            tests/test_core_functionality.py::TestSimulation::test_sequence_generation \
            -v --tb=short --timeout=5
        ;;
    "network-only")
        echo "🕸️  Network tests only"
        $PYTEST_CMD tests/test_networks_*.py -v --tb=short --timeout=20
        ;;
    "ml-only")
        echo "🤖 Machine learning tests only"
        $PYTEST_CMD tests/test_ml_*.py -v --tb=short --timeout=20
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
