#!/bin/bash
# UV-based test runner with comprehensive coverage
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

echo "🧪 Running tests with UV..."
echo "  → Filesystem: $(python3 -c "from pathlib import Path; import subprocess; result = subprocess.run(['df', '-T', str(Path.cwd())], capture_output=True, text=True, timeout=5); print(result.stdout.split('\n')[1].split()[1] if result.returncode == 0 else 'unknown')" 2>/dev/null || echo "unknown")"
echo "  → UV cache: ${UV_CACHE_DIR:-.uv-cache}"
echo "  → Virtual env: $VENV_DIR"
case "${1:-all}" in
    "fast")
        $PYTEST_CMD -x -q --tb=short
        ;;
    "coverage")
        $PYTEST_CMD --cov=src/metainformant --cov-report=term-missing --cov-report=html:output/coverage_html
        ;;
    "parallel")
        $PYTEST_CMD -n auto --tb=short
        ;;
    "benchmark")
        $PYTEST_CMD tests/ -m benchmark --benchmark-json=output/benchmarks/results.json
        ;;
    "network")
        $PYTEST_CMD tests/ -m network --tb=short
        ;;
    "external")
        $PYTEST_CMD tests/ -m external_tool --tb=short
        ;;
    "integration")
        $PYTEST_CMD tests/ -m integration --tb=short
        ;;
    "all"|*)
        $PYTEST_CMD --cov=src/metainformant --cov-report=term-missing --cov-report=html:output/coverage_html --cov-report=xml:output/coverage.xml
        ;;
esac
