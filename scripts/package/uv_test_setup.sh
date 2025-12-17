#!/bin/bash
# UV-based test environment setup for METAINFORMANT (DEPRECATED)
# This script is deprecated. Use scripts/package/setup.sh and scripts/package/verify.sh instead.
set -euo pipefail

# Source common utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_common.sh"

# Show deprecation warning
show_deprecation_warning "uv_test_setup.sh" "scripts/package/setup.sh --dev && scripts/package/verify.sh --mode deps --fix"

# Parse command line arguments
TEST_TYPE="all"
VERIFY_ONLY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --test-type)
            TEST_TYPE="${2:-all}"
            shift 2 || true
            ;;
        --verify-only)
            VERIFY_ONLY=1
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage
            exit 2
            ;;
    esac
done

# Check UV availability
check_uv() {
    echo "🔍 Checking UV availability..."
    if ! command -v uv >/dev/null 2>&1; then
        echo "❌ uv is not installed or not in PATH" >&2
        echo "   Please install uv first:" >&2
        echo "   curl -LsSf https://astral.sh/uv/install.sh | sh" >&2
        echo "   Or visit: https://github.com/astral.sh/uv" >&2
        echo ""
        echo "   After installation, add uv to your PATH:"
        echo "   export PATH=\"\$HOME/.cargo/bin:\$PATH\"" >&2
        exit 1
    fi
    echo "✅ uv $(uv --version) found"
}

# Detect filesystem type and configure UV
detect_filesystem() {
    echo "🔍 Detecting filesystem type..."
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
        echo "📁 FAT filesystem detected - using /tmp locations"
        echo "   UV Cache: $UV_CACHE_DIR"
        echo "   Virtual Env: $VENV_DIR"
    else
        VENV_DIR=".venv"
        echo "📁 Standard filesystem detected - using repo locations"
        echo "   UV Cache: ${UV_CACHE_DIR:-.uv-cache}"
        echo "   Virtual Env: $VENV_DIR"
    fi
}

# Create/update virtual environment
setup_virtual_environment() {
    echo "🏗️  Setting up virtual environment..."

    if [[ ! -d "$VENV_DIR" ]]; then
        echo "   Creating virtual environment at $VENV_DIR..."
        uv venv "$VENV_DIR"
    else
        echo "   Virtual environment already exists at $VENV_DIR"
    fi

    echo "✅ Virtual environment ready"
}

# Sync dependencies based on test type
sync_dependencies() {
    local test_type="$1"
    echo "🔄 Syncing dependencies for test type: $test_type..."

    case "$test_type" in
        "fast")
            echo "   Installing fast test dependencies..."
            uv sync --group test-fast
            ;;
        "network")
            echo "   Installing network test dependencies..."
            uv sync --group test-network
            ;;
        "external")
            echo "   Installing external tool test dependencies..."
            uv sync --group test-external
            ;;
        "all")
            echo "   Installing all test dependencies..."
            uv sync --group test-all
            ;;
        *)
            echo "❌ Unknown test type: $test_type" >&2
            exit 1
            ;;
    esac

    echo "✅ Dependencies synced"
}

# Verify test environment
verify_test_environment() {
    local test_type="$1"
    echo "🔍 Verifying test environment for type: $test_type..."

    # Check basic Python availability
    if ! uv run python --version >/dev/null 2>&1; then
        echo "❌ Python not available in virtual environment"
        return 1
    fi
    echo "   ✅ Python available"

    # Check pytest availability
    if ! uv run python -c "import pytest; print('pytest available')" >/dev/null 2>&1; then
        echo "❌ pytest not available"
        return 1
    fi
    echo "   ✅ pytest available"

    # Check test data directory
    if [[ ! -d "tests/data" ]]; then
        echo "⚠️  Test data directory not found at tests/data"
        echo "   Some tests may fail - run from repository root"
    else
        echo "   ✅ Test data directory available"
    fi

    # Check for output directory
    if [[ ! -d "output" ]]; then
        echo "   Creating output directory..."
        mkdir -p output
    fi
    echo "   ✅ Output directory available"

    # Type-specific checks
    case "$test_type" in
        "network")
            # Check network connectivity
            if ! uv run python -c "
import requests
try:
    response = requests.get('https://httpbin.org/status/200', timeout=5)
    print('Network available')
except:
    raise Exception('Network not available')
" >/dev/null 2>&1; then
                echo "⚠️  Network connectivity check failed"
                echo "   Network tests may be skipped"
            else
                echo "   ✅ Network connectivity available"
            fi
            ;;
        "external")
            # Check for common external tools
            echo "   Checking external tools..."
            tools=("amalgkit" "muscle")
            for tool in "${tools[@]}"; do
                if command -v "$tool" >/dev/null 2>&1; then
                    echo "   ✅ $tool available"
                else
                    echo "   ⚠️  $tool not available - some tests may be skipped"
                fi
            done
            ;;
    esac

    echo "✅ Test environment verification complete"
    return 0
}

# Set up test data if needed
setup_test_data() {
    echo "📁 Setting up test data..."

    # Ensure test data directory exists
    mkdir -p tests/data

    # Basic check for essential test files
    if [[ ! -f "tests/data/sample_dna.fasta" ]]; then
        echo "   Creating sample test data..."
        cat > tests/data/sample_dna.fasta << 'EOF'
>sample_seq_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>sample_seq_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
EOF
    fi

    echo "✅ Test data ready"
}

# Configure environment variables
configure_environment() {
    echo "⚙️  Configuring environment variables..."

    # Set UV environment variables
    if [[ -n "${UV_CACHE_DIR:-}" ]]; then
        export UV_CACHE_DIR
        echo "   UV_CACHE_DIR=$UV_CACHE_DIR"
    fi

    # Ensure PYTHONPATH includes src
    export PYTHONPATH="${PWD}/src:${PYTHONPATH:-}"
    echo "   PYTHONPATH includes src/"

    # Set test-specific environment variables
    export PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
    export PYTHONUNBUFFERED=1

    echo "✅ Environment configured"
}

# Main execution
main() {
    echo "🧪 METAINFORMANT Test Environment Setup"
    echo "======================================"
    echo ""

    # Check prerequisites
    check_uv
    detect_filesystem

    if [[ $VERIFY_ONLY -eq 1 ]]; then
        echo "🔍 Running verification only..."
        verify_test_environment "$TEST_TYPE"
        echo ""
        echo "✅ Verification complete!"
        exit 0
    fi

    # Setup steps
    setup_virtual_environment
    sync_dependencies "$TEST_TYPE"
    setup_test_data
    configure_environment

    # Final verification
    echo ""
    echo "🔍 Running final verification..."
    if verify_test_environment "$TEST_TYPE"; then
        echo ""
        echo "🎉 Test environment setup complete!"
        echo ""
        echo "Run tests with:"
        echo "  bash scripts/package/uv_test.sh $TEST_TYPE"
        echo "  bash scripts/package/uv_test_optimized.sh fast"
        echo ""
        echo "Verify setup with:"
        echo "  bash scripts/package/uv_test_setup.sh --verify-only"
    else
        echo ""
        echo "❌ Test environment setup failed!"
        exit 1
    fi
}

# Run main function
main "$@"


