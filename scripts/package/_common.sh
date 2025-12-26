#!/bin/bash
# METAINFORMANT Common Utilities Library
# Shared functions for package scripts
# This file should be sourced by other scripts, not executed directly

set -euo pipefail

# Color codes for output (standardized across all scripts)
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Status printing functions
print_status() {
    local status="$1"
    local message="$2"
    case "$status" in
        "OK")
            echo -e "${GREEN}✅${NC} $message"
            ;;
        "WARN")
            echo -e "${YELLOW}⚠️${NC}  $message"
            ;;
        "ERROR")
            echo -e "${RED}❌${NC} $message"
            ;;
        "INFO")
            echo -e "${BLUE}ℹ️${NC}  $message"
            ;;
        "FIXED")
            echo -e "${GREEN}🔧${NC} $message"
            ;;
        *)
            echo "$message"
            ;;
    esac
}

# Check UV availability
check_uv() {
    if ! command -v uv >/dev/null 2>&1; then
        print_status "ERROR" "uv is not installed or not in PATH"
        echo "       Install uv from: https://github.com/astral-sh/uv"
        echo "       curl -LsSf https://astral.sh/uv/install.sh | sh"
        return 1
    fi

    local uv_version
    uv_version=$(uv --version 2>/dev/null || echo "unknown")
    print_status "OK" "uv available: $uv_version"
    return 0
}

# Detect filesystem type and symlink support
# Returns: sets FS_TYPE and HAS_SYMLINKS variables
# Exit code: 0 for standard filesystem, 1 for FAT filesystem
detect_filesystem() {
    # Temporarily disable exit on error to capture Python script exit code
    set +e
    local FS_DETECT
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

print(f"{fs_type}:{has_symlinks}")
PYTHON
)

    local FS_EXIT_CODE=$?
    set -e

    FS_TYPE=$(echo "$FS_DETECT" | cut -d: -f1)
    HAS_SYMLINKS=$(echo "$FS_DETECT" | cut -d: -f2)

    if [[ "$HAS_SYMLINKS" == "False" ]]; then
        return 1  # FAT filesystem
    else
        return 0  # Standard filesystem
    fi
}

# Configure UV cache directory based on filesystem
# Sets UV_CACHE_DIR environment variable
configure_uv_cache() {
    if detect_filesystem; then
        # Standard filesystem
        UV_CACHE_DIR="${UV_CACHE_DIR:-.uv-cache}"
        print_status "INFO" "Using standard filesystem - cache: $UV_CACHE_DIR"
    else
        # FAT filesystem
        UV_CACHE_DIR="/tmp/uv-cache"
        export UV_CACHE_DIR
        mkdir -p "$UV_CACHE_DIR"
        print_status "WARN" "FAT filesystem detected - cache: $UV_CACHE_DIR"
    fi
}

# Get appropriate virtual environment directory
get_venv_dir() {
    if detect_filesystem; then
        # Standard filesystem
        echo ".venv"
    else
        # FAT filesystem
        echo "/tmp/metainformant_venv"
    fi
}

# Get appropriate pytest command based on filesystem
get_pytest_cmd() {
    local venv_dir
    venv_dir=$(get_venv_dir)

    if [[ -d "$venv_dir" ]]; then
        if [[ -f "$venv_dir/bin/pytest" ]]; then
            echo "$venv_dir/bin/pytest"
        elif [[ -f "$venv_dir/bin/python3" ]]; then
            echo "$venv_dir/bin/python3 -m pytest"
        else
            # Fallback to uv run
            echo "uv run pytest"
        fi
    else
        echo "uv run pytest"
    fi
}

# Sync dependencies based on test type
sync_dependencies() {
    local test_type="${1:-fast}"
    print_status "INFO" "Syncing dependencies for test type: $test_type"

    case "$test_type" in
        "fast")
            uv sync --group test-fast
            ;;
        "network")
            uv sync --group test-network
            ;;
        "external")
            uv sync --group test-external
            ;;
        "all")
            uv sync --group test-all --extra all
            ;;
        *)
            uv sync --group test-fast
            ;;
    esac

    print_status "OK" "Dependencies synced"
}

# Get repository root directory
get_repo_root() {
    local script_dir
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    # Go up from scripts/package to repo root
    cd "$script_dir/../.." && pwd
}

# Setup basic environment
setup_environment() {
    local repo_root
    repo_root=$(get_repo_root)
    cd "$repo_root"

    # Ensure output directory exists
    mkdir -p output

    # Configure UV cache
    configure_uv_cache

    # Set PYTHONPATH
    export PYTHONPATH="${repo_root}/src:${PYTHONPATH:-}"
}

# Show deprecation warning and redirect to new command
show_deprecation_warning() {
    local old_script="$1"
    local new_command="$2"

    print_status "WARN" "Script '$old_script' is deprecated"
    echo "       Please use: $new_command"
    echo "       This script will be removed in a future version."
    echo ""

    # Execute the new command
    eval "$new_command"
}

# Export functions so they can be used when sourced
export -f print_status
export -f check_uv
export -f detect_filesystem
export -f configure_uv_cache
export -f get_venv_dir
export -f get_pytest_cmd
export -f sync_dependencies
export -f get_repo_root
export -f setup_environment
export -f show_deprecation_warning
