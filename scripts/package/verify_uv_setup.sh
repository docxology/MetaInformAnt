#!/bin/bash
# Verification script for UV setup on FAT filesystems
# Checks uv installation, filesystem detection, cache configuration, and venv handling

set -euo pipefail

cd "$(dirname "$0")/.."

echo "🔍 Verifying UV setup for METAINFORMANT..."
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

PASSED=0
FAILED=0
WARNINGS=0

check_pass() {
    echo -e "${GREEN}✓${NC} $1"
    ((PASSED++))
}

check_fail() {
    echo -e "${RED}✗${NC} $1"
    ((FAILED++))
}

check_warn() {
    echo -e "${YELLOW}⚠${NC} $1"
    ((WARNINGS++))
}

# 1. Check if uv is installed
echo "[1/7] Checking UV installation..."
if command -v uv >/dev/null 2>&1; then
    UV_VERSION=$(uv --version)
    check_pass "UV is installed: $UV_VERSION"
else
    check_fail "UV is not installed"
    echo "  Install with: curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
fi

# 2. Detect filesystem type
echo ""
echo "[2/7] Detecting filesystem type..."
FS_INFO=$(python3 <<'PYTHON'
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

FS_TYPE=$(echo "$FS_INFO" | cut -d: -f1)
HAS_SYMLINKS=$(echo "$FS_INFO" | cut -d: -f2)

if [[ "$HAS_SYMLINKS" == "False" ]]; then
    check_warn "FAT filesystem detected: $FS_TYPE (no symlink support)"
    IS_FAT=1
else
    check_pass "Filesystem: $FS_TYPE (symlink support: yes)"
    IS_FAT=0
fi

# 3. Check UV cache directory configuration
echo ""
echo "[3/7] Checking UV cache directory configuration..."
if [[ -n "${UV_CACHE_DIR:-}" ]]; then
    CACHE_DIR="$UV_CACHE_DIR"
    check_pass "UV_CACHE_DIR is set: $CACHE_DIR"
else
    if [[ $IS_FAT -eq 1 ]]; then
        CACHE_DIR="/tmp/uv-cache"
        check_warn "UV_CACHE_DIR not set, but FAT filesystem detected"
        check_warn "Should be set to: $CACHE_DIR"
    else
        CACHE_DIR=".uv-cache"
        check_pass "UV_CACHE_DIR not set (using default: $CACHE_DIR)"
    fi
fi

# Check if cache directory is writable
if mkdir -p "$CACHE_DIR" 2>/dev/null && [[ -w "$CACHE_DIR" ]]; then
    check_pass "Cache directory is writable: $CACHE_DIR"
else
    check_fail "Cache directory is not writable: $CACHE_DIR"
fi

# 4. Check virtual environment location
echo ""
echo "[4/7] Checking virtual environment location..."
if [[ $IS_FAT -eq 1 ]]; then
    EXPECTED_VENV="/tmp/metainformant_venv"
    ALT_VENV=".venv"
else
    EXPECTED_VENV=".venv"
    ALT_VENV="/tmp/metainformant_venv"
fi

if [[ -d "$EXPECTED_VENV" ]]; then
    if [[ -f "$EXPECTED_VENV/bin/python3" ]]; then
        PYTHON_VERSION=$("$EXPECTED_VENV/bin/python3" --version 2>&1 || echo "unknown")
        check_pass "Virtual environment found at: $EXPECTED_VENV ($PYTHON_VERSION)"
        VENV_DIR="$EXPECTED_VENV"
        VENV_EXISTS=1
    else
        check_warn "Directory exists but Python not found: $EXPECTED_VENV"
        VENV_EXISTS=0
    fi
elif [[ -d "$ALT_VENV" ]] && [[ $IS_FAT -eq 1 ]]; then
    check_warn "Using alternative venv location: $ALT_VENV (expected $EXPECTED_VENV for FAT)"
    if [[ -f "$ALT_VENV/bin/python3" ]]; then
        VENV_DIR="$ALT_VENV"
        VENV_EXISTS=1
    else
        VENV_EXISTS=0
    fi
else
    check_warn "Virtual environment not found (expected: $EXPECTED_VENV)"
    VENV_EXISTS=0
fi

# 5. Test venv creation (if needed)
echo ""
echo "[5/7] Testing venv creation capability..."
if [[ $VENV_EXISTS -eq 0 ]]; then
    TEST_VENV="/tmp/uv_verify_test_venv"
    rm -rf "$TEST_VENV" 2>/dev/null || true
    
    if [[ $IS_FAT -eq 1 ]]; then
        export UV_CACHE_DIR="/tmp/uv-cache"
    fi
    
    if uv venv "$TEST_VENV" >/dev/null 2>&1; then
        if [[ -f "$TEST_VENV/bin/python3" ]]; then
            check_pass "Venv creation test passed"
            rm -rf "$TEST_VENV"
        else
            check_fail "Venv created but Python not found"
        fi
    else
        check_fail "Venv creation test failed"
        if [[ $IS_FAT -eq 1 ]]; then
            check_warn "This may be due to FAT filesystem limitations"
            check_warn "Try creating venv in /tmp: uv venv /tmp/metainformant_venv"
        fi
    fi
else
    check_pass "Venv already exists, skipping creation test"
fi

# 6. Test basic uv commands
echo ""
echo "[6/7] Testing basic UV commands..."
if [[ $VENV_EXISTS -eq 1 ]]; then
    PYTHON_CMD="$VENV_DIR/bin/python3"
else
    PYTHON_CMD="python3"
fi

# Test uv pip list
if [[ $IS_FAT -eq 1 ]]; then
    export UV_CACHE_DIR="/tmp/uv-cache"
fi

if uv pip list --python "$PYTHON_CMD" >/dev/null 2>&1; then
    check_pass "uv pip list command works"
else
    check_warn "uv pip list command failed (may need venv setup)"
fi

# Test uv run
if [[ $VENV_EXISTS -eq 1 ]]; then
    if "$VENV_DIR/bin/python3" -c "import sys; sys.exit(0)" 2>/dev/null; then
        check_pass "Python in venv is executable"
    else
        check_fail "Python in venv is not executable"
    fi
else
    if uv run python -c "import sys; sys.exit(0)" >/dev/null 2>&1; then
        check_pass "uv run python command works"
    else
        check_warn "uv run python command failed (venv may need setup)"
    fi
fi

# 7. Check package installation capability
echo ""
echo "[7/7] Testing package installation capability..."
TEST_PACKAGE="packaging"  # Small, common package

if [[ $VENV_EXISTS -eq 1 ]]; then
    INSTALL_CMD=("uv" "pip" "install" "$TEST_PACKAGE" "--python" "$VENV_DIR/bin/python3")
else
    INSTALL_CMD=("uv" "pip" "install" "$TEST_PACKAGE")
fi

# Check if package is already installed
if [[ $VENV_EXISTS -eq 1 ]]; then
    if "$VENV_DIR/bin/python3" -c "import packaging" 2>/dev/null; then
        check_pass "Test package already installed"
        SKIP_INSTALL=1
    else
        SKIP_INSTALL=0
    fi
else
    SKIP_INSTALL=0
fi

if [[ $SKIP_INSTALL -eq 0 ]]; then
    if "${INSTALL_CMD[@]}" >/dev/null 2>&1; then
        check_pass "Package installation test passed"
        # Uninstall test package
        if [[ $VENV_EXISTS -eq 1 ]]; then
            uv pip uninstall -y "$TEST_PACKAGE" --python "$VENV_DIR/bin/python3" >/dev/null 2>&1 || true
        else
            uv pip uninstall -y "$TEST_PACKAGE" >/dev/null 2>&1 || true
        fi
    else
        check_warn "Package installation test failed (may need venv setup)"
        if [[ $IS_FAT -eq 1 ]]; then
            check_warn "Ensure UV_CACHE_DIR=/tmp/uv-cache is set"
        fi
    fi
fi

# Summary
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Verification Summary:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo -e "${GREEN}Passed:${NC} $PASSED"
if [[ $WARNINGS -gt 0 ]]; then
    echo -e "${YELLOW}Warnings:${NC} $WARNINGS"
fi
if [[ $FAILED -gt 0 ]]; then
    echo -e "${RED}Failed:${NC} $FAILED"
fi
echo ""

if [[ $FAILED -eq 0 ]]; then
    if [[ $IS_FAT -eq 1 ]]; then
        echo -e "${GREEN}✓ UV setup verified for FAT filesystem${NC}"
        echo ""
        echo "Configuration:"
        echo "  Filesystem: $FS_TYPE (FAT - no symlink support)"
        echo "  UV cache: ${UV_CACHE_DIR:-/tmp/uv-cache}"
        echo "  Virtual env: ${VENV_DIR:-/tmp/metainformant_venv}"
    else
        echo -e "${GREEN}✓ UV setup verified${NC}"
    fi
    exit 0
else
    echo -e "${RED}✗ Some checks failed${NC}"
    echo ""
    echo "Recommendations:"
    if [[ $IS_FAT -eq 1 ]]; then
        echo "  1. Run: bash scripts/package/setup_uv.sh"
        echo "  2. Ensure UV_CACHE_DIR=/tmp/uv-cache is set"
        echo "  3. Create venv in /tmp: uv venv /tmp/metainformant_venv"
    else
        echo "  1. Run: bash scripts/package/setup_uv.sh"
        echo "  2. Check virtual environment setup"
    fi
    exit 1
fi

