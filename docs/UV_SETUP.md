# UV Setup Guide for METAINFORMANT

Complete guide for setting up and using `uv` (fast Python package manager) with METAINFORMANT, including support for FAT filesystems and external drives.

## Overview

METAINFORMANT uses `uv` for fast, reliable dependency management. This guide covers:

- Installation and basic setup
- FAT filesystem support (exFAT, FAT32)
- Cache directory configuration
- Virtual environment handling
- Testing and verification
- Troubleshooting

## Quick Start

### Standard Setup (ext4, NTFS, etc.)

```bash
# Install uv if not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Run automated setup
bash scripts/package/setup_uv.sh

# Verify setup
bash scripts/package/verify_uv_setup.sh
```

### FAT Filesystem Setup (exFAT, FAT32)

On FAT filesystems, setup is **automatic** - the scripts detect the filesystem and configure appropriately:

```bash
# Same command works on FAT filesystems
bash scripts/package/setup_uv.sh

# The script automatically:
# - Detects FAT filesystem
# - Sets UV_CACHE_DIR=/tmp/uv-cache
# - Creates venv at /tmp/metainformant_venv
```

## Installation

### Installing UV

**macOS (Homebrew)**:
```bash
brew install uv
```

**Official Installer**:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Verify Installation**:
```bash
uv --version
```

## Filesystem Detection

METAINFORMANT automatically detects filesystem type and configures UV accordingly.

### Supported Filesystems

| Filesystem | Symlinks | UV Cache | Venv Location | Notes |
|------------|----------|----------|---------------|-------|
| **ext4** | ✅ Yes | `.uv-cache` | `.venv` | Standard Linux |
| **NTFS** | ✅ Yes | `.uv-cache` | `.venv` | Cross-platform |
| **exFAT** | ❌ No | `/tmp/uv-cache` | `/tmp/metainformant_venv` | Auto-detected |
| **FAT32** | ❌ No | `/tmp/uv-cache` | `/tmp/metainformant_venv` | Auto-detected |
| **Btrfs** | ✅ Yes | `.uv-cache` | `.venv` | Advanced Linux |
| **APFS** | ✅ Yes | `.uv-cache` | `.venv` | macOS |

### Automatic Detection

All setup scripts automatically detect filesystem type:

```bash
# Filesystem detection happens automatically
bash scripts/package/setup_uv.sh

# Output shows detection:
# [1.5/5] Detecting filesystem type and configuring UV cache...
#   → FAT filesystem detected - using /tmp/uv-cache for UV cache
#   → Virtual environment will be created at /tmp/metainformant_venv
```

## Cache Directory Configuration

### Standard Filesystems

On filesystems with symlink support (ext4, NTFS, etc.), UV uses `.uv-cache` in the repository root:

```bash
# Default behavior
cache-dir = ".uv-cache"  # In pyproject.toml
```

### FAT Filesystems

On FAT filesystems (exFAT, FAT32), UV cache **must** be in a location that supports symlinks:

```bash
# Automatically set by setup scripts
export UV_CACHE_DIR="/tmp/uv-cache"
```

**Why?** UV requires symlink support for its cache directory. FAT filesystems don't support symlinks, so the cache must be on a filesystem that does (typically `/tmp`).

### Manual Override

You can manually set the cache directory:

```bash
# For FAT filesystems
export UV_CACHE_DIR="/tmp/uv-cache"

# For standard filesystems (optional)
export UV_CACHE_DIR="/path/to/custom/cache"
```

## Virtual Environment Handling

### Standard Filesystems

Virtual environments are created in the repository root:

```bash
# Standard location
.venv/
```

### FAT Filesystems

On FAT filesystems, virtual environments are created in `/tmp`:

```bash
# FAT filesystem location
/tmp/metainformant_venv/
```

**Why?** Python virtual environments require symlink support, which FAT filesystems don't provide.

**Note**: The `/tmp` venv persists until reboot. After reboot, recreate it:

```bash
uv venv /tmp/metainformant_venv
uv pip install -e . --python /tmp/metainformant_venv/bin/python3
```

### Venv Location Detection

Scripts automatically detect and use the correct venv location:

```bash
# Test scripts automatically find venv
bash scripts/package/uv_test.sh

# Output shows detection:
#   → Filesystem: exfat
#   → UV cache: /tmp/uv-cache
#   → Virtual env: /tmp/metainformant_venv
```

## Setup Scripts

### Main Setup Script

**`scripts/package/setup_uv.sh`** - Complete environment setup:

```bash
# Basic setup (installs dev + scientific deps + amalgkit by default)
bash scripts/package/setup_uv.sh

# Skip amalgkit installation
bash scripts/package/setup_uv.sh --skip-amalgkit

# With all optional dependencies (database, networks, etc.)
bash scripts/package/setup_uv.sh --with-all

# With external CLI dependencies (seqkit, sra-tools, etc.)
bash scripts/package/setup_uv.sh --with-deps

# Skip tests during setup
bash scripts/package/setup_uv.sh --skip-tests
```

**Features**:
- ✅ Automatic filesystem detection
- ✅ Cache directory configuration
- ✅ Virtual environment creation
- ✅ Dependency installation
- ✅ FAT filesystem support

### Development Setup Script

**`scripts/package/uv_dev_setup.sh`** - Development environment:

```bash
bash scripts/package/uv_dev_setup.sh
```

**Features**:
- ✅ Filesystem detection
- ✅ Cache configuration
- ✅ Dependency syncing
- ✅ Pre-commit hooks
- ✅ Task runner scripts

## Testing

### Test Scripts

All test scripts handle FAT filesystems automatically:

```bash
# Standard test runner
bash scripts/package/uv_test.sh

# Optimized test runner
bash scripts/package/uv_test_optimized.sh fast

# Test modes
bash scripts/package/uv_test.sh coverage
bash scripts/package/uv_test.sh parallel
bash scripts/package/uv_test.sh integration
```

**Automatic Handling**:
- ✅ Filesystem detection
- ✅ Cache directory configuration
- ✅ Venv location detection
- ✅ Proper pytest command selection

### Verification Script

**`scripts/package/verify_uv_setup.sh`** - Verify setup:

```bash
bash scripts/package/verify_uv_setup.sh
```

**Checks**:
- ✅ UV installation
- ✅ Filesystem type detection
- ✅ Cache directory configuration
- ✅ Virtual environment location
- ✅ Venv creation capability
- ✅ Basic UV commands
- ✅ Package installation

**Example Output**:
```
🔍 Verifying UV setup for METAINFORMANT...

[1/7] Checking UV installation...
✓ UV is installed: uv 0.1.0

[2/7] Detecting filesystem type...
⚠ FAT filesystem detected: exfat (no symlink support)

[3/7] Checking UV cache directory configuration...
✓ UV_CACHE_DIR is set: /tmp/uv-cache
✓ Cache directory is writable: /tmp/uv-cache

[4/7] Checking virtual environment location...
✓ Virtual environment found at: /tmp/metainformant_venv (Python 3.11.5)

[5/7] Testing venv creation capability...
✓ Venv already exists, skipping creation test

[6/7] Testing basic UV commands...
✓ uv pip list command works
✓ Python in venv is executable

[7/7] Testing package installation capability...
✓ Package installation test passed

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Verification Summary:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Passed: 7
Warnings: 1

✓ UV setup verified for FAT filesystem

Configuration:
  Filesystem: exfat (FAT - no symlink support)
  UV cache: /tmp/uv-cache
  Virtual env: /tmp/metainformant_venv
```

## Usage Examples

### Running Commands

**Standard Filesystem**:
```bash
# Use uv run
uv run pytest
uv run python -m metainformant --help

# Or activate venv
source .venv/bin/activate
pytest
```

**FAT Filesystem**:
```bash
# Use venv Python directly
/tmp/metainformant_venv/bin/python3 -m pytest
/tmp/metainformant_venv/bin/python3 -m metainformant --help

# Or activate venv
source /tmp/metainformant_venv/bin/activate
pytest
```

### Installing Packages

**Standard Filesystem**:
```bash
uv pip install <package>
```

**FAT Filesystem**:
```bash
# Cache directory automatically set by scripts
export UV_CACHE_DIR="/tmp/uv-cache"
uv pip install <package> --python /tmp/metainformant_venv/bin/python3
```

### Running Tests

**All Filesystems** (automatic detection):
```bash
# Test scripts handle everything automatically
bash scripts/package/uv_test.sh
bash scripts/package/uv_test_optimized.sh fast
```

## Troubleshooting

### Problem: "Operation not permitted" when creating venv

**Cause**: FAT filesystem doesn't support symlinks

**Solution**: Scripts automatically handle this. If manual setup needed:

```bash
export UV_CACHE_DIR="/tmp/uv-cache"
uv venv /tmp/metainformant_venv
```

### Problem: "failed to symlink file" errors

**Cause**: UV cache directory on FAT filesystem

**Solution**: Set cache directory to `/tmp`:

```bash
export UV_CACHE_DIR="/tmp/uv-cache"
```

### Problem: Virtual environment not found

**Cause**: Venv location differs on FAT filesystems

**Solution**: Check both locations:

```bash
# Check standard location
ls -la .venv/bin/python3

# Check FAT filesystem location
ls -la /tmp/metainformant_venv/bin/python3

# Or run verification script
bash scripts/package/verify_uv_setup.sh
```

### Problem: Packages not installing

**Cause**: Cache directory not configured for FAT filesystem

**Solution**: Ensure UV_CACHE_DIR is set:

```bash
export UV_CACHE_DIR="/tmp/uv-cache"
uv pip install <package> --python /path/to/venv/bin/python3
```

### Problem: Tests fail with import errors

**Cause**: Venv not activated or wrong Python used

**Solution**: Use correct Python:

```bash
# FAT filesystem
/tmp/metainformant_venv/bin/python3 -m pytest

# Standard filesystem
uv run pytest
```

### Problem: Venv deleted after reboot

**Cause**: `/tmp` is cleared on reboot (FAT filesystem limitation)

**Solution**: Recreate venv:

```bash
export UV_CACHE_DIR="/tmp/uv-cache"
uv venv /tmp/metainformant_venv
uv pip install -e . --python /tmp/metainformant_venv/bin/python3
```

Or use setup script:

```bash
bash scripts/package/setup_uv.sh
```

## Filesystem Detection API

METAINFORMANT provides Python utilities for filesystem detection:

```python
from metainformant.core.filesystem import (
    detect_filesystem_type,
    supports_symlinks,
    get_uv_cache_dir,
    get_venv_location,
    is_fat_filesystem,
)

# Detect filesystem type
fs_type = detect_filesystem_type(".")
print(f"Filesystem: {fs_type}")

# Check symlink support
has_symlinks = supports_symlinks(".")
print(f"Symlinks supported: {has_symlinks}")

# Get appropriate UV cache directory
cache_dir = get_uv_cache_dir(".")
print(f"UV cache: {cache_dir}")

# Get appropriate venv location
venv_dir = get_venv_location(".")
print(f"Venv location: {venv_dir}")

# Check if FAT filesystem
is_fat = is_fat_filesystem(".")
print(f"Is FAT: {is_fat}")
```

## Best Practices

### 1. Always Use Setup Scripts

Setup scripts handle all filesystem detection automatically:

```bash
# Recommended
bash scripts/package/setup_uv.sh
```

### 2. Verify Setup

Always verify setup after installation:

```bash
bash scripts/package/verify_uv_setup.sh
```

### 3. Check Filesystem Type

If experiencing issues, check filesystem type:

```bash
df -T .
```

### 4. Use Verification Script

When troubleshooting, run verification script:

```bash
bash scripts/package/verify_uv_setup.sh
```

### 5. FAT Filesystem Considerations

On FAT filesystems:
- ✅ Cache directory automatically uses `/tmp/uv-cache`
- ✅ Venv automatically uses `/tmp/metainformant_venv`
- ⚠️ Venv deleted on reboot (recreate with setup script)
- ⚠️ No manual configuration needed

## Related Documentation

- **[Setup Guide](setup.md)** - General setup instructions
- **[Testing Guide](testing.md)** - Testing documentation
- **[External Drive Setup](rna/EXTERNAL_DRIVE_SETUP.md)** - Detailed external drive guide
- **[Quick Start](../QUICKSTART.md)** - Quick start guide

## Summary

METAINFORMANT's UV setup is **fully automatic** on all filesystems:

✅ **Automatic filesystem detection** - No manual configuration needed  
✅ **Automatic cache directory** - Uses `/tmp/uv-cache` on FAT filesystems  
✅ **Automatic venv location** - Uses `/tmp/metainformant_venv` on FAT filesystems  
✅ **All scripts handle FAT** - Setup, testing, and verification scripts work automatically  
✅ **Comprehensive verification** - Verification script confirms proper setup  

**No manual `UV_CACHE_DIR` exports needed** - everything is handled automatically!

---

*Last Updated: January 2025*  
*Tested on: ext4, NTFS, exFAT, FAT32 filesystems*

