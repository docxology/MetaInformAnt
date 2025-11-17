# Testing Setup Guide for METAINFORMANT

This guide documents the complete setup process for running METAINFORMANT tests, particularly when working with external drives.

## Quick Summary

**Problem**: External drive formatted as exFAT doesn't support symlinks → Python venv creation fails  
**Solution**: Either reformat to ext4/NTFS or use `/tmp` for venv  
**Status**: ✅ Working setup established in `/tmp/metainformant-venv`

## Current Working Setup

As of the latest session, we have a fully functional test environment:

```bash
# Virtual environment location
/tmp/metainformant-venv

# Dependencies installed
- pytest 9.0.1
- pytest-cov 7.0.0
- metainformant 0.2.0 (editable install)
- All core dependencies (numpy, pandas, matplotlib, biopython, etc.)
- Scientific dependencies (scipy, scikit-learn, seaborn, networkx)
- Database dependencies (psycopg2-binary)
```

## Setup Commands

### Initial Setup (One-time)

```bash
# Navigate to repository
cd /media/q/ext6/github/MetaInformAnt

# Create venv in /tmp where symlinks work
python3 -m venv /tmp/metainformant-venv

# Install pytest and coverage
/tmp/metainformant-venv/bin/pip install pytest pytest-cov

# Install metainformant in editable mode
/tmp/metainformant-venv/bin/pip install -e .

# Install optional dependencies
/tmp/metainformant-venv/bin/pip install scipy scikit-learn seaborn networkx psycopg2-binary
```

### Running Tests

```bash
# Run all tests
/tmp/metainformant-venv/bin/pytest

# Run with coverage
/tmp/metainformant-venv/bin/pytest --cov=src/metainformant

# Run specific test file
/tmp/metainformant-venv/bin/pytest tests/test_rna_amalgkit.py -v

# Run tests matching a pattern
/tmp/metainformant-venv/bin/pytest -k "test_workflow" -v
```

## Filesystem Issues

### The Problem

Your external drive is mounted as:
```
/dev/sdb2 on /media/q/ext6 type exfat
```

**exFAT does NOT support symlinks**, which are required for:
- Python virtual environments (venv)
- UV package manager cache
- Many Python build tools

### The Solutions

#### Option 1: Reformat to ext4 (Recommended for Linux-only)

```bash
# ⚠️ BACKUP ALL DATA FIRST! ⚠️
sudo umount /dev/sdb2
sudo mkfs.ext4 -L "MetaInformAnt" /dev/sdb2
sudo mount /dev/sdb2 /media/q/ext6
```

**Pros**: Native Linux, full symlink support, best performance  
**Cons**: Not readable on Windows/macOS without third-party tools

#### Option 2: Reformat to NTFS (Cross-platform)

```bash
# ⚠️ BACKUP ALL DATA FIRST! ⚠️
sudo apt-get install ntfs-3g
sudo umount /dev/sdb2
sudo mkfs.ntfs -f -L "MetaInformAnt" /dev/sdb2
sudo mount -t ntfs-3g /dev/sdb2 /media/q/ext6
```

**Pros**: Works on Windows, Linux, macOS (read-only)  
**Cons**: Slightly slower than ext4 on Linux

#### Option 3: Use /tmp (No reformatting - Current Solution)

```bash
# Create venv in /tmp instead of .venv
python3 -m venv /tmp/metainformant-venv
```

**Pros**: No reformatting needed, works immediately  
**Cons**: Venv deleted on reboot, must recreate after restart

## UV Package Manager

The project standardizes on **UV** for dependency management, but UV also requires symlink support.

### UV with exFAT Workaround

```bash
# Set UV cache to /tmp
export UV_CACHE_DIR=/tmp/uv-cache

# Install packages with UV
uv pip install <package> --python /tmp/metainformant-venv/bin/python3
```

### UV with Reformatted Drive (ext4/NTFS)

```bash
# After reformatting to ext4/NTFS, UV works normally
cd /media/q/ext6/github/MetaInformAnt

# Create venv in repo
uv venv

# Install dependencies
uv pip install -e .
```

## Test Configuration

The project's `pyproject.toml` includes pytest configuration:

```toml
[tool.pytest.ini_options]
addopts = [
    "--cov=src/metainformant",
    "--cov-report=term-missing",
    "--cov-report=html",
    "--cov-report=xml",
    "--cov-fail-under=85",
]
```

This requires `pytest-cov` to be installed (which we did).

## Testing Best Practices

### 1. Always Use UV

The project standardizes on UV. When reformatted to ext4/NTFS:

```bash
# Good
uv pip install pytest

# Avoid
pip install pytest
```

### 2. Run Tests Regularly

```bash
# Quick check
/tmp/metainformant-venv/bin/pytest --tb=short -q

# Full check with coverage
/tmp/metainformant-venv/bin/pytest --cov=src/metainformant
```

### 3. Test Markers

```bash
# Skip slow tests
pytest -m "not slow"

# Skip network tests (when offline)
pytest -m "not network"

# Skip tests requiring external tools
pytest -m "not external_tool"
```

## Troubleshooting

### Problem: "Operation not permitted" when creating venv

**Cause**: exFAT doesn't support symlinks  
**Solution**: Use `/tmp` or reformat drive

### Problem: pytest-cov not found

**Cause**: Plugin not installed  
**Solution**: 
```bash
/tmp/metainformant-venv/bin/pip install pytest-cov
```

### Problem: Tests fail with import errors

**Cause**: Optional dependencies not installed  
**Solution**:
```bash
/tmp/metainformant-venv/bin/pip install scipy scikit-learn seaborn networkx
```

### Problem: Cannot write system packages

**Cause**: Externally-managed environment  
**Solution**: Always use venv, never system Python

## Documentation Links

- [Quick Start Guide](../QUICKSTART.md) - Initial setup options
- [External Drive Setup](rna/EXTERNAL_DRIVE_SETUP.md) - Detailed filesystem guide
- [Testing Documentation](testing.md) - Test philosophy and organization

## Summary

✅ **Current Setup**: Fully working test environment in `/tmp/metainformant-venv`  
⚠️ **Limitation**: Venv deleted on reboot (exFAT limitation)  
🔧 **Permanent Fix**: Reformat to ext4 (Linux) or NTFS (cross-platform)  
📝 **Alternative**: Use UV with `UV_CACHE_DIR=/tmp/uv-cache`  

---

*Last Updated*: January 2025  
*Tested On*: Parrot OS with exFAT external drive

