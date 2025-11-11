# External Drive Setup and Filesystem Considerations

This document covers setup considerations, known issues, and solutions when running METAINFORMANT RNA workflows on external drives or filesystems with limitations (e.g., ext6 without symlink support).

## Overview

METAINFORMANT is designed to work seamlessly on external drives, but some filesystems (particularly ext6) have limitations that require special handling. This guide documents all known issues and their solutions.

## Filesystem Limitations

### Symlink Support

**Issue**: Some filesystems (e.g., ext6) don't support symlinks in certain locations, which can cause failures when:
- Creating virtual environments
- Using `uv` package manager cache
- Building Python packages

**Symptoms**:
```
Operation not permitted (os error 1)
failed to symlink file from ... to ...
```

**Solution**: The codebase automatically handles this by:
1. Using alternative venv location (`/tmp/metainformant_venv`) when repo `.venv` fails
2. Using `/tmp/uv-cache` for uv cache directory instead of `output/.uv-cache`
3. Detecting filesystem capabilities and falling back gracefully

### Virtual Environment Location

**Automatic Detection**: The setup utilities (`scripts/rna/_setup_utils.py`) automatically:
1. Try to create `.venv` in repo root first
2. If that fails due to symlink issues, use `/tmp/metainformant_venv`
3. Log which location is being used

**Manual Override**: You can set the venv location via environment variable:
```bash
export METAINFORMANT_VENV=/path/to/your/venv
```

## UV Package Manager Configuration

### Cache Directory

**Issue**: `uv` requires symlink support for its cache directory. On ext6 filesystems, using `output/.uv-cache` fails.

**Solution**: The setup utilities automatically use `/tmp/uv-cache` for uv cache:
- Set via `UV_CACHE_DIR=/tmp/uv-cache` environment variable
- Automatically configured in `_setup_utils.py`

**Manual Override**:
```bash
export UV_CACHE_DIR=/tmp/uv-cache
uv pip install <package> --python /path/to/venv/bin/python3
```

### Installing Dependencies

When installing packages with `uv`, always use the cache directory that supports symlinks:

```bash
# Correct: Use /tmp/uv-cache
export UV_CACHE_DIR=/tmp/uv-cache
uv pip install biopython --python /tmp/metainformant_venv/bin/python3

# Incorrect: Will fail on ext6
uv pip install biopython --python /tmp/metainformant_venv/bin/python3
# (without UV_CACHE_DIR set, defaults to output/.uv-cache which may fail)
```

## Temporary Directory Selection

### Automatic Selection

The codebase uses `get_recommended_temp_dir()` from `metainformant.core.disk` to:
1. Check if `TMPDIR` environment variable is set (use that)
2. Check if `output/` directory is on a large/medium drive (use `output/.tmp`)
3. Fall back to system temp directory (`/tmp`)

### Manual Override

Set `TMPDIR` environment variable:
```bash
export TMPDIR=/path/to/large/temp/directory
```

## Known Issues and Solutions

### Issue 1: Biopython Not Available Warning

**Symptom**: 
```
WARNING | metainformant.rna.discovery | Biopython not available - discovery functions will be limited
```

**Cause**: Biopython is listed in `pyproject.toml` but wasn't installed in the virtual environment, often due to:
- Venv created before biopython was added to dependencies
- Installation failed due to symlink issues

**Solution**:
```bash
# Install biopython in venv using proper cache location
export UV_CACHE_DIR=/tmp/uv-cache
uv pip install biopython --python /tmp/metainformant_venv/bin/python3

# Or reinstall metainformant with all dependencies
uv pip install -e . --python /tmp/metainformant_venv/bin/python3
```

### Issue 2: NCBI Datasets Library Not Available

**Symptom**:
```
WARNING | metainformant.rna.discovery | ncbi-datasets-pylib not available - genome info will be limited
```

**Solution**:
```bash
export UV_CACHE_DIR=/tmp/uv-cache
uv pip install ncbi-datasets-pylib --python /tmp/metainformant_venv/bin/python3
```

### Issue 3: UV Cache Symlink Errors

**Symptom**:
```
× Failed to download `biopython==1.86`
├─▶ Failed to read from the distribution cache
╰─▶ failed to symlink file from ... to ...
   Operation not permitted (os error 1)
```

**Solution**: Set `UV_CACHE_DIR` to a location that supports symlinks:
```bash
export UV_CACHE_DIR=/tmp/uv-cache
```

### Issue 4: Virtual Environment Creation Fails

**Symptom**:
```
Failed to create venv at .venv: Operation not permitted
```

**Solution**: The codebase automatically falls back to `/tmp/metainformant_venv`. This is handled transparently.

## Best Practices

### 1. Environment Setup

Always use the provided setup utilities which handle filesystem limitations automatically:

```bash
# Run environment check (auto-handles venv setup)
python3 scripts/rna/check_environment.py

# Run workflow (auto-handles all setup)
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_<species>.yaml
```

### 2. Manual Package Installation

When manually installing packages, always set the cache directory:

```bash
export UV_CACHE_DIR=/tmp/uv-cache
uv pip install <package> --python /tmp/metainformant_venv/bin/python3
```

### 3. Temporary Directory

For large workflows, ensure sufficient space:

```bash
# Check available space
df -h /tmp
df -h /media/q/ext6/github/MetaInformAnt/output

# Set TMPDIR if needed
export TMPDIR=/path/to/large/drive/tmp
```

### 4. Disk Space Management

RNA-seq workflows can use significant disk space:
- Genome downloads: 100MB - 5GB per species
- FASTQ files: 1-10GB per sample
- Quantification results: 50-500MB per sample

**Recommendation**: Use `keep_fastq: no` in config to delete FASTQ files after quantification.

## Configuration Examples

### Standard Configuration (Works on Any Filesystem)

```yaml
# config/amalgkit/amalgkit_species.yaml
work_dir: output/amalgkit/species/work
log_dir: output/amalgkit/species/logs
threads: 12

# Automatically handles external drive issues
auto_install_amalgkit: true

steps:
  getfastq:
    keep_fastq: no  # Delete FASTQ after quantification to save space
  quant:
    keep_fastq: no
```

### Large Dataset Configuration

```yaml
# For workflows with many samples on external drive
work_dir: output/amalgkit/species/work
log_dir: output/amalgkit/species/logs
threads: 24

steps:
  getfastq:
    num_download_workers: 10  # Parallel downloads
    keep_fastq: no  # Critical for disk space
  quant:
    keep_fastq: no
    threads: 24
```

## Verification

### Check Environment

```bash
# Verify all dependencies including biopython
python3 scripts/rna/check_environment.py
```

### Verify Virtual Environment

```bash
# Check venv location
/tmp/metainformant_venv/bin/python3 -c "import sys; print(sys.prefix)"

# Check biopython
/tmp/metainformant_venv/bin/python3 -c "from Bio import Entrez; print('✓ biopython available')"

# Check ncbi-datasets
/tmp/metainformant_venv/bin/python3 -c "from ncbi.datasets import GenomeApi; print('✓ ncbi-datasets available')"
```

### Verify UV Cache

```bash
# Check cache location
echo $UV_CACHE_DIR
# Should be: /tmp/uv-cache

# Verify cache is writable
ls -la /tmp/uv-cache
```

## Troubleshooting

### Problem: Workflow fails with symlink errors

**Solution**: Ensure `UV_CACHE_DIR` is set:
```bash
export UV_CACHE_DIR=/tmp/uv-cache
```

### Problem: Virtual environment not found

**Solution**: The setup utilities automatically create it. Check:
```bash
ls -la /tmp/metainformant_venv/bin/python3
```

### Problem: Out of disk space

**Solution**: 
1. Check space: `df -h`
2. Clean up old outputs: `rm -rf output/amalgkit/old_species/`
3. Use `keep_fastq: no` in config
4. Set `TMPDIR` to a drive with more space

### Problem: Packages not installing

**Solution**: 
1. Set `UV_CACHE_DIR=/tmp/uv-cache`
2. Verify venv Python: `/tmp/metainformant_venv/bin/python3 --version`
3. Try manual install: `uv pip install <package> --python /tmp/metainformant_venv/bin/python3`

## Implementation Details

### Automatic Handling

The following files automatically handle external drive issues:

1. **`scripts/rna/_setup_utils.py`**:
   - Detects filesystem limitations
   - Uses alternative venv location when needed
   - Sets `UV_CACHE_DIR` to `/tmp/uv-cache`
   - Handles symlink failures gracefully

2. **`src/metainformant/core/disk.py`**:
   - `get_recommended_temp_dir()`: Selects best temp directory
   - `detect_drive_size_category()`: Detects drive size for optimization
   - `get_recommended_batch_size()`: Optimizes batch processing

3. **`src/metainformant/rna/discovery.py`**:
   - Gracefully handles missing biopython/ncbi-datasets
   - Logs warnings but continues operation
   - Functions degrade gracefully when optional deps missing

### Code Locations

- Venv location detection: `scripts/rna/_setup_utils.py::_find_venv_location()`
- UV cache configuration: `scripts/rna/_setup_utils.py::setup_venv_and_dependencies()`
- Temp directory selection: `src/metainformant/core/disk.py::get_recommended_temp_dir()`

## Related Documentation

- [RNA Getting Started Guide](GETTING_STARTED.md) - Initial setup instructions
- [RNA Configuration Guide](CONFIGURATION.md) - Configuration options
- [RNA Workflow Guide](workflow.md) - Workflow execution
- [Core Path Management](../core/paths.md) - Path handling utilities

## Summary

METAINFORMANT automatically handles external drive and filesystem limitations:

✅ **Automatic venv location selection** - Falls back to `/tmp/metainformant_venv` when needed  
✅ **Automatic UV cache configuration** - Uses `/tmp/uv-cache` to avoid symlink issues  
✅ **Automatic temp directory selection** - Chooses best location based on drive size  
✅ **Graceful degradation** - Optional dependencies (biopython, ncbi-datasets) are optional  
✅ **Comprehensive error handling** - Clear messages and automatic fallbacks  

No manual configuration needed - the system handles everything automatically!

---

*Last Updated: November 2025*  
*Tested on: ext6 filesystem, external drives, standard Linux filesystems*

