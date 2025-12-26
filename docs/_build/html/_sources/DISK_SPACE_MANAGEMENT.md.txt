# Disk Space Management Guide

Critical guide for managing temporary files and disk space when running METAINFORMANT, especially on external drives.

## The Problem

On external drives or systems with limited `/tmp` space, the default system temporary directory may:
- Be a small RAM-based filesystem (tmpfs)
- Fill up quickly during large data processing
- Cause "No space left on device" errors
- Interrupt long-running workflows

## Solution: Repository-Local Temp Directories

METAINFORMANT uses repository-local temp directories instead of system `/tmp`:

```
metainformant/
├── .tmp/              # Repository-local temp directory
│   ├── python/        # Python tempfile operations
│   └── bash/          # Shell script heredocs and temp files
├── .cache/            # Cache directory for expensive operations
└── output/            # All program outputs
```

## Configuration

### Environment Variables

Set these before running any METAINFORMANT commands:

```bash
# Add to ~/.bashrc or ~/.zshrc for persistence
export TMPDIR="$(pwd)/.tmp/bash"
export TEMP=$TMPDIR
export TMP=$TMPDIR

# Create directories
mkdir -p "$TMPDIR"

# Verify configuration
echo $TMPDIR
# Should show: /path/to/metainformant/.tmp/bash
```

### Python Configuration

For Python scripts, configure at the start of your script:

```python
import os
import tempfile
from pathlib import Path

# Set temp directory to repository root
REPO_ROOT = Path(__file__).resolve().parent.parent
TEMP_DIR = REPO_ROOT / ".tmp" / "python"
TEMP_DIR.mkdir(parents=True, exist_ok=True)

# Override Python temp directory
os.environ["TMPDIR"] = str(TEMP_DIR)
os.environ["TEMP"] = str(TEMP_DIR)
os.environ["TMP"] = str(TEMP_DIR)
tempfile.tempdir = str(TEMP_DIR)

# Now all temp operations use the repository temp directory
with tempfile.NamedTemporaryFile() as tmp:
    # tmp.name is in .tmp/python/
    pass
```

### Shell Script Pattern

For shell scripts, add this at the start:

```bash
#!/bin/bash
# At start of every shell script:

# Get repository root from script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Set temp directory to external drive
export TMPDIR="$REPO_ROOT/.tmp/bash"
mkdir -p "$TMPDIR"

# Now all bash operations use external drive temp
```

## Git Configuration

Ensure `.tmp/` and `.cache/` are ignored by git:

```gitignore
# Temporary files (repository-local)
.tmp/
.cache/
*.tmp
*.temp
```

These should already be in the repository `.gitignore`.

## Troubleshooting

### Error: "No space left on device"

**Cause**: System `/tmp` is full or too small.

**Solution**:
1. Set `TMPDIR` environment variable as described above
2. Verify: `echo $TMPDIR` should show repository path
3. Clear old temp files: `rm -rf .tmp/* .cache/*`
4. Restart your workflow

### Error: Temp files in wrong location

**Cause**: Environment variables not set correctly.

**Diagnosis**:
```bash
# Check current temp directory
python3 -c "import tempfile; print(tempfile.gettempdir())"

# Should show: /path/to/metainformant/.tmp/python
# NOT: /tmp or /var/tmp
```

**Solution**: Ensure environment variables are set before running Python.

### Error: Permission denied on temp directory

**Cause**: Incorrect permissions on `.tmp/` directory.

**Solution**:
```bash
chmod -R 755 .tmp/
chmod -R 755 .cache/
```

### Large temporary files accumulating

**Cause**: Failed workflows leave behind temp files.

**Solution**:
```bash
# Check temp directory size
du -sh .tmp/

# Clean up safely (preserves directory structure)
find .tmp/ -type f -mtime +1 -delete
find .cache/ -type f -mtime +7 -delete
```

## Disk Space Monitoring

### Check available space

```bash
# On external drive
df -h /path/to/metainformant

# Check temp directory usage
du -sh .tmp/ .cache/ output/
```

### Automated cleanup script

Add to `scripts/core/cleanup_temp.sh`:

```bash
#!/bin/bash
# Clean temporary files older than 1 day

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "Cleaning temp files in $REPO_ROOT..."

# Clean Python temp (files older than 1 day)
find "$REPO_ROOT/.tmp/python" -type f -mtime +1 -delete 2>/dev/null

# Clean bash temp (files older than 1 day)
find "$REPO_ROOT/.tmp/bash" -type f -mtime +1 -delete 2>/dev/null

# Clean cache (files older than 7 days)
find "$REPO_ROOT/.cache" -type f -mtime +7 -delete 2>/dev/null

# Report remaining usage
echo "Remaining usage:"
du -sh "$REPO_ROOT/.tmp" "$REPO_ROOT/.cache" 2>/dev/null
```

## External Drive Setup

For best results on external drives:

### FAT/exFAT Filesystems

FAT filesystems don't support symlinks, which affects `uv` virtual environments.

**Solution**: Use the automated setup script:
```bash
bash scripts/package/setup.sh
```

This automatically:
- Detects FAT filesystem
- Sets `UV_CACHE_DIR=/tmp/uv-cache`
- Creates venv at `/tmp/metainformant_venv`

**Note**: The `/tmp` venv is deleted on reboot. Re-run setup after restart.

### Better Alternative: Reformat Drive

For best compatibility, reformat your external drive to:
- **ext4** (Linux-only, best performance)
- **NTFS** (cross-platform, Windows-compatible)

These filesystems support symlinks natively.

## Best Practices

1. **Always set TMPDIR** before running long workflows
2. **Monitor disk space** during large data processing
3. **Clean temp files regularly** to prevent accumulation
4. **Use output/ for results**, not temp directories
5. **Never store important data** in `.tmp/` or `.cache/`
6. **Reformat FAT drives** to ext4/NTFS if possible

## Related Documentation

- **[UV Setup Guide](UV_SETUP.md)** - Package management and virtual environment setup
- **[External Drive Setup](rna/EXTERNAL_DRIVE_SETUP.md)** - RNA workflow setup on external drives
- **[Architecture](architecture.md)** - Directory structure and conventions

---

*Proper temp directory management ensures reliable execution of METAINFORMANT workflows on any system configuration.*





