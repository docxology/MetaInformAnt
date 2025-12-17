# Core Utility Scripts

Core infrastructure and utility scripts for METAINFORMANT development and deployment.

## Directory Structure

```
scripts/core/
├── run_demo.py                 # Complete workflow demonstration
├── fix_disk_space.sh           # External drive temp directory setup
├── setup_temp_dirs.sh          # Repository-local temp directory management
├── README.md                   # This file
└── AGENTS.md                   # AI agent contribution documentation
```

## Scripts Overview

### Workflow Demonstration (`run_demo.py`)

Complete workflow demonstration showing best practices for METAINFORMANT usage.

**Features:**
- Configuration loading and validation
- Data processing with statistics
- Visualization generation with informative filenames
- Comprehensive output organization
- Error handling and logging

**Usage:**
```bash
# Run complete workflow demonstration
python3 scripts/core/run_demo.py
```

**Output Structure:**
```
output/demo/
├── workflow_configuration.json    # Workflow parameters and settings
├── input_samples.json             # Generated input data
├── processed_normalized_data.json # Processed results with statistics
├── workflow_summary_report.json   # Summary metrics and status
└── visualizations/
    ├── input_sample_values_lineplot.png
    ├── normalized_values_barplot.png
    ├── sample_correlation_heatmap.png
    └── visualization_metadata.json
```

### Disk Space Management (`fix_disk_space.sh`)

Fixes "out of space" errors by setting up external drive temp directories.

**Purpose:**
- Creates repository-local temp directories (`.tmp/`)
- Sets up separate temp spaces for bash, python, git, and downloads
- Configures proper permissions for temp directories
- Works on external drives and limited filesystems

**Usage:**
```bash
# Setup temp directories for external drive
bash scripts/core/fix_disk_space.sh
```

**What it creates:**
- `.tmp/bash/` - Bash temp directory
- `.tmp/python/` - Python temp directory
- `.tmp/git/` - Git operations temp space
- `.tmp/downloads/` - Download cache directory
- `.cache/` - General cache directory

### Temp Directory Setup (`setup_temp_dirs.sh`)

Generalized temp directory setup for any repository location and filesystem.

**Features:**
- Automatic filesystem detection
- Works on external drives, home directories, or any filesystem
- Creates comprehensive temp directory structure
- Handles filesystem limitations transparently
- Sets appropriate permissions

**Usage:**
```bash
# Setup temp directories for current repository
bash scripts/core/setup_temp_dirs.sh
```

**Filesystem Support:**
- **Regular filesystems**: Uses `.tmp/` in repository root
- **Limited filesystems** (e.g., FAT, some external drives): Uses `/tmp/metainformant_*`
- **Automatic detection**: No manual configuration needed

## Integration

These scripts integrate with:
- **METAINFORMANT core utilities**: `core.io`, `core.paths`, `core.logging`
- **External drive support**: Automatic detection and handling
- **Development workflows**: Setup and demonstration tools
- **Package management**: Temp directory management for builds

## Dependencies

- **Python scripts**: Require `metainformant.core` modules
- **Shell scripts**: Standard bash, no external dependencies
- **Filesystem**: Works on all common filesystems (ext4, NTFS, FAT, etc.)

## Development Notes

### Best Practices Demonstrated
- **Configuration management**: Save all parameters to output directory
- **Informative filenames**: Use descriptive names for output files
- **Comprehensive logging**: Log all operations with context
- **Error handling**: Graceful failure with clear error messages
- **Output organization**: Structured output in `output/` directory

### File System Handling
- **External drives**: Automatic temp directory setup
- **Limited filesystems**: Fallback to `/tmp/` directories
- **Permissions**: Proper temp directory permissions (1777/777)
- **Cleanup**: Scripts handle cleanup automatically

## Related Documentation

- [Core Utilities Documentation](../../docs/core/README.md)
- [Disk Space Management](../../docs/DISK_SPACE_MANAGEMENT.md)
- [UV Setup Guide](../../docs/UV_SETUP.md)
- [METAINFORMANT CLI](../../docs/cli.md)
