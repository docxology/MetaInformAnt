# RNA-seq Workflow Setup Guide

## System Requirements

This guide covers setup for Linux systems (Debian/Ubuntu-based). The workflow requires:
- Python 3.11+
- uv (Python package manager)
- SRA Toolkit (for downloading sequencing data)
- kallisto or salmon (for RNA-seq quantification)
- amalgkit (RNA-seq metadata and workflow management)

## Quick Start

```bash
cd /home/q/Documents/GitHub/MetaInformAnt

# 1. Create virtual environment and install Python packages
uv venv
source .venv/bin/activate
uv pip install -e .
uv pip install git+https://github.com/kfuku52/amalgkit

# 2. Install SRA Toolkit (for downloading FASTQ files)
sudo apt-get update
sudo apt-get install -y sra-toolkit

# 3. Install kallisto (for RNA-seq quantification)
sudo apt-get install -y kallisto

# 4. Install fastp and seqkit (required for amalgkit getfastq)
sudo apt-get install -y fastp
sudo apt-get install -y seqkit

# 5. Verify installations
amalgkit --version
fasterq-dump --version
kallisto version
fastp --version
seqkit version

# 6. Run multi-species workflow (auto-activates venv if .venv exists)
# Prerequisites: .venv must exist with amalgkit installed (see steps 1-5 above)
python3 scripts/rna/run_multi_species.py

# With configurable threads:
export AK_THREADS=12
python3 scripts/rna/run_multi_species.py
```

**Note**: The script automatically activates the virtual environment if `.venv` exists. The virtual environment must be set up first (steps 1-5 above).

## Detailed Setup Instructions

### 1. Python Environment Setup

The project uses `uv` for fast, reliable Python package management:

```bash
# Install uv if not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create virtual environment
cd /home/q/Documents/GitHub/MetaInformAnt
uv venv

# Activate virtual environment
source .venv/bin/activate

# Install metainformant in editable mode
uv pip install -e .

# Install amalgkit
uv pip install git+https://github.com/kfuku52/amalgkit
```

**Verification:**
```bash
python3 -c "import metainformant; print('metainformant OK')"
amalgkit --version  # Should show: amalgkit version 0.12.19 or higher
```

### 2. SRA Toolkit Installation

Required for downloading RNA-seq data from NCBI SRA:

```bash
# Debian/Ubuntu
sudo apt-get update
sudo apt-get install -y sra-toolkit

# Verify installation
fasterq-dump --version
prefetch --version
```

**Important Configuration:**
The multi-species script automatically configures SRA Toolkit for large downloads:
- Creates wrapper script: `output/sra_temp/fasterq-dump` with `--size-check off`
- Sets temp directory to project location instead of `/tmp`
- Configures environment variables: `TMPDIR`, `TEMP`, `TMP`
- No manual SRA configuration needed
vdb-config --version

# Configure SRA Toolkit (first time only)
vdb-config --interactive  # Set cache location and preferences
```

**Alternative: Manual Installation**
```bash
# Download latest version
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
export PATH=$PATH:$(pwd)/sratoolkit.*/bin
```

### 3. RNA-seq Quantification Tools

Install at least one of these tools for transcript quantification:

#### Option A: kallisto (Recommended - Fast and Accurate)
```bash
# Debian/Ubuntu
sudo apt-get install -y kallisto

# Verify
kallisto version  # Should show: kallisto, version 0.48.0 or higher
```

#### Option B: salmon (Alternative)
```bash
# Download latest release
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
tar -xzf salmon-1.10.0_linux_x86_64.tar.gz
export PATH=$PATH:$(pwd)/salmon-latest_linux_x86_64/bin

# Verify
salmon --version
```

### 4. fastp and seqkit Installation

Required by amalgkit for FASTQ quality control and sequence processing:

```bash
# Debian/Ubuntu
sudo apt-get update
sudo apt-get install -y fastp seqkit

# Verify installations
fastp --version  # Should show 0.23.2 or higher
seqkit version   # Should show 2.3.0 or higher
```

**Note**: These tools are required for the `getfastq` step. The script automatically checks for them.

### 5. System Dependencies

```bash
# Install basic build tools and libraries
sudo apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    python3-dev \
    libz-dev \
    libbz2-dev \
    liblzma-dev
```

## Environment Configuration

### Required Environment Variables

```bash
# Add to ~/.bashrc or ~/.zshrc
export NCBI_EMAIL="your.email@example.com"  # Required for NCBI API

# Optional: Configure data directories
export METAINFORMANT_OUTPUT="/home/q/Documents/GitHub/MetaInformAnt/output"
export METAINFORMANT_DATA="/home/q/Documents/GitHub/MetaInformAnt/data"
```

### SRA Toolkit Configuration

Configure caching and download preferences:

```bash
# Run interactive configuration
vdb-config --interactive

# Or configure via command line
vdb-config --set /repository/user/main/public/root=output/sra_cache
```

## Verification Checklist

Run these commands to verify complete setup:

```bash
cd /home/q/Documents/GitHub/MetaInformAnt
source .venv/bin/activate

# Check Python packages
python3 -c "import metainformant; print('✓ metainformant')"
python3 -c "import amalgkit; print('✓ amalgkit')"

# Check command-line tools
command -v amalgkit && echo "✓ amalgkit CLI"
command -v fasterq-dump && echo "✓ fasterq-dump"
command -v prefetch && echo "✓ prefetch"
command -v kallisto && echo "✓ kallisto"

# Check environment
[ -n "$NCBI_EMAIL" ] && echo "✓ NCBI_EMAIL set" || echo "⚠ NCBI_EMAIL not set"

# Run quick test
python3 -c "from metainformant.rna.workflow import load_workflow_config; print('✓ RNA workflow module')"
```

## Platform-Specific Notes

### Linux (Current System)
- Uses system package manager (apt) for most tools
- Native performance, recommended for production
- All tools have precompiled binaries

### macOS (Previous System)
- Use Homebrew for package management:
  ```bash
  brew install sratoolkit kallisto
  ```
- May need Xcode Command Line Tools

### Migration from macOS to Linux
If migrating from macOS:
1. Remove old `.venv` and recreate: `rm -rf .venv && uv venv`
2. Reinstall system tools using apt instead of brew
3. Verify all paths in configs use absolute paths (no `/Users/...`)
4. Update any hardcoded Mac-specific paths

## Troubleshooting

### "amalgkit: command not found"
```bash
# Ensure virtual environment is activated
source .venv/bin/activate

# Reinstall if needed
uv pip install --force-reinstall git+https://github.com/kfuku52/amalgkit
```

### "fasterq-dump: command not found"
```bash
# Install SRA Toolkit
sudo apt-get install -y sra-toolkit

# Or add to PATH if manually installed
export PATH=$PATH:/path/to/sratoolkit/bin
```

### "Permission denied" errors
```bash
# Ensure output directory is writable
chmod -R u+w output/

# Run script from repository root
# Prerequisites: .venv must exist with amalgkit installed
# If not set up, run:
#   python3 -m venv .venv
#   source .venv/bin/activate
#   pip install -e .
#   pip install git+https://github.com/kfuku52/amalgkit

cd /home/q/Documents/GitHub/MetaInformAnt
python3 scripts/rna/run_multi_species.py

# With configurable threads:
export AK_THREADS=12
python3 scripts/rna/run_multi_species.py
```

### PEP 668 "externally-managed-environment" error
This is expected on modern Linux systems. Always use virtual environments:
```bash
source .venv/bin/activate  # Always activate before installing packages
```

## Running the Workflow

### Multi-Species RNA-seq Analysis

The script provides automatic environment management:

```bash
cd /home/q/Documents/GitHub/MetaInformAnt

# No manual venv activation needed - script handles it automatically (if .venv exists)
# Prerequisites: .venv must exist with amalgkit installed
# If not set up, run:
#   python3 -m venv .venv
#   source .venv/bin/activate
#   pip install -e .
#   pip install git+https://github.com/kfuku52/amalgkit

python3 scripts/rna/run_multi_species.py

# With configurable threads:
export AK_THREADS=12
python3 scripts/rna/run_multi_species.py

# Output will be in:
# - output/amalgkit/<species>/work/     - Intermediate files
# - output/amalgkit/<species>/logs/     - Log files  
# - output/amalgkit/<species>/quant/    - Quantification results
# - output/amalgkit/<species>/merged/   - Merged expression matrices
```

**Features:**
- Auto-activation: Detects and activates virtual environment automatically
- Auto-discovery: Finds all species configs in `config/amalgkit/`
- Batched processing: 10 samples at a time (~20-50 GB peak disk usage)
- SRA optimization: Automatic wrapper and environment setup
- Resume capability: Skips already-quantified samples

### Single Species Analysis

```bash
# Example: Pogonomyrmex barbatus
python3 -m metainformant.rna.workflow \
    --config config/amalgkit/amalgkit_pbarbatus.yaml
```

## Performance Optimization

### Parallel Processing
```bash
# Adjust threads in config files
# config/amalgkit/amalgkit_*.yaml
threads: 10  # Default: 10 threads for downloads and quantification
```

### Disk Space Management
The script automatically manages disk space:
```yaml
# Batched processing configuration (automatic)
# - Downloads 10 samples at a time
# - Quantifies batch
# - Deletes FASTQs after quantification
# - Peak usage: ~20-50 GB (10 samples of FASTQs)

# SRA temp directory (automatic)
# - Uses output/sra_temp instead of /tmp
# - Wrapper script disables conservative size checks
# - Symlinks for fastp, kallisto, seqkit
```

No manual configuration needed - script handles all disk management automatically.

### Download Acceleration
```bash
# Use cloud mirrors for faster downloads
# In config YAML:
steps:
  getfastq:
    aws: yes  # Use AWS Open Data Program
    gcp: yes  # Use Google Cloud
    ncbi: yes # Use NCBI directly
```

## Getting Help

- **amalgkit Issues**: https://github.com/kfuku52/amalgkit/issues
- **SRA Toolkit**: https://github.com/ncbi/sra-tools/wiki
- **kallisto**: https://pachterlab.github.io/kallisto/

## References

- amalgkit: https://github.com/kfuku52/amalgkit
- SRA Toolkit: https://github.com/ncbi/sra-tools
- kallisto: https://pachterlab.github.io/kallisto/
- uv: https://github.com/astral-sh/uv

