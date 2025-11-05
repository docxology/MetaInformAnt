# Getting Started with RNA-seq Analysis

Complete guide for setting up and running your first RNA-seq workflow with METAINFORMANT.

## Quick Start

### 1. Installation

```bash
cd /path/to/MetaInformAnt

# Create virtual environment
uv venv
source .venv/bin/activate

# Install metainformant
uv pip install -e .

# Install amalgkit
uv pip install git+https://github.com/kfuku52/amalgkit

# Install system dependencies
sudo apt-get install -y sra-toolkit kallisto fastp seqkit
```

### 2. Verify Installation

```bash
# Check Python packages
python3 -c "import metainformant; print('✓ metainformant')"
amalgkit --version

# Check command-line tools
command -v amalgkit && echo "✓ amalgkit CLI"
command -v fasterq-dump && echo "✓ fasterq-dump"
command -v kallisto && echo "✓ kallisto"

# Set NCBI email (required)
export NCBI_EMAIL="your.email@example.com"
```

### 3. Run Your First Workflow

```bash
# Single species workflow (ENA-based, recommended)
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 12 \
  --threads 12
```

## System Requirements

- **Python 3.11+** with virtual environment (`.venv/`)
- **uv** (Python package manager) or pip
- **amalgkit** (installed in venv)
- **wget** (for ENA downloads) or **SRA Toolkit** (for SRA downloads)
- **kallisto** or **salmon** (for quantification)
- **fastp** and **seqkit** (for quality control)

## Detailed Setup

### Python Environment Setup

The project uses `uv` for fast, reliable Python package management:

```bash
# Install uv if not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create virtual environment
cd /path/to/MetaInformAnt
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

### System Dependencies

#### SRA Toolkit (for SRA downloads)
```bash
# Debian/Ubuntu
sudo apt-get update
sudo apt-get install -y sra-toolkit

# Verify installation
fasterq-dump --version
prefetch --version
```

#### RNA-seq Quantification Tools

**kallisto (Recommended):**
```bash
sudo apt-get install -y kallisto
kallisto version  # Should show: kallisto, version 0.48.0 or higher
```

**salmon (Alternative):**
```bash
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz
tar -xzf salmon-1.10.0_linux_x86_64.tar.gz
export PATH=$PATH:$(pwd)/salmon-latest_linux_x86_64/bin
salmon --version
```

#### Quality Control Tools
```bash
sudo apt-get install -y fastp seqkit

# Verify installations
fastp --version  # Should show 0.23.2 or higher
seqkit version   # Should show 2.3.0 or higher
```

### Environment Configuration

#### Required Environment Variables

```bash
# Add to ~/.bashrc or ~/.zshrc
export NCBI_EMAIL="your.email@example.com"  # Required for NCBI API

# Optional: Configure data directories
export METAINFORMANT_OUTPUT="/path/to/MetaInformAnt/output"
export METAINFORMANT_DATA="/path/to/MetaInformAnt/data"
```

## Running Workflows

### Single Species Workflow (Recommended)

**ENA-based workflow** (100% reliability):
```bash
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 12 \
  --threads 12
```

**Features:**
- Direct ENA downloads with 100% reliability
- Batched processing (12 samples at a time)
- Automatic resume with `wget --continue`
- Auto-cleanup (FASTQs deleted after quantification)

### Multi-Species Workflows

**Batch download** (configurable parallelism):
```bash
# Default: 3 species × 10 threads = 30 total downloads
python3 scripts/rna/batch_download_species.py

# Custom: 4 species × 12 threads = 48 total downloads
python3 scripts/rna/batch_download_species.py \
  --species-count 4 \
  --threads-per-species 12
```

**Multi-species SRA workflow** (legacy):
```bash
# Prerequisites: .venv must exist with amalgkit installed (see Installation section above)
# If not set up, run:
#   uv venv .venv  # or /tmp/metainformant_venv on ext6 filesystems
#   source .venv/bin/activate  # or /tmp/metainformant_venv/bin/activate
#   uv pip install -e .
#   uv pip install git+https://github.com/kfuku52/amalgkit

# Scripts auto-discover venv location (.venv or /tmp/metainformant_venv)
python3 scripts/rna/run_multi_species.py

# With configurable threads:
export AK_THREADS=12
python3 scripts/rna/run_multi_species.py
```

See [MULTI_SPECIES_QUICK_START.md](MULTI_SPECIES_QUICK_START.md) for detailed production workflows.

## Understanding the Output

### Directory Structure

```
output/amalgkit/
├── {species}/
│   ├── quant/              # Quantification results (~2 MB per sample)
│   │   ├── SRR1234567/
│   │   │   ├── abundance.tsv
│   │   │   └── run_info.json
│   │   └── ...
│   ├── fastq/              # Temporary FASTQs (auto-deleted)
│   ├── work/
│   │   ├── metadata/
│   │   └── index/          # Kallisto index
│   └── logs/               # Workflow logs
```

### Key Files

- **Per-sample quantification**: `output/amalgkit/{species}/quant/{SRR_ID}/abundance.tsv`
- **Metadata**: `output/amalgkit/{species}/work/metadata/metadata.tsv`
- **Kallisto index**: `output/amalgkit/{species}/work/index/transcriptome.idx`
- **Logs**: `output/amalgkit/{species}/logs/`

## Next Steps

1. **Learn about workflows**: See [WORKFLOW.md](WORKFLOW.md) for planning and execution
2. **Understand steps**: See [STEPS.md](STEPS.md) for individual workflow steps
3. **Configure workflows**: See [CONFIGURATION.md](CONFIGURATION.md) for configuration options
4. **Production workflows**: See [MULTI_SPECIES_QUICK_START.md](MULTI_SPECIES_QUICK_START.md)
5. **Discover species**: See [DISCOVERY/GUIDE.md](DISCOVERY/GUIDE.md) for species discovery

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
cd /path/to/MetaInformAnt
python3 scripts/rna/workflow_ena_integrated.py ...
```

### PEP 668 "externally-managed-environment" error
This is expected on modern Linux systems. Always use virtual environments:
```bash
source .venv/bin/activate  # Always activate before installing packages
```

### Downloads Failing
- Check network connectivity: `ping -c 4 8.8.8.8`
- Verify NCBI_EMAIL is set: `echo $NCBI_EMAIL`
- Try reducing threads: `--threads 8`
- Check disk space: `df -h /`

### Disk Space Issues
- Clean up partial downloads: `python3 scripts/rna/cleanup_partial_downloads.py --execute`
- Reduce batch size: `--batch-size 6`
- Check for large temp files: `du -sh output/amalgkit/*/fastq`

## Getting Help

- **Documentation**: See [README.md](README.md) for navigation
- **Workflow Issues**: See [WORKFLOW.md](WORKFLOW.md#common-issues-and-solutions)
- **Performance**: See [MULTI_SPECIES_QUICK_START.md](MULTI_SPECIES_QUICK_START.md#performance-characteristics)
- **Orchestration**: See [ORCHESTRATION/README.md](ORCHESTRATION/README.md)

## Related Documentation

- **[SETUP.md](SETUP.md)**: Detailed setup instructions (superseded by this guide)
- **[WORKFLOW.md](WORKFLOW.md)**: Workflow planning and execution
- **[STEPS.md](STEPS.md)**: Individual workflow steps
- **[CONFIGURATION.md](CONFIGURATION.md)**: Configuration management
- **[ORCHESTRATION/README.md](ORCHESTRATION/README.md)**: Orchestrator overview

## See Also

- **Source Code**: `src/metainformant/rna/` - Implementation details
- **Module Documentation**: `src/metainformant/rna/README.md` - API reference
- **Tests**: `tests/test_rna_*.py` - Comprehensive test coverage

