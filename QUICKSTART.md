# MetaInformAnt Quick Start

Get started with MetaInformAnt in minutes.

## Prerequisites

- Python 3.11 or higher
- Git

## Installation

### Option 1: Quick Setup with UV (Recommended)

```bash
# Clone repository
git clone https://github.com/q/MetaInformAnt.git
cd MetaInformAnt

# Run automated setup script
bash scripts/package/setup_uv.sh

# Activate virtual environment
source .venv/bin/activate
```

### Option 2: Manual Installation

```bash
# Clone repository
git clone https://github.com/q/MetaInformAnt.git
cd MetaInformAnt

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install package
pip install -e .
```

## Verify Installation

```bash
# Check Python version
python --version

# Run tests
pytest tests/ -v

# Check package installation
python -c "import metainformant; print(metainformant.__version__)"
```

## Basic Usage Examples

### DNA Analysis

```python
from metainformant.dna import sequences, composition

# Read FASTA file
seqs = sequences.read_fasta("data/sequences.fasta")

# Calculate GC content
for name, seq in seqs.items():
    gc = composition.gc_content(seq)
    print(f"{name}: GC content = {gc:.2%}")
```

### RNA-seq Workflow

```bash
# Check if amalgkit is available
python -c "from metainformant.rna import check_cli_available; print(check_cli_available())"

# Run multi-species RNA-seq workflow
python scripts/rna/run_multi_species.py

# Or use CLI
uv run metainformant rna run --work-dir output/rna --threads 8 --species Apis_mellifera
```

### CLI Workflows

All modules are accessible via the unified CLI:

```bash
# DNA analysis
uv run metainformant dna fetch --assembly GCF_000001405.40

# Network analysis
uv run metainformant networks run --input data/interactions.tsv --output output/networks

# Multi-omics integration
uv run metainformant multiomics run --genomics data/genomics.tsv --transcriptomics data/rna.tsv --output output/multiomics

# Single-cell analysis
uv run metainformant singlecell run --input data/counts.h5ad --output output/singlecell --qc --normalize

# Quality control
uv run metainformant quality run --fastq data/reads.fq --output output/quality --analyze-fastq

# Machine learning
uv run metainformant ml run --features data/features.csv --labels data/labels.csv --output output/ml --classify

# See all commands
uv run metainformant --help
```

See [`docs/cli.md`](docs/cli.md) for complete CLI reference.

### Complete Demonstration

```bash
# Run comprehensive workflow demo
python scripts/run_complete_demo.py

# View results in output/demo/ directory
ls output/demo/
```

## Directory Structure

MetaInformAnt follows a clean directory policy:

- **`config/`** - Configuration files (YAML/TOML)
- **`data/`** - Input datasets and databases (read-only)
- **`output/`** - All analysis outputs (safe to delete/regenerate)
- **`src/metainformant/`** - Main package source code
- **`tests/`** - Test suite
- **`scripts/`** - Workflow orchestration scripts
- **`docs/`** - Comprehensive documentation

## Optional External Tools

Some workflows require external tools:

### RNA Analysis
- **amalgkit**: RNA-seq workflow orchestration
  ```bash
  # Install via UV (automatic with setup_uv.sh --with-amalgkit)
  uv pip install amalgkit
  ```

### GWAS Analysis
- **SRA Toolkit**: For downloading sequencing data
  ```bash
  sudo apt-get install sra-toolkit
  ```
- **BWA, samtools, bcftools**: For alignment and variant calling
  ```bash
  sudo apt-get install bwa samtools bcftools
  ```

### Sequence Alignment
- **MUSCLE or ClustalO**: For multiple sequence alignment
  ```bash
  sudo apt-get install muscle clustalo
  ```

## Environment Variables

Set optional environment variables for enhanced functionality:

```bash
# NCBI E-utilities (for data download)
export NCBI_EMAIL="your.email@example.com"

# Add to ~/.bashrc for persistence
echo 'export NCBI_EMAIL="your.email@example.com"' >> ~/.bashrc
```

## Next Steps

### Documentation
- **[Documentation Guide](docs/DOCUMENTATION_GUIDE.md)** - Complete navigation guide
- **[Architecture Overview](docs/architecture.md)** - System design
- **[Testing Guide](docs/testing.md)** - Running tests

### Module-Specific Guides
- **[DNA Analysis](docs/dna/index.md)** - Sequence analysis workflows
- **[RNA-seq](docs/rna/index.md)** - Transcriptomics pipelines
- **[GWAS](docs/gwas/index.md)** - Association studies
- **[Single-Cell](docs/singlecell/index.md)** - scRNA-seq analysis
- **[Machine Learning](docs/ml/index.md)** - ML methods

### Workflow Scripts
- **[Scripts README](scripts/README.md)** - All available scripts
- **[RNA Workflows](scripts/rna/README.md)** - RNA-seq pipelines
- **[GWAS Workflows](scripts/gwas/)** - Genome-wide studies

## Common Commands

```bash
# Run all tests
bash scripts/package/run_tests.sh

# Run fast tests only
bash scripts/package/run_tests.sh --fast

# Check code quality
bash scripts/package/uv_quality.sh

# Update documentation
bash scripts/package/uv_docs.sh
```

## Troubleshooting

### Virtual Environment Not Activated
```bash
source .venv/bin/activate
```

### Missing Dependencies
```bash
# Reinstall all dependencies
pip install -e .
```

### Import Errors
```bash
# Ensure package is installed in development mode
pip install -e .
```

### Permission Errors
```bash
# Make scripts executable
chmod +x scripts/**/*.sh
```

## Getting Help

- **Documentation**: See [docs/](docs/) directory for comprehensive guides
- **Examples**: Check [scripts/](scripts/) for working examples
- **Issues**: Report issues at https://github.com/q/MetaInformAnt/issues

---

**Ready to start analyzing? Pick a module and dive into the documentation!**
