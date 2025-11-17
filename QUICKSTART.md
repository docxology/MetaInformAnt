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

### Option 2: Manual Installation with UV

```bash
# Clone repository
git clone https://github.com/q/MetaInformAnt.git
cd MetaInformAnt

# Install uv if not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create virtual environment
uv venv

# Activate virtual environment
source .venv/bin/activate

# Install package
uv pip install -e .
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

### Protein Analysis

```python
from metainformant.protein import sequences, alignment

# Read protein sequences
proteins = sequences.read_fasta("data/proteins.fasta")

# Pairwise alignment
align_result = alignment.global_align(proteins["seq1"], proteins["seq2"])
print(f"Alignment score: {align_result.score}")
```

### Quality Control

```python
from metainformant.quality import fastq

# Assess FASTQ quality
qc_report = fastq.assess_quality("data/reads.fastq")
print(f"Mean quality: {qc_report['mean_quality']}")
print(f"Total reads: {qc_report['total_reads']}")
```

### Visualization

```python
from metainformant.visualization import lineplot
import matplotlib.pyplot as plt

# Create a simple line plot
data = [1, 2, 3, 4, 5]
ax = lineplot(None, data)
ax.set_ylabel("Values")
ax.set_title("Example Plot")
plt.savefig("output/example_plot.png", dpi=300)
```

### RNA-seq Workflow

```bash
# Check if amalgkit is available
python -c "from metainformant.rna import check_cli_available; print(check_cli_available())"

# Run end-to-end workflow for a single species (recommended)
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Check workflow status
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status

# Or use CLI
uv run metainformant rna run --work-dir output/rna --threads 8 --species Apis_mellifera
```

### CLI Workflows

All modules are accessible via the unified CLI:

```bash
# DNA analysis
uv run metainformant dna fetch --assembly GCF_000001405.40
uv run metainformant dna align --input data/sequences.fasta --output output/dna/alignment
uv run metainformant dna variants --input data/variants.vcf --format vcf --output output/dna/variants

# RNA analysis (see RNA-seq Workflow section above for more details)
uv run metainformant rna run --work-dir output/rna --threads 8 --species Apis_mellifera
uv run metainformant rna run-config --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Protein analysis
uv run metainformant protein taxon-ids --file data/taxon_ids.txt
uv run metainformant protein comp --fasta data/proteins.fasta
uv run metainformant protein rmsd-ca --pdb-a data/structure1.pdb --pdb-b data/structure2.pdb

# Epigenome analysis
uv run metainformant epigenome run --methylation data/methylation.tsv --output output/epigenome

# Ontology analysis
uv run metainformant ontology run --go data/go.obo --output output/ontology

# Phenotype analysis
uv run metainformant phenotype run --input data/traits.csv --output output/phenotype --analyze-statistics

# Ecology analysis
uv run metainformant ecology run --input data/species.csv --output output/ecology --diversity

# Mathematical biology
uv run metainformant math popgen --input data/sequences.fasta --output output/math/popgen
uv run metainformant math coalescent --n-samples 10 --output output/math/coalescent

# GWAS analysis
uv run metainformant gwas run --config config/gwas/gwas_template.yaml

# Information theory
uv run metainformant information entropy --input data/seqs.fasta --output output/information
uv run metainformant information mutual-information --x data/x.txt --y data/y.txt --output output/information

# Life events analysis
uv run metainformant life-events embed --input data/events.json --output output/life_events/embeddings
uv run metainformant life-events predict --events data/events.json --model output/life_events/model.json --output output/life_events/predictions

# Visualization
uv run metainformant visualization run --input data/matrix.csv --plot-type heatmap --output output/visualization

# Simulation
uv run metainformant simulation run --model sequences --output output/simulation

# Network analysis
uv run metainformant networks run --input data/interactions.tsv --output output/networks --analyze-metrics

# Multi-omics integration
uv run metainformant multiomics run --genomics data/genomics.tsv --transcriptomics data/rna.tsv --output output/multiomics --joint-pca

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
python scripts/core/run_demo.py

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
- **[Core Utilities](docs/core/README.md)** - Shared infrastructure and utilities
- **[DNA Analysis](docs/dna/index.md)** - Sequence analysis workflows
- **[RNA-seq](docs/rna/index.md)** - Transcriptomics pipelines
- **[Protein Analysis](docs/protein/index.md)** - Protein sequences and structures
- **[Epigenome](docs/epigenome/index.md)** - Epigenetic modification analysis
- **[Ontology](docs/ontology/index.md)** - Functional annotation and ontologies
- **[Phenotype](docs/phenotype/index.md)** - Phenotypic trait analysis
- **[Ecology](docs/ecology/index.md)** - Ecological metadata and community analysis
- **[Mathematical Biology](docs/math/index.md)** - Mathematical and theoretical biology
- **[GWAS](docs/gwas/index.md)** - Association studies
- **[Information Theory](docs/information/index.md)** - Information-theoretic analysis
- **[Life Events](docs/life_events/index.md)** - Life course event analysis
- **[Visualization](docs/visualization/index.md)** - Plotting and visualization
- **[Simulation](docs/simulation/index.md)** - Synthetic data generation
- **[Single-Cell](docs/singlecell/index.md)** - scRNA-seq analysis
- **[Quality Control](docs/quality/index.md)** - Data quality assessment
- **[Network Analysis](docs/networks/index.md)** - Biological network analysis
- **[Machine Learning](docs/ml/index.md)** - ML methods
- **[Multi-Omics](docs/multiomics/index.md)** - Multi-omic data integration

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
# Reinstall all dependencies with uv
uv pip install -e . --python .venv/bin/python3
```

### Import Errors
```bash
# Ensure package is installed in development mode with uv
uv pip install -e . --python .venv/bin/python3
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
