# Amalgkit RNA Analysis Documentation

This directory contains comprehensive documentation for METAINFORMANT's integration with the amalgkit RNA-seq analysis toolkit.

## Overview

The amalgkit integration provides a complete transcriptomic analysis pipeline from raw SRA data to publication-ready results, with intelligent error recovery and progress monitoring. Key features:

- **Direct ENA Downloads**: Production workflow bypasses SRA Toolkit with 100% reliability
- **Robust Retry Logic**: wget-based downloads with automatic resume capability
- **Auto-activation**: Scripts automatically activate virtual environments
- **Immediate Processing**: Per-sample processing (download → immediately quantify → immediately delete FASTQs)
- **Thread Allocation**: Total threads distributed across all species (default: 24 total, not per species)
- **Zero-Read Detection**: Enhanced validation detects when getfastq produces 0 reads but reports success, preventing false positives

## Documentation Files

### Core Amalgkit Integration
- **`amalgkit.md`**: Complete transcriptomic analysis pipeline documentation (includes advanced usage)
- **`README.md`**: This file - overview and quick start
- **`FUNCTIONS.md`**: Quick function lookup table

### Workflow Guides
- **Getting Started**: See [GETTING_STARTED.md](../GETTING_STARTED.md) for complete setup and workflow guide
- **`R_INSTALLATION.md`**: R installation and setup guide
- **`r_packages.md`**: R package setup and troubleshooting

### Genome Setup
- **`genome_preparation.md`**: Technical documentation for genome download and kallisto index building
- **`genome_setup_guide.md`**: Complete step-by-step guide for setting up genomes and indexes
- **`commands.md`**: Command reference for all genome setup scripts

**Genome Setup Scripts** (located in `scripts/rna/`):
- **`verify_genomes_and_indexes.py`**: Check status of genome downloads and kallisto indexes
- **`download_missing_genomes.py`**: Download missing genome packages from NCBI
- **`prepare_transcriptomes.py`**: Extract and prepare RNA FASTA files from genomes
- **`build_kallisto_indexes.py`**: Build kallisto indexes from transcriptome FASTA files
- **`orchestrate_genome_setup.py`**: Master orchestrator for complete genome setup pipeline
- **`run_genome_setup.sh`**: Shell script to run all steps sequentially

### Step Documentation
- **`steps/`**: Complete documentation for all 11 amalgkit steps
  - metadata, integrate, config, select, getfastq, quant, merge, cstmm, curate, csca, sanity

### Validation
- **`testing_coverage.md`**: Testing coverage, validation, and production results

## Related Source Code

- See `src/metainformant/rna/` for implementation details
- See `src/metainformant/rna/amalgkit.py` for CLI wrapper implementation
- See `src/metainformant/rna/workflow.py` for workflow orchestration
- See `tests/test_rna_*.py` for comprehensive test coverage

## Quick Links

- **[API Reference](../API.md)** - Complete function documentation
- **[Function Index](FUNCTIONS.md)** - Quick function lookup
- **[Step Documentation](steps/README.md)** - All 11 step guides
- **[Pipeline Overview](amalgkit.md)** - Complete pipeline documentation
- **[Main Index](../README.md)** - RNA domain master index

## See Also

### Documentation
- **[API Reference](../API.md)** - Complete function documentation
- **[Function Index](FUNCTIONS.md)** - Quick function lookup
- **[Step Documentation](steps/README.md)** - All 11 step guides
- **[Pipeline Overview](amalgkit.md)** - Complete pipeline documentation

### Related Topics
- **[Workflow Guide](../workflow.md)** - Workflow planning and execution
- **[Configuration Guide](../CONFIGURATION.md)** - Configuration management
- **[Orchestration Guide](../ORCHESTRATION.md)** - Orchestrator overview
- **[Getting Started](../GETTING_STARTED.md)** - Setup and installation

### Genome Setup
- **[genome_setup_guide.md](genome_setup_guide.md)** - User guide (step-by-step)
- **[genome_preparation.md](genome_preparation.md)** - Technical API reference
- **[commands.md](commands.md)** - Command reference

## Usage Examples

The amalgkit integration supports complete RNA-seq workflows:

```python
from metainformant.rna import amalgkit, workflow

# Check amalgkit availability
ok, help_text = amalgkit.check_cli_available()

# Execute complete workflow (12 threads default)
cfg = workflow.AmalgkitWorkflowConfig(
    work_dir="output/amalgkit/work",
    species_list=["Homo sapiens"],
    threads=12
)
results = workflow.execute_workflow(cfg)
```

**Command-line usage** (recommended for end-to-end workflows):
```bash
# Full end-to-end workflow (all steps)
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml

# Specific steps only
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --steps getfastq quant merge

# Check status
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --status

# Print planned steps + exact amalgkit commands (does not execute)
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml --plan
```

The `run_workflow.py` script provides complete end-to-end execution via `execute_workflow()`, including automatic genome setup, per-sample processing, and all 11 amalgkit steps.

See [scripts/rna/README.md](../../scripts/rna/README.md) for complete usage documentation.

## Integration

Amalgkit integration connects with:
- **ENA database** for direct FASTQ retrieval (production)
- **NCBI SRA database** for metadata and legacy downloads
- **Kallisto** for pseudoalignment quantification
- **R environment** for statistical analysis
- **Visualization tools** for publication figures

## Testing

Tests cover:
- CLI availability checks
- Step runner execution (real `amalgkit` when available)
- Error handling paths and manifest/report generation

External-tool tests are skipped when `amalgkit` is not on `PATH`.
To attempt an install during a test run, set:
`METAINFORMANT_AK_AUTO_INSTALL=1`.

**Disk space management:**
- Direct ENA downloads with automatic retry and resume
- Immediate per-sample processing: only one sample's FASTQs exist at a time
- Automatic FASTQ deletion immediately after quantification (maximum disk efficiency)
- Final results: ~40-55 GB for 20 species (4,548 samples total)
- No /tmp partition limitations

**Performance:**
- Total threads distributed across all species (default: 24 total threads)
- Thread allocation: evenly distributed with minimum 1 thread per species
- Dynamic redistribution as species complete
- Immediate processing: download → quant → delete per sample (not batched)
- Direct ENA downloads: 100% reliability vs 0% SRA Toolkit
- Multi-species coordination and cross-species analysis
- Automatic virtual environment activation

## Contributing

When enhancing amalgkit integration:
1. Update workflow and pipeline documentation
2. Add comprehensive integration tests
3. Ensure compatibility with amalgkit ecosystem
4. Update success metrics and benchmarks

## Example Dataset

Complete working example with *Pogonomyrmex barbatus*:
- **Location**: `output/amalgkit/pbarbatus/`
- **Samples**: 83 brain RNA-seq runs
- **Status**: Complete workflow with all visualizations
- **Results**: Quantification files, expression matrices, and QC reports

**Note**: Any markdown files in `output/` are program-generated execution outputs, not documentation. See [EXAMPLES.md](../EXAMPLES.md) for documentation of this analysis.

## Quick Start for New Species

1. **Setup genomes and indexes**: See [genome_setup_guide.md](genome_setup_guide.md) for complete genome setup
   - Verify genomes: `python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_<species>.yaml --verify-only`
   - Full genome setup: `python3 scripts/rna/setup_genome.py --config config/amalgkit/amalgkit_<species>.yaml`
   - Note: Genome setup happens automatically when running workflows if genome config exists
2. **Review workflow guide**: [GETTING_STARTED.md](../GETTING_STARTED.md)
3. **Setup R environment**: `R_INSTALLATION.md` and `r_packages.md`
4. **Run end-to-end workflow**: `python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_<species>.yaml`
5. **Validate outputs**: Use `--status` flag or `scripts/rna/amalgkit/verify_workflow.sh <species>`

**All setup uses `uv` for package management** (required):
- Create venv: `uv venv`
- Install packages: `uv pip install -e . --python .venv/bin/python3`
- Install amalgkit: `uv pip install git+https://github.com/kfuku52/amalgkit --python .venv/bin/python3`

**No manual venv activation needed** - scripts automatically discover and activate virtual environments (`.venv` or `/tmp/metainformant_venv`). Setup uses `uv` for reliable package management, with automatic fallback to `/tmp/metainformant_venv` on ext6 filesystems that don't support symlinks.

**Disk space requirements:**
- Minimal peak usage: only one sample's FASTQs exist at a time
- Immediate deletion after quantification ensures maximum disk efficiency
- Direct ENA downloads with automatic resume
- Final results: ~2-3 GB per species (expression matrices + QC)

## Related Documentation

- See `docs/rna/README.md` for RNA domain overview
- See `docs/rna/workflow.md` for workflow orchestration details
- See `docs/rna/amalgkit/steps/README.md` for individual step documentation
- See `output/amalgkit/pbarbatus/` for complete working example
- See `scripts/rna/README.md` for production workflow scripts

This documentation provides complete coverage of METAINFORMANT's amalgkit integration capabilities for any species.
