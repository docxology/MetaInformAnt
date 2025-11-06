# Amalgkit RNA Analysis Documentation

This directory contains comprehensive documentation for METAINFORMANT's integration with the amalgkit RNA-seq analysis toolkit.

## Overview

The amalgkit integration provides a complete transcriptomic analysis pipeline from raw SRA data to publication-ready results, with intelligent error recovery and progress monitoring. Key features:

- **Direct ENA Downloads**: Production workflow bypasses SRA Toolkit with 100% reliability
- **Robust Retry Logic**: wget-based downloads with automatic resume capability
- **Auto-activation**: Scripts automatically activate virtual environments
- **Immediate Processing**: Per-sample processing (download → immediately quantify → immediately delete FASTQs)
- **Thread Allocation**: Total threads distributed across all species (default: 24 total, not per species)

## Documentation Files

### Core Amalgkit Integration
- **`amalgkit.md`**: Complete transcriptomic analysis pipeline documentation (includes advanced usage)
- **`README.md`**: This file - overview and quick start
- **`FUNCTIONS.md`**: Quick function lookup table

### Workflow Guides
- **`quick_start.md`**: Quick start guide for sanity and curate steps (all species)
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
- **[Pipeline Overview](amalgkit.md)** - Complete pipeline documentation (includes advanced usage)

### Related Topics
- **[Workflow Guide](../workflow.md)** - Workflow planning and execution
- **[Configuration Guide](../CONFIGURATION.md)** - Configuration management
- **[Orchestration Guide](../ORCHESTRATION/README.md)** - Orchestrator overview
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

**Command-line usage** (production ENA workflow):
```bash
# Production workflow with ENA direct downloads
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 12 \
  --threads 12
```

**Legacy SRA-based workflow**:
```bash
# Prerequisites: .venv must exist with amalgkit installed
# If not set up, run:
#   uv venv .venv  # or /tmp/metainformant_venv on ext6 filesystems
#   source .venv/bin/activate  # or /tmp/metainformant_venv/bin/activate
#   uv pip install -e .
#   uv pip install git+https://github.com/kfuku52/amalgkit

# Alternative workflow using SRA Toolkit (scripts auto-discover venv location)
python3 scripts/rna/run_multi_species.py

# With configurable threads:
export AK_THREADS=12
python3 scripts/rna/run_multi_species.py
```

## Integration

Amalgkit integration connects with:
- **ENA database** for direct FASTQ retrieval (production)
- **NCBI SRA database** for metadata and legacy downloads
- **Kallisto** for pseudoalignment quantification
- **R environment** for statistical analysis
- **Visualization tools** for publication figures

## Testing

Comprehensive tests ensure workflow reliability:
- CLI tool availability validation
- Workflow execution testing
- Error handling and recovery verification
- Integration testing with real data

**Disk Space Management:**
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
- **Documentation**: Species-specific quick reference and analysis report

See `output/amalgkit/pbarbatus/QUICK_REFERENCE.md` for immediate usage.

## Quick Start for New Species

1. **Setup genomes and indexes**: See [genome_setup_guide.md](genome_setup_guide.md) for complete genome setup
   - Verify genomes: `python3 scripts/rna/verify_genomes_and_indexes.py`
   - Download missing genomes: `python3 scripts/rna/download_missing_genomes.py`
   - Prepare transcriptomes: `python3 scripts/rna/prepare_transcriptomes.py`
   - Build kallisto indexes: `python3 scripts/rna/build_kallisto_indexes.py`
2. **Review workflow guide**: `quick_start.md`
3. **Setup R environment**: `R_INSTALLATION.md` and `r_packages.md`
4. **Run ENA workflow**: See `scripts/rna/workflow_ena_integrated.py`
5. **For batch processing**: Use `scripts/rna/run_all_species_parallel.py` for parallel execution of all species
6. **Validate outputs**: Run sanity and curate steps

**No manual venv activation needed** - scripts automatically discover and activate virtual environments (`.venv` or `/tmp/metainformant_venv`). Setup uses `uv` for reliable package management, with automatic fallback to `/tmp/metainformant_venv` on ext6 filesystems that don't support symlinks.

**Disk space requirements:**
- Minimal peak usage: only one sample's FASTQs exist at a time
- Immediate deletion after quantification ensures maximum disk efficiency
- Direct ENA downloads with automatic resume
- Final results: ~2-3 GB per species (expression matrices + QC)

## Related Documentation

- See `docs/rna/README.md` for RNA domain overview
- See `docs/rna/workflow.md` for workflow orchestration details
- See `docs/rna/steps.md` for individual step documentation
- See `output/amalgkit/pbarbatus/` for complete working example
- See `scripts/rna/README.md` for production workflow scripts

This documentation provides complete coverage of METAINFORMANT's amalgkit integration capabilities for any species.
