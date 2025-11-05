# Amalgkit RNA Analysis Documentation

This directory contains comprehensive documentation for METAINFORMANT's integration with the amalgkit RNA-seq analysis toolkit.

## Overview

The amalgkit integration provides a complete transcriptomic analysis pipeline from raw SRA data to publication-ready results, with intelligent error recovery and progress monitoring. Key features:

- **Direct ENA Downloads**: Production workflow bypasses SRA Toolkit with 100% reliability
- **Robust Retry Logic**: wget-based downloads with automatic resume capability
- **Auto-activation**: Scripts automatically activate virtual environments
- **Batched Processing**: 12-sample batches for disk-friendly workflows
- **Multi-threading**: 12 parallel threads by default

## Documentation Files

### Core Amalgkit Integration
- **`amalgkit.md`**: Complete transcriptomic analysis pipeline documentation
- **`comprehensive_guide.md`**: Comprehensive amalgkit workflow guide
- **`README.md`**: This file - overview and quick start

### Workflow Guides
- **`quick_start.md`**: Quick start guide for sanity and curate steps (all species)
- **`R_INSTALLATION.md`**: R installation and setup guide
- **`r_packages.md`**: R package setup and troubleshooting

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

## See Also

- **[../WORKFLOW.md](../WORKFLOW.md)**: Workflow planning and execution
- **[../STEPS.md](../STEPS.md)**: Individual step documentation
- **[../CONFIGURATION.md](../CONFIGURATION.md)**: Configuration management
- **[../ORCHESTRATION/README.md](../ORCHESTRATION/README.md)**: Orchestrator overview
- **[../GETTING_STARTED.md](../GETTING_STARTED.md)**: Setup and installation

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
- Batched processing: ~50-100 GB peak disk usage (temporary, auto-cleaned)
- Automatic FASTQ cleanup after quantification (saves ~2-10 GB per sample)
- Final results: ~40-55 GB for 20 species (4,548 samples total)
- No /tmp partition limitations

**Performance:**
- 12 parallel threads for downloads and quantification
- Direct ENA downloads: 100% reliability vs 0% SRA Toolkit
- Multi-species coordination and cross-species analysis
- Automatic virtual environment activation
- Currently processing 3,820 samples in Batch 1 (10 ant species, November 2025)
- Batch 2 queued: 728 samples (10 additional species)

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

1. **Review workflow guide**: `quick_start.md`
2. **Setup R environment**: `R_INSTALLATION.md` and `r_packages.md`
3. **Run ENA workflow**: See `scripts/rna/workflow_ena_integrated.py`
4. **For batch processing**: Use `scripts/rna/run_all_species_parallel.py` for parallel execution of all species
5. **Validate outputs**: Run sanity and curate steps

**No manual venv activation needed** - scripts automatically discover and activate virtual environments (`.venv` or `/tmp/metainformant_venv`). Setup uses `uv` for reliable package management, with automatic fallback to `/tmp/metainformant_venv` on ext6 filesystems that don't support symlinks.

**Disk space requirements:**
- ~50-100 GB peak usage (temporary, during download/quant)
- Direct ENA downloads with automatic resume
- Automatic cleanup after quantification (FASTQs deleted immediately)
- Final results: ~2-3 GB per species (expression matrices + QC)

## Related Documentation

- See `docs/rna/README.md` for RNA domain overview
- See `docs/rna/workflow.md` for workflow orchestration details
- See `docs/rna/steps.md` for individual step documentation
- See `output/amalgkit/pbarbatus/` for complete working example
- See `scripts/rna/README.md` for production workflow scripts

This documentation provides complete coverage of METAINFORMANT's amalgkit integration capabilities for any species.
