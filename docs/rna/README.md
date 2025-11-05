# RNA Domain Documentation

Comprehensive documentation for METAINFORMANT's RNA analysis capabilities, including transcriptomic workflows, species discovery, and multi-species orchestration.

## Current Status (November 2025)

🎯 **Ant Species Discovery Complete**: 55 species discovered, 20 validated genomes  
🔄 **Batch 1 Running**: 10 species (3,820 samples) processing  
⏳ **Batch 2 Queued**: 10 species (728 samples) ready to launch  
📊 **Total Pipeline**: 4,548 samples across 20 ant species  

See [Discovery Results](DISCOVERY/README.md) for complete details.

---

## Getting Started

New to RNA-seq analysis with METAINFORMANT? Start here:

- **[Getting Started Guide](GETTING_STARTED.md)** ⭐ - Installation, setup, and first workflow
- **[Multi-Species Quick Start](MULTI_SPECIES_QUICK_START.md)** ⭐⭐ - Production workflow guide
- **[Run All Species](RUN_ALL_SPECIES.md)** ⭐⭐⭐ - Complete guide for all 24 species with configurable threads
- **[Workflow Orchestration](ORCHESTRATION/README.md)** - Choose the right orchestrator for your needs

### Quick Commands

```bash
# All species, end-to-end workflow (configurable threads)
# Prerequisites: venv must exist with amalgkit installed (see RUN_ALL_SPECIES.md)
# See docs/rna/RUN_ALL_SPECIES.md for complete setup instructions
export AK_THREADS=12  # Optional: set threads globally
python3 scripts/rna/run_all_species_parallel.py  # Recommended: parallel execution
# Or: python3 scripts/rna/run_multi_species.py  # Sequential execution

# Single species workflow (ENA-based, recommended)
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 12 --threads 12

# Multi-species batch download (getfastq + quant only)
python3 scripts/rna/batch_download_species.py \
  --species-count 3 --threads-per-species 10

# Discover ant species with RNA-seq data
python3 scripts/rna/discover_ant_rnaseq_by_genus.py
```

**Quick Start for All Species**: See [RUN_ALL_SPECIES.md](RUN_ALL_SPECIES.md) for complete guide.

---

## Core Concepts

### Workflow Execution
- **[WORKFLOW.md](WORKFLOW.md)** - Planning and executing RNA-seq workflows
- **[STEPS.md](STEPS.md)** - Individual workflow steps (metadata, getfastq, quant, etc.)
- **[CONFIGURATION.md](CONFIGURATION.md)** - Configuration management and species profiles

### Amalgkit Integration
- **[amalgkit/README.md](amalgkit/README.md)** - Overview of amalgkit wrapper
- **[amalgkit/steps/](amalgkit/steps/)** - Detailed step documentation (11 steps)
- **[amalgkit/testing_coverage.md](amalgkit/testing_coverage.md)** - Test coverage and validation

### Module Reference
- **[index.md](index.md)** - RNA domain overview and module index
- **Source Code**: `src/metainformant/rna/` - Implementation details
- **Tests**: `tests/test_rna_*.py` - Comprehensive test coverage

---

## Advanced Topics

### Species Discovery
- **[DISCOVERY/GUIDE.md](DISCOVERY/GUIDE.md)** - Comprehensive discovery system documentation
- **[DISCOVERY/QUICK_REF.md](DISCOVERY/QUICK_REF.md)** - Quick reference for discovery
- **[DISCOVERY/README.md](DISCOVERY/README.md)** - Discovery results and current status

### Orchestration
- **[ORCHESTRATION/README.md](ORCHESTRATION/README.md)** - Overview of all orchestrators
- **[ORCHESTRATION/ENA_WORKFLOW.md](ORCHESTRATION/ENA_WORKFLOW.md)** - ENA-based workflow documentation
- **[ORCHESTRATION/MULTI_SPECIES.md](ORCHESTRATION/MULTI_SPECIES.md)** - Multi-species SRA workflow
- **[ORCHESTRATION/BATCH_DOWNLOAD.md](ORCHESTRATION/BATCH_DOWNLOAD.md)** - Configurable batch download

### Examples
- **[examples/](examples/)** - Real-world analysis examples and case studies

---

## Reference

### API Documentation
- **Source Module**: `src/metainformant/rna/README.md` - Module-specific API reference
- **Workflow Functions**: See [WORKFLOW.md](WORKFLOW.md) for `plan_workflow`, `execute_workflow`
- **Step Runners**: See [STEPS.md](STEPS.md) for individual step functions

### Troubleshooting
- **Common Issues**: See [WORKFLOW.md](WORKFLOW.md#common-issues-and-solutions)
- **Setup Problems**: See [GETTING_STARTED.md](GETTING_STARTED.md#troubleshooting)
- **Performance**: See [MULTI_SPECIES_QUICK_START.md](MULTI_SPECIES_QUICK_START.md#performance-characteristics)

### Testing
- **Test Coverage**: [amalgkit/testing_coverage.md](amalgkit/testing_coverage.md)
- **Running Tests**: `pytest tests/test_rna_*.py -v`
- **Integration Tests**: See test files for comprehensive coverage
- **End-to-End Validation**: [END_TO_END_VALIDATION.md](END_TO_END_VALIDATION.md) - Validate all 24 species ready for workflows

---

## Overview

The RNA domain provides tools for transcriptomic analysis, workflow orchestration, and integration with external bioinformatics tools like amalgkit. Key features include:

- **Species Discovery**: Automated identification of species with RNA-seq data
- **Genome Integration**: Validated genome assemblies from NCBI
- **Direct ENA Downloads**: Robust FASTQ retrieval using ENA API, bypassing problematic SRA Toolkit
- **Retry Logic**: wget-based downloads with --continue for automatic resumption
- **Batched Processing**: Disk-friendly processing of large RNA-seq cohorts (download→quantify→delete)
- **Auto-activation**: Scripts automatically detect and activate virtual environments
- **Multi-species Support**: Coordinated analysis across multiple species
- **Production-Ready**: Currently processing 3,820 samples in Batch 1 (10 ant species)

## Integration

RNA analysis integrates with:
- **DNA sequences** for genomic context
- **Protein annotations** for functional analysis
- **Statistical methods** for differential expression
- **Visualization tools** for expression plots

See [WORKFLOW.md](WORKFLOW.md) for integration examples and patterns.

## Navigation Map

### By Use Case
- **New User**: [GETTING_STARTED.md](GETTING_STARTED.md) → [WORKFLOW.md](WORKFLOW.md) → [STEPS.md](STEPS.md)
- **Production Workflows**: [MULTI_SPECIES_QUICK_START.md](MULTI_SPECIES_QUICK_START.md) → [ORCHESTRATION/README.md](ORCHESTRATION/README.md)
- **Species Discovery**: [DISCOVERY/GUIDE.md](DISCOVERY/GUIDE.md) → [DISCOVERY/QUICK_REF.md](DISCOVERY/QUICK_REF.md)
- **Configuration**: [CONFIGURATION.md](CONFIGURATION.md) → [ORCHESTRATION/README.md](ORCHESTRATION/README.md)
- **API Reference**: [index.md](index.md) → [WORKFLOW.md](WORKFLOW.md) → [STEPS.md](STEPS.md)

### By Component
- **Core Module**: [index.md](index.md), [WORKFLOW.md](WORKFLOW.md), [STEPS.md](STEPS.md), [CONFIGURATION.md](CONFIGURATION.md)
- **Amalgkit**: [amalgkit/README.md](amalgkit/README.md), [amalgkit/steps/](amalgkit/steps/)
- **Orchestration**: [ORCHESTRATION/README.md](ORCHESTRATION/README.md)
- **Discovery**: [DISCOVERY/](DISCOVERY/)
- **Examples**: [examples/](examples/)

## Current Project Timeline (November 2025)

### Phase 1: Discovery ✅ COMPLETE
- **Duration**: ~80 minutes
- **Outcome**: 55 species discovered, 20 genomes validated, 20 configs generated
- **Date**: November 3, 2025

### Phase 2: Batch 1 Execution 🔄 IN PROGRESS
- **Species**: 10 (top by sample count)
- **Samples**: 3,820 
- **Started**: November 3, 2025 16:18 PST
- **Current Stage**: Downloading FASTQs
- **ETA**: November 4-5, 2025 (24-48 hours)

### Phase 3: Batch 2 Execution ⏳ QUEUED
- **Species**: 10 (remaining validated)
- **Samples**: 728
- **Launch**: After Batch 1 completes
- **ETA**: 12-24 hours after launch

### Phase 4: Analysis & Validation ⏳ FUTURE
- **Tasks**: Expression matrix review, QC analysis, comparative genomics
- **ETA**: After both batches complete

**Total Project Duration**: ~36-72 hours from discovery to complete analysis

## See Also

- **Source Code**: `src/metainformant/rna/` - Implementation details
- **Module Documentation**: `src/metainformant/rna/README.md` - API reference
- **Tests**: `tests/test_rna_*.py` - Comprehensive test coverage
- **Production Validation**: [amalgkit/testing_coverage.md](amalgkit/testing_coverage.md)

## Contributing

When adding new RNA analysis functionality:
1. Update workflow and step documentation
2. Add comprehensive integration tests
3. Update configuration templates
4. Ensure compatibility with amalgkit ecosystem
5. Update this README if adding new documentation sections

---

**Last Updated**: November 2025  
**Documentation Status**: Complete and organized  
**Production Status**: Validated with 4,548 samples across 20 ant species
