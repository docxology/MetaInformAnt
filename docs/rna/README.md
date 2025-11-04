# RNA Domain Documentation

This directory contains comprehensive documentation for METAINFORMANT's RNA analysis capabilities.

## Current Status (November 2025)

🎯 **Ant Species Discovery Complete**: 55 species discovered, 20 validated genomes  
🔄 **Batch 1 Running**: 10 species (3,820 samples) processing  
⏳ **Batch 2 Queued**: 10 species (728 samples) ready to launch  
📊 **Total Pipeline**: 4,548 samples across 20 ant species  

See `discovery/` for complete details.

## Quick Start

**→ [Multi-Species Quick Start Guide](MULTI_SPECIES_QUICK_START.md)** ⭐

Complete guide for starting and monitoring multi-species RNA-seq workflows:
- Step-by-step workflow launching
- Real-time progress monitoring
- Troubleshooting and optimization
- Expected timelines and resource usage

**→ [Ant Species Discovery](ANT_SPECIES_DISCOVERY.md)** 🐜  
Automated discovery system for finding ant species with RNA-seq data

**→ [Batched Processing Strategy](discovery/batched_processing.md)** 📦  
Efficient batch processing approach for large-scale workflows

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

## Documentation Files

### Getting Started
- **`MULTI_SPECIES_QUICK_START.md`**: Complete guide for production workflows ⭐
- **`ANT_SPECIES_DISCOVERY.md`**: Comprehensive discovery system documentation 🐜
- **`ANT_DISCOVERY_QUICK_REF.md`**: Quick reference for discovery system
- **`SETUP.md`**: Installation and environment setup
- **`workflow.md`**: Complete workflow planning and execution

### Discovery & Batch Processing
- **`discovery/README.md`**: Discovery results and current status
- **`discovery/batched_processing.md`**: Batch processing strategy and timeline
- Raw discovery data: `output/ant_discovery/*.json`

### Core RNA Analysis
- **`index.md`**: RNA domain overview and module index
- **`configs.md`**: Configuration management and species profiles
- **`steps.md`**: Individual workflow step documentation

### Amalgkit Integration
- **`amalgkit/`**: Amalgkit CLI wrapper documentation
  - **`README.md`**: Overview and quick start
  - **`amalgkit.md`**: Complete transcriptomic analysis pipeline
  - **`comprehensive_guide.md`**: Detailed workflow documentation
  - **`testing_coverage.md`**: Testing, validation, and production results
  - **`R_INSTALLATION.md`**: R dependency installation guide
  - **`steps/`**: Individual step documentation (11 steps)
- **`workflow.md`**: Workflow orchestration with troubleshooting and optimizations

## Related Source Code

- See `src/metainformant/rna/` for implementation details
- See `tests/test_rna_*.py` for comprehensive test coverage
- See `src/metainformant/rna/README.md` for module-specific documentation

## Usage Examples

The RNA domain supports complete transcriptomic workflows:

```python
from metainformant.rna import workflow

# Execute complete RNA analysis workflow
cfg = AmalgkitWorkflowConfig(
    work_dir="output/amalgkit/work",
    species_list=["Homo sapiens"],
    threads=12  # Default thread count for parallel processing
)
results = workflow.execute_workflow(cfg)
```

**Command-line usage** (production ENA-based workflow with auto-activation):
```bash
# Production workflow - no manual venv activation needed
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 12 \
  --threads 12
```

**Legacy SRA-based workflow** (alternative):
```bash
# Legacy workflow using SRA Toolkit
python3 scripts/rna/run_multi_species.py
```

## Integration

RNA analysis integrates with:
- **DNA sequences** for genomic context
- **Protein annotations** for functional analysis
- **Statistical methods** for differential expression
- **Visualization tools** for expression plots

## Testing

Comprehensive tests ensure workflow reliability:
- Workflow execution validation
- Configuration parsing and validation
- External tool integration testing
- Error handling and recovery
- Production validation with 4,548 samples across 20 ant species (November 2025)

Performance characteristics:
- 12 parallel threads for download and quantification
- Batched processing: ~50-100 GB peak disk usage (temporary)
- Direct ENA downloads with 100% reliability (vs 0% SRA Toolkit)
- Virtual environment auto-activation
- wget-based downloads with automatic resume
- Auto-cleanup: FASTQs deleted after quantification

See `amalgkit/testing_coverage.md` for production validation results.

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

## Contributing

When adding new RNA analysis functionality:
1. Update workflow and step documentation
2. Add comprehensive integration tests
3. Update configuration templates
4. Ensure compatibility with amalgkit ecosystem

This documentation provides complete coverage of METAINFORMANT's RNA analysis capabilities.
