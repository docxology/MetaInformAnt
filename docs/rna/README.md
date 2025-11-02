# RNA Domain Documentation

This directory contains comprehensive documentation for METAINFORMANT's RNA analysis capabilities.

## Quick Start

**→ [Multi-Species Quick Start Guide](MULTI_SPECIES_QUICK_START.md)** ⭐

Complete guide for starting and monitoring multi-species RNA-seq workflows:
- Step-by-step workflow launching
- Real-time progress monitoring
- Troubleshooting and optimization
- Expected timelines and resource usage

## Overview

The RNA domain provides tools for transcriptomic analysis, workflow orchestration, and integration with external bioinformatics tools like amalgkit. Key features include:

- **Direct ENA Downloads**: Robust FASTQ retrieval using ENA API, bypassing problematic SRA Toolkit
- **Retry Logic**: wget-based downloads with --continue for automatic resumption
- **Batched Processing**: Disk-friendly processing of large RNA-seq cohorts (download→quantify→delete)
- **Auto-activation**: Scripts automatically detect and activate virtual environments
- **Multi-species Support**: Coordinated analysis across multiple species
- **Production-Ready**: Currently processing 844 samples across 4 ant species with 100% reliability

## Documentation Files

### Getting Started
- **`MULTI_SPECIES_QUICK_START.md`**: Complete guide for production workflows ⭐
- **`SETUP.md`**: Installation and environment setup
- **`workflow.md`**: Complete workflow planning and execution

### Core RNA Analysis
- **`index.md`**: RNA domain overview and module index
- **`configs.md`**: Configuration management and species profiles
- **`steps.md`**: Individual workflow step documentation

### Amalgkit Integration
- **`amalgkit/`**: Amalgkit CLI wrapper documentation
  - **`amalgkit.md`**: Complete transcriptomic analysis pipeline
  - **`comprehensive_guide.md`**: Detailed workflow documentation
  - **`testing_coverage.md`**: Testing, validation, and production results
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
- Production validation with 844 samples across 4 ant species (November 2025)

Performance characteristics:
- 12 parallel threads for download and quantification
- Batched processing: ~18 GB peak disk usage per batch
- Direct ENA downloads with 100% reliability (vs 0% SRA Toolkit)
- Virtual environment auto-activation
- wget-based downloads with automatic resume

See `amalgkit/testing_coverage.md` for production validation results.

## Contributing

When adding new RNA analysis functionality:
1. Update workflow and step documentation
2. Add comprehensive integration tests
3. Update configuration templates
4. Ensure compatibility with amalgkit ecosystem

This documentation provides complete coverage of METAINFORMANT's RNA analysis capabilities.
