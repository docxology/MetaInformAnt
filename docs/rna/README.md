# RNA Domain Documentation

**Master Index** - Complete documentation for METAINFORMANT's RNA analysis capabilities.

## Quick Start

Choose your path:

### 🆕 New User
1. **[Getting Started Guide](GETTING_STARTED.md)** - Installation and setup
2. **[Workflow Guide](workflow.md)** - Understanding RNA-seq workflows
3. **[Step Documentation](amalgkit/steps/README.md)** - Learn individual steps

### 🚀 Production Workflows
1. **[Multi-Species Quick Start](MULTI_SPECIES_QUICK_START.md)** - Production workflow guide
2. **[Orchestration Guide](ORCHESTRATION/README.md)** - Choose the right orchestrator
3. **[Configuration Guide](CONFIGURATION.md)** - Configure workflows

### 📚 API Reference
1. **[API Reference](API.md)** - Complete function documentation
2. **[Function Index](amalgkit/FUNCTIONS.md)** - Quick function lookup
3. **[Step Documentation](amalgkit/steps/README.md)** - Detailed step guides

---

## Documentation Map

### By Use Case

| Use Case | Start Here | Next Steps |
|----------|------------|------------|
| **First-time user** | [GETTING_STARTED.md](GETTING_STARTED.md) | [workflow.md](workflow.md) → [amalgkit/steps/](amalgkit/steps/) |
| **Production analysis** | [MULTI_SPECIES_QUICK_START.md](MULTI_SPECIES_QUICK_START.md) | [ORCHESTRATION/README.md](ORCHESTRATION/README.md) → [CONFIGURATION.md](CONFIGURATION.md) |
| **Species discovery** | [discovery/GUIDE.md](discovery/GUIDE.md) | [discovery/QUICK_REF.md](discovery/QUICK_REF.md) |
| **API development** | [API.md](API.md) | [amalgkit/FUNCTIONS.md](amalgkit/FUNCTIONS.md) → [amalgkit/steps/](amalgkit/steps/) |
| **Genome setup** | [amalgkit/genome_setup_guide.md](amalgkit/genome_setup_guide.md) | [amalgkit/genome_preparation.md](amalgkit/genome_preparation.md) |

### By Component

#### Core Workflow
- **[workflow.md](workflow.md)** - Workflow planning and execution
- **[steps.md](steps.md)** - Step overview and integration
- **[CONFIGURATION.md](CONFIGURATION.md)** - Configuration management
- **[index.md](index.md)** - Module overview

#### Amalgkit Integration
- **[amalgkit/README.md](amalgkit/README.md)** - Amalgkit wrapper overview
- **[amalgkit/amalgkit.md](amalgkit/amalgkit.md)** - Complete pipeline documentation
- **[amalgkit/steps/](amalgkit/steps/)** - All 11 step guides
- **[amalgkit/FUNCTIONS.md](amalgkit/FUNCTIONS.md)** - Function quick reference
- **[amalgkit/quick_start.md](amalgkit/quick_start.md)** - Quick start guide

#### Genome Setup
- **[amalgkit/genome_setup_guide.md](amalgkit/genome_setup_guide.md)** - User guide (step-by-step)
- **[amalgkit/genome_preparation.md](amalgkit/genome_preparation.md)** - Technical API reference
- **[amalgkit/commands.md](amalgkit/commands.md)** - Command reference

#### Orchestration
- **[ORCHESTRATION/README.md](ORCHESTRATION/README.md)** - Orchestrator overview
- **[ORCHESTRATION/ENA_WORKFLOW.md](ORCHESTRATION/ENA_WORKFLOW.md)** - ENA-based workflow
- **[ORCHESTRATION/MULTI_SPECIES.md](ORCHESTRATION/MULTI_SPECIES.md)** - Multi-species workflow
- **[ORCHESTRATION/BATCH_DOWNLOAD.md](ORCHESTRATION/BATCH_DOWNLOAD.md)** - Batch download

#### Discovery
- **[discovery/GUIDE.md](discovery/GUIDE.md)** - Comprehensive guide
- **[discovery/QUICK_REF.md](discovery/QUICK_REF.md)** - Quick reference
- **[discovery/README.md](discovery/README.md)** - Current status

#### Examples
- **[examples/README.md](examples/README.md)** - Example index
- **[examples/pbarbatus_analysis.md](examples/pbarbatus_analysis.md)** - Complete analysis example

---

## Core Concepts

### Workflow Execution

The RNA module provides a complete workflow system for transcriptomic analysis:

1. **Planning**: Use `plan_workflow()` to generate step sequence
2. **Execution**: Use `execute_workflow()` to run complete pipelines
3. **Configuration**: YAML configs with environment variable overrides

**Key Documents**:
- [workflow.md](workflow.md) - Workflow planning and execution
- [CONFIGURATION.md](CONFIGURATION.md) - Configuration management
- [API.md](API.md#workflow-functions) - Workflow function reference

### Amalgkit Steps

11 workflow steps from metadata retrieval to quality assurance:

1. **metadata** - Retrieve RNA-seq metadata from NCBI SRA/ENA
2. **config** - Generate configuration files
3. **select** - Filter samples by quality criteria
4. **getfastq** - Download and convert SRA to FASTQ
5. **integrate** - Integrate FASTQ paths into metadata
6. **quant** - Quantify transcript abundances (kallisto/salmon)
7. **merge** - Merge per-sample results into expression matrices
8. **cstmm** - Cross-species TMM normalization
9. **curate** - Quality control and batch correction
10. **csca** - Cross-species correlation analysis
11. **sanity** - Validate workflow outputs

**Key Documents**:
- [amalgkit/steps/README.md](amalgkit/steps/README.md) - Step index with workflow diagram
- [amalgkit/steps/](amalgkit/steps/) - Individual step documentation
- [API.md](API.md#amalgkit-step-functions) - Step function reference

### Genome Preparation

Automated genome download, transcriptome extraction, and kallisto index building:

**Key Documents**:
- [amalgkit/genome_setup_guide.md](amalgkit/genome_setup_guide.md) - User guide
- [amalgkit/genome_preparation.md](amalgkit/genome_preparation.md) - Technical API
- [API.md](API.md#genome-preparation-functions) - Function reference

---

## API Reference

### Quick Access

- **[API.md](API.md)** - Complete function reference with signatures
- **[amalgkit/FUNCTIONS.md](amalgkit/FUNCTIONS.md)** - Quick function lookup table

### Function Categories

| Category | Functions | Documentation |
|----------|-----------|---------------|
| **Amalgkit Steps** | `metadata`, `quant`, `merge`, etc. | [API.md](API.md#amalgkit-step-functions) |
| **Step Runners** | `run_metadata`, `run_quant`, etc. | [API.md](API.md#step-runner-functions) |
| **Workflow** | `plan_workflow`, `execute_workflow` | [API.md](API.md#workflow-functions) |
| **Genome Prep** | `prepare_genome_for_quantification`, etc. | [API.md](API.md#genome-preparation-functions) |
| **Orchestration** | `discover_species_configs`, etc. | [API.md](API.md#orchestration-functions) |
| **Utilities** | `check_cli_available`, `build_cli_args` | [API.md](API.md#utility-functions) |

---

## Quick Commands

### Single Species Workflow

```bash
# ENA-based workflow (recommended for production)
python3 scripts/rna/workflow_ena_integrated.py \
  --config config/amalgkit/amalgkit_cfloridanus.yaml \
  --batch-size 12 --threads 12
```

### Multi-Species Workflow

```bash
# Immediate processing with 24 threads TOTAL distributed across all species
python3 scripts/rna/batch_download_species.py --total-threads 24

# Or sequential execution (one species at a time)
export AK_THREADS=24
python3 scripts/rna/run_multi_species.py
```

### Genome Setup

```bash
# Check status
python3 scripts/rna/verify_genomes_and_indexes.py

# Complete setup (download → prepare → index)
python3 scripts/rna/orchestrate_genome_setup.py

# Or step-by-step:
python3 scripts/rna/download_missing_genomes.py
python3 scripts/rna/prepare_transcriptomes.py
python3 scripts/rna/build_kallisto_indexes.py
```

### Species Discovery

```bash
# Discover ant species with RNA-seq data
python3 scripts/rna/discover_ant_rnaseq_by_genus.py

# Generate configs for discovered species
python3 scripts/rna/generate_ant_configs_with_genomes.py
```

---

## Advanced Topics

### Species Discovery

Automated discovery of species with RNA-seq data and genome assembly retrieval:

- **[discovery/GUIDE.md](discovery/GUIDE.md)** - Comprehensive guide
- **[discovery/QUICK_REF.md](discovery/QUICK_REF.md)** - Quick reference
- **[discovery/README.md](discovery/README.md)** - Current status and results

### Orchestration

Multiple orchestrators for different use cases:

- **[ORCHESTRATION/README.md](ORCHESTRATION/README.md)** - Choose the right orchestrator
- **[ORCHESTRATION/ENA_WORKFLOW.md](ORCHESTRATION/ENA_WORKFLOW.md)** - ENA-based (recommended)
- **[ORCHESTRATION/MULTI_SPECIES.md](ORCHESTRATION/MULTI_SPECIES.md)** - Multi-species SRA workflow
- **[ORCHESTRATION/BATCH_DOWNLOAD.md](ORCHESTRATION/BATCH_DOWNLOAD.md)** - Configurable batch downloads

### Examples

Real-world analysis examples:

- **[examples/pbarbatus_analysis.md](examples/pbarbatus_analysis.md)** - Complete Pogonomyrmex barbatus analysis
- **[examples/pbarbatus_quick_reference.md](examples/pbarbatus_quick_reference.md)** - Quick reference

---

## Troubleshooting

| Issue | Documentation |
|-------|---------------|
| **Setup problems** | [GETTING_STARTED.md](GETTING_STARTED.md#troubleshooting) |
| **Workflow errors** | [workflow.md](workflow.md#troubleshooting) |
| **Step-specific issues** | [amalgkit/steps/](amalgkit/steps/) - See individual step docs |
| **Configuration** | [CONFIGURATION.md](CONFIGURATION.md#troubleshooting) |
| **Performance** | [MULTI_SPECIES_QUICK_START.md](MULTI_SPECIES_QUICK_START.md#performance-characteristics) |

---

## Testing

- **Test Coverage**: [amalgkit/testing_coverage.md](amalgkit/testing_coverage.md)
- **Running Tests**: `pytest tests/test_rna_*.py -v`
- **End-to-End Validation**: [END_TO_END_VALIDATION.md](END_TO_END_VALIDATION.md)

---

## Current Status (November 2025)

- **Ant Species Discovery**: 55 species discovered, 20 validated genomes
- **Batch 1**: 10 species (3,820 samples) processing
- **Batch 2**: 10 species (728 samples) queued
- **Total Pipeline**: 4,548 samples across 20 ant species

See [discovery/README.md](discovery/README.md) for complete details.

---

## Module Overview

The RNA domain provides:

- **Complete Workflow System**: From metadata retrieval to expression matrices
- **Species Discovery**: Automated identification of species with RNA-seq data
- **Genome Integration**: Validated genome assemblies from NCBI
- **Direct ENA Downloads**: Robust FASTQ retrieval (100% reliability)
- **Multi-Species Support**: Coordinated analysis across multiple species
- **Production-Ready**: Tested with 4,548+ samples

**Source Code**: `src/metainformant/rna/`  
**Module Documentation**: `src/metainformant/rna/README.md`  
**Tests**: `tests/test_rna_*.py`

---

## See Also

- **[index.md](index.md)** - RNA domain overview
- **[API.md](API.md)** - Complete API reference
- **[amalgkit/FUNCTIONS.md](amalgkit/FUNCTIONS.md)** - Function quick reference
- **[workflow.md](workflow.md)** - Workflow planning and execution
- **[CONFIGURATION.md](CONFIGURATION.md)** - Configuration management

---

**Last Updated**: November 2025  
**Documentation Status**: Complete and organized  
**Production Status**: Validated with 4,548 samples across 20 ant species
