# RNA Analysis Examples

This directory contains educational examples demonstrating METAINFORMANT's RNA transcriptomic analysis capabilities, focusing on amalgkit workflow integration and expression quantification.

## Overview

These examples showcase RNA-seq analysis workflows, from raw read processing to expression quantification and differential analysis.

## Examples

### Amalgkit Workflow Basics (`example_amalgkit.py`)

Learn amalgkit RNA-seq pipeline integration and workflow planning.

**Demonstrates:**
- Amalgkit workflow configuration and validation
- Step-by-step pipeline planning and execution
- Species-specific workflow management
- Integration with external RNA-seq tools

```bash
# Learn amalgkit workflow basics
python examples/rna/example_amalgkit.py
```

**Output:** `output/examples/rna/amalgkit_workflow.json`

### Expression Quantification (`example_quantification.py`)

Master gene expression quantification and analysis.

**Demonstrates:**
- RNA-seq count data processing
- Normalization techniques (TPM, FPKM)
- Expression level calculations
- Quality assessment metrics

```bash
# Learn expression quantification
python examples/rna/example_quantification.py
```

**Output:** `output/examples/rna/expression_analysis.json`

## Learning Progression

1. **Workflow Setup**: `example_amalgkit.py` - Learn pipeline orchestration
2. **Expression Analysis**: `example_quantification.py` - Master quantification techniques

## Related Documentation

- **RNA Module Docs**: [`docs/rna/`](../../docs/rna/) - Complete RNA analysis documentation
- **Amalgkit Integration**: [`docs/rna/amalgkit/`](../../docs/rna/amalgkit/) - Amalgkit-specific documentation
- **Core Examples**: [`examples/core/`](../core/) - Foundational METAINFORMANT concepts
- **Scripts**: [`scripts/rna/`](../../scripts/rna/) - Production RNA analysis workflows

## Data Sources

Examples use simulated RNA-seq data for demonstration. For real data analysis, see:
- **Amalgkit**: Full RNA-seq pipeline with real data
- **NCBI SRA**: Raw RNA-seq data downloads
- **Production Scripts**: `scripts/rna/run_workflow.py` for complete workflows
