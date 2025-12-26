# GWAS Analysis Examples

This directory contains educational examples demonstrating METAINFORMANT's genome-wide association study (GWAS) capabilities, from variant calling to statistical analysis and visualization.

## Overview

These examples showcase the complete GWAS workflow, from genotype data processing through association testing to publication-quality visualizations.

## Examples

### Association Testing (`example_association.py`)

Learn basic GWAS association test implementation and statistical analysis.

**Demonstrates:**
- Linear and logistic regression for association testing
- P-value calculation and significance assessment
- Multiple testing correction (Bonferroni, FDR)
- Effect size estimation

```bash
# Learn GWAS association testing
python examples/gwas/example_association.py
```

**Output:** `output/examples/gwas/association_results.json`

### GWAS Visualization (`example_visualization.py`)

Master GWAS visualization techniques including Manhattan and Q-Q plots.

**Demonstrates:**
- Manhattan plot generation with significance thresholds
- Q-Q plot for distribution checking
- Regional association plots
- Genome-wide visualization

```bash
# Learn GWAS visualization
python examples/gwas/example_visualization.py
```

**Output:** `output/examples/gwas/gwas_plots.json`

## Learning Progression

1. **Statistical Testing**: `example_association.py` - Learn association testing fundamentals
2. **Visualization**: `example_visualization.py` - Master GWAS data visualization

## Related Documentation

- **GWAS Module Docs**: [`docs/gwas/`](../../docs/gwas/) - Complete GWAS analysis documentation
- **Core Examples**: [`examples/core/`](../core/) - Foundational METAINFORMANT concepts
- **Scripts**: [`scripts/gwas/`](../../scripts/gwas/) - Production GWAS analysis workflows

## Data Sources

Examples use simulated genotype and phenotype data for demonstration. For real GWAS analysis, see:
- **SRA Integration**: Raw sequence data downloads from NCBI
- **Variant Calling**: BWA + samtools + bcftools pipeline
- **Production Scripts**: `scripts/gwas/run_genome_scale_gwas.py` for complete workflows
