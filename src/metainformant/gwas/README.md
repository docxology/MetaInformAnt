# GWAS Module

Genome-Wide Association Studies (GWAS) module for identifying genetic variants associated with phenotypic traits.

## Overview

This module provides comprehensive GWAS functionality integrated with METAINFORMANT's bioinformatics ecosystem.

## Components

- **quality.py** - Quality control filters (MAF, missingness, HWE, quality scores)
- **structure.py** - Population structure analysis (PCA, kinship matrices)
- **association.py** - Association testing (linear/logistic regression)
- **correction.py** - Multiple testing correction (Bonferroni, FDR, genomic control)
- **visualization.py** - Result visualization (Manhattan plots, Q-Q plots)
- **calling.py** - Variant calling integration (bcftools, GATK) ✅ Integrated into workflow
- **download.py** - Data acquisition (reference genomes ✅, variant databases ⚠️ placeholder)
- **config.py** - Configuration management and validation
- **workflow.py** - Complete workflow orchestration

## Quick Start

```python
from metainformant.gwas import load_gwas_config, execute_gwas_workflow

# Load and execute workflow
config = load_gwas_config("config/gwas_config.yaml")
results = execute_gwas_workflow(config)
```

## CLI Usage

```bash
python -m metainformant gwas run --config config/gwas_config.yaml
```

## Documentation

Complete documentation available in [`docs/gwas/`](../../docs/gwas/):
- [Overview](../../docs/gwas/index.md)
- [Workflow Guide](../../docs/gwas/workflow.md)
- [Configuration Reference](../../docs/gwas/config.md)

## Integration

This module integrates with:
- `dna.variants` - VCF parsing
- `dna.population` - Population genetics
- `math.popgen` - Theoretical calculations
- `ml.regression` - Regression models
- `visualization.plots` - Plotting utilities

## Scientific Methods

Based on established GWAS methods:
- Price et al. (2006) - PCA for population stratification
- VanRaden (2008) - Kinship matrix computation
- Benjamini & Hochberg (1995) - FDR correction
- Yang et al. (2011) - GCTA methods

## Related Modules

- **DNA**: Variant parsing and population genetics
- **Math**: Theoretical population genetics
- **ML**: Regression and classification models
- **Visualization**: Manhattan and Q-Q plots

