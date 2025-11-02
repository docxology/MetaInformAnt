# GWAS (Genome-Wide Association Studies) Module

Comprehensive end-to-end GWAS pipeline for identifying genetic variants associated with phenotypic traits.

## Overview

The GWAS module provides a complete workflow for genome-wide association studies, from variant data acquisition through quality control, population structure analysis, association testing, multiple testing correction, and visualization.

## Features

- **Variant Data Acquisition**: Support for pre-existing VCF files, variant calling from BAM/CRAM files (bcftools, GATK), and reference genome downloads. Note: Direct database downloads (dbSNP, 1000 Genomes) are placeholders for future implementation.
- **Quality Control**: MAF filtering, missing data filtering, Hardy-Weinberg equilibrium testing, indel exclusion
- **Population Structure**: Principal component analysis (PCA) and kinship matrix computation
- **Association Testing**: Linear and logistic regression models with covariate adjustment
- **Multiple Testing Correction**: Bonferroni, FDR (Benjamini-Hochberg), and genomic control
- **Visualization**: Manhattan plots, Q-Q plots, and regional association plots

## Quick Start

### Basic Usage

```python
from metainformant.gwas import load_gwas_config, execute_gwas_workflow

# Load configuration
config = load_gwas_config("config/gwas_template.yaml")

# Run workflow
results = execute_gwas_workflow(config)
```

### Command Line

```bash
# Run GWAS workflow from config file
python -m metainformant gwas run --config config/gwas_config.yaml

# Validate configuration only
python -m metainformant gwas run --config config/gwas_config.yaml --check
```

## Workflow Steps

1. **Genome Preparation**: Download reference genome (if needed)**
2. **Variant Acquisition**: Obtain variant data (VCF files, download, or calling)
3. **Quality Control**: Filter variants by MAF, missingness, HWE, quality
4. **Population Structure**: Compute PCA and kinship matrix
5. **Association Testing**: Test each variant for association with trait
6. **Multiple Testing Correction**: Apply correction for multiple comparisons
7. **Visualization**: Generate Manhattan plots, Q-Q plots
8. **Results Export**: Write results tables and summary reports

## Configuration

See [Configuration Reference](./config.md) for detailed configuration options.

### Example Configuration

```yaml
work_dir: output/gwas/work
threads: 8

genome:
  accession: GCF_000001405.40
  dest_dir: output/gwas/genome

variants:
  vcf_files:
    - data/variants/cohort.vcf.gz

qc:
  min_maf: 0.01
  max_missing: 0.05
  hwe_pval: 1e-6

samples:
  phenotype_file: data/phenotypes.tsv
  covariates_file: data/covariates.tsv

association:
  model: linear
  trait: height
  covariates: [age, sex, PC1, PC2]

correction:
  method: bonferroni
  alpha: 0.05
```

## Input Files

### Phenotype File (TSV)

```tsv
sample_id    height    weight    disease_status
sample1      175.5     70.2      0
sample2      168.3     65.1      1
...
```

### Covariates File (TSV)

```tsv
sample_id    age    sex    PC1        PC2        PC3
sample1      45     1      0.0234     -0.0123    0.0045
sample2      38     2      -0.0156    0.0089     -0.0023
...
```

## Output Structure

```
output/gwas/
├── genome/              # Downloaded reference genome
├── variants/            # VCF files (raw and QC'd)
├── structure/           # PCA results, kinship matrices
│   ├── pca_components.tsv
│   ├── kinship_matrix.tsv
│   └── structure_summary.json
├── results/             # Association test results
│   ├── association_results.tsv
│   └── summary.json
└── plots/               # Visualizations
    ├── manhattan.png
    ├── qq_plot.png
    └── regional_*.png
```

## API Reference

### Core Functions

- `load_gwas_config(config_path)`: Load GWAS configuration from file
- `execute_gwas_workflow(config)`: Run complete GWAS workflow
- `parse_vcf_full(vcf_path)`: Parse VCF and extract genotype matrix
- `apply_qc_filters(vcf_path, config)`: Apply quality control filters
- `estimate_population_structure(vcf_path, config)`: Compute PCA and kinship
- `run_gwas(vcf_path, phenotype_path, config)`: Run association tests
- `manhattan_plot(results, output_path)`: Generate Manhattan plot
- `qq_plot(pvalues, output_path)`: Generate Q-Q plot

See module documentation for full API details.

## Examples

See [Examples](./examples/) directory for complete usage examples.

## Related Modules

- `dna.variants`: Basic VCF parsing
- `dna.population`: Population genetics statistics
- `math.popgen`: Population genetics theory
- `ml.regression`: Regression models
- `visualization.plots`: Plotting utilities

## Limitations

See [Comprehensive Review](./comprehensive_review.md) for detailed assessment of functionality and gaps.

**Current Status**:
- ✅ Core GWAS workflow: Fully functional
- ✅ Variant calling: Integrated and working
- ⚠️ Database downloads: Placeholder (use pre-existing VCF or calling)
- ❌ Functional annotation: Not yet implemented
- ⚠️ Advanced methods: Basic implementation (MLM planned)

## References

1. Price et al. (2006). Principal components analysis corrects for stratification in genome-wide association studies. *Nature Genetics*.
2. VanRaden (2008). Efficient methods to compute genomic predictions. *Journal of Dairy Science*.
3. Benjamini & Hochberg (1995). Controlling the false discovery rate. *Journal of the Royal Statistical Society*.

