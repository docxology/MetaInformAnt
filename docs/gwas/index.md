# GWAS (Genome-Wide Association Studies)

Comprehensive pipeline for identifying genetic variants associated with phenotypic traits through genome-wide association analysis.

## Overview

The GWAS module provides end-to-end functionality for genome-wide association studies:

- Variant data acquisition (VCF files, database downloads, variant calling)
- Quality control (MAF, missingness, Hardy-Weinberg equilibrium)
- Population structure analysis (PCA, kinship matrices)
- Association testing (linear and logistic regression)
- Multiple testing correction (Bonferroni, FDR, genomic control)
- Visualization (Manhattan plots, Q-Q plots, regional plots)

## Quick Start

```python
from metainformant.gwas import load_gwas_config, execute_gwas_workflow

# Load configuration
config = load_gwas_config("config/gwas_config.yaml")

# Execute workflow
results = execute_gwas_workflow(config)
```

## Command Line

```bash
# Run GWAS from configuration file
python -m metainformant gwas run --config config/gwas_config.yaml

# Validate configuration
python -m metainformant gwas run --config config/gwas_config.yaml --check
```

## Workflow Steps

1. **Genome Preparation** - Download reference genome if needed
2. **Variant Acquisition** - Obtain VCF files (existing, download, or call)
3. **Quality Control** - Filter by MAF, missingness, HWE, quality
4. **Population Structure** - Compute PCA and kinship matrix
5. **Association Testing** - Test variants for trait associations
6. **Correction** - Apply multiple testing correction
7. **Visualization** - Generate Manhattan and Q-Q plots
8. **Export** - Write results tables and summaries

## Configuration

Example YAML configuration:

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

## Module Components

### Quality Control (`quality.py`)
- MAF filtering
- Missing data filtering
- Hardy-Weinberg equilibrium testing
- Quality score filtering
- Indel exclusion

### Population Structure (`structure.py`)
- Principal component analysis (PCA)
- Kinship matrix computation (VanRaden, Astle-Balding, Yang)
- Population stratification correction

### Association Testing (`association.py`)
- Linear regression (continuous traits)
- Logistic regression (binary traits)
- Covariate adjustment
- Effect size estimation

### Multiple Testing Correction (`correction.py`)
- Bonferroni correction
- False discovery rate (FDR/Benjamini-Hochberg)
- Genomic control (lambda_GC calculation)

### Visualization (`visualization.py`)
- Manhattan plots
- Q-Q plots
- Regional association plots

### Variant Calling (`calling.py`)
- bcftools integration
- GATK integration
- Multi-sample variant calling

### Data Download (`download.py`)
- dbSNP integration
- 1000 Genomes Project downloads
- Public GWAS datasets

## Output Structure

```
output/gwas/
├── genome/              # Reference genome
├── variants/            # VCF files (raw and filtered)
├── structure/           # Population structure results
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

## Input File Formats

### Phenotype File (TSV)
```tsv
sample_id    trait_name    covariate1    covariate2
sample1      175.5         45            1
sample2      168.3         38            2
```

### VCF Files
Standard VCF format (version 4.2+), optionally compressed with bgzip and indexed with tabix.

## API Reference

Core functions:
- `load_gwas_config()` - Load workflow configuration
- `execute_gwas_workflow()` - Run complete workflow
- `parse_vcf_full()` - Parse VCF to genotype matrix
- `apply_qc_filters()` - Apply quality control
- `estimate_population_structure()` - Compute PCA/kinship
- `run_gwas()` - Run association tests
- `manhattan_plot()` - Generate Manhattan plot
- `qq_plot()` - Generate Q-Q plot

## Related Modules

- `dna.variants` - VCF parsing utilities
- `dna.population` - Population genetics statistics
- `math.popgen` - Theoretical population genetics
- `ml.regression` - Regression models
- `visualization.plots` - Plotting utilities

## Documentation

- [Complete README](./README.md) - Module overview and features
- [Workflow Guide](./workflow.md) - Step-by-step workflow execution
- [Configuration Reference](./config.md) - Detailed configuration options
- [P. barbatus Example](./pbarbatus_config.md) - Real-world example configuration
- [Verification Report](./verification_report.md) - Validation and testing

## References

1. Price et al. (2006). Principal components analysis corrects for stratification in genome-wide association studies. *Nature Genetics*.
2. VanRaden (2008). Efficient methods to compute genomic predictions. *Journal of Dairy Science*.
3. Benjamini & Hochberg (1995). Controlling the false discovery rate. *Journal of the Royal Statistical Society*.
4. Yang et al. (2011). GCTA: A tool for genome-wide complex trait analysis. *American Journal of Human Genetics*.

## Performance Considerations

- VCF file size impacts memory usage
- Parallel processing available for association testing
- Consider pre-filtering variants by region for large datasets
- Use indexed VCF files for faster random access

## Integration Examples

### With DNA Module
```python
from metainformant.dna.variants import parse_vcf
from metainformant.gwas import apply_qc_filters

# Parse VCF
variants = parse_vcf("data/variants.vcf")

# Apply GWAS-specific QC
filtered = apply_qc_filters(variants, min_maf=0.01)
```

### With Math Module
```python
from metainformant.math.popgen import hardy_weinberg_test
from metainformant.gwas import run_gwas

# Test Hardy-Weinberg equilibrium
hwe_results = hardy_weinberg_test(genotypes)

# Run association with HWE-passing variants
gwas_results = run_gwas(hwe_passing_variants, phenotypes)
```

### With Visualization
```python
from metainformant.gwas import run_gwas
from metainformant.visualization.plots import manhattan_plot

# Run GWAS
results = run_gwas(vcf_path, phenotype_path, config)

# Visualize results
manhattan_plot(results, output_path="plots/manhattan.png")
```

---

This module provides production-ready GWAS functionality integrated with the METAINFORMANT ecosystem for comprehensive genomic analysis.

