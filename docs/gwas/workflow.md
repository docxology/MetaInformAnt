# GWAS Workflow Guide

Step-by-step guide to running a GWAS workflow.

## Workflow Overview

The GWAS workflow executes the following steps in order:

1. **Configuration Loading** - Load and validate configuration file
2. **Genome Preparation** - Download reference genome if needed
3. **Variant Acquisition** - Obtain variant data (VCF files, download, or calling)
4. **Quality Control** - Apply QC filters (MAF, missingness, HWE)
5. **Population Structure** - Compute PCA and kinship matrix
6. **Association Testing** - Run association tests for each variant
7. **Multiple Testing Correction** - Apply correction methods
8. **Visualization** - Generate Manhattan plots, Q-Q plots
9. **Results Export** - Write results tables and summary reports

## Prerequisites

- Python 3.11+
- NumPy
- Matplotlib (for visualization)
- Optional: scipy (for advanced statistics)
- Optional: pandas (for data handling)
- Optional: bcftools (for variant calling/extraction)
- Optional: GATK (for variant calling)

## Step-by-Step Workflow

### Step 1: Prepare Configuration

Create a configuration file (YAML recommended):

```yaml
work_dir: output/gwas/my_study
threads: 8

genome:
  accession: GCF_000001405.40
  dest_dir: output/gwas/genome

variants:
  vcf_files:
    - data/variants/my_cohort.vcf.gz

qc:
  min_maf: 0.01
  max_missing: 0.05

samples:
  phenotype_file: data/phenotypes.tsv

association:
  model: linear
  trait: height
  covariates: [age, sex]
```

### Step 2: Prepare Input Files

#### Phenotype File Format

Create a TSV file with sample IDs and trait values:

```tsv
sample_id    height    weight
sample1      175.5     70.2
sample2      168.3     65.1
sample3      180.1     75.8
```

#### Covariates File Format (Optional)

```tsv
sample_id    age    sex
sample1      45     1
sample2      38     2
sample3      52     1
```

### Step 3: Run Workflow

#### Command Line

```bash
python -m metainformant gwas run --config config/my_gwas_config.yaml
```

#### Python API

```python
from metainformant.gwas import load_gwas_config, execute_gwas_workflow

config = load_gwas_config("config/my_gwas_config.yaml")
results = execute_gwas_workflow(config)
```

### Step 4: Review Results

Results are written to `output/gwas/my_study/`:

- **Results**: `results/association_results.tsv` - Association test results
- **Plots**: `plots/manhattan.png`, `plots/qq_plot.png` - Visualizations
- **Structure**: `structure/pca_components.tsv`, `structure/kinship_matrix.tsv` - Population structure

## Workflow Steps Explained

### Genome Preparation

If genome configuration is provided, the reference genome is downloaded using NCBI infrastructure (reusing `dna.ncbi` module).

### Variant Acquisition

Three options available:

1. **Pre-existing VCF**: Specify paths to existing VCF files ✅ Fully supported
2. **Download**: Download from public databases (dbSNP, 1000 Genomes) ⚠️ Placeholder (requires external tools)
3. **Calling**: Call variants from BAM/CRAM files using bcftools or GATK ✅ Fully supported

### Quality Control

Applied filters:
- **MAF**: Remove variants with minor allele frequency < threshold
- **Missing Data**: Remove variants with > max_missing missing genotypes
- **Hardy-Weinberg**: Remove variants deviating from HWE (p < threshold)
- **Quality**: Remove variants with quality score < min_qual
- **Indels**: Optionally exclude insertion/deletion variants

### Population Structure

- **PCA**: Principal component analysis on genotype matrix
- **Kinship**: Compute kinship/relatedness matrix (VanRaden, Astle-Balding, or Yang method)

Principal components can be included as covariates in association testing to control for population stratification.

### Association Testing

Per-variant association tests:
- **Linear Regression**: For continuous traits
- **Logistic Regression**: For binary traits (case/control)

Model: `phenotype ~ genotype + covariates`

Outputs: beta (effect size), standard error, p-value, R² (for linear models), odds ratio (for logistic models)

### Multiple Testing Correction

- **Bonferroni**: Conservative correction (alpha / n_tests)
- **FDR (Benjamini-Hochberg)**: False discovery rate control
- **Genomic Control**: Calculate genomic inflation factor (lambda_GC)

### Visualization

- **Manhattan Plot**: -log10(p-value) vs. genomic position, colored by chromosome
- **Q-Q Plot**: Observed vs. expected p-values under null, shows genomic inflation
- **Regional Plot**: Zoomed view of specific genomic regions

## Troubleshooting

### Common Issues

1. **VCF file not found**: Ensure VCF paths in configuration are correct
2. **Phenotype file missing samples**: Check that sample IDs match between VCF and phenotype files
3. **Memory errors**: Reduce number of variants or samples, or use subset
4. **bcftools/GATK not found**: Install required tools or use pre-existing VCF files

### Performance Tips

- Use pre-filtered VCF files when possible
- Limit variant set to region of interest
- Use subset of samples for initial testing
- Run QC filters before association testing to reduce computation

## Advanced Usage

### Custom Association Models

You can use the low-level API to run custom association tests:

```python
from metainformant.gwas import association_test_linear, parse_vcf_full

vcf_data = parse_vcf_full("variants.vcf")
for variant_idx in range(len(vcf_data["variants"])):
    genotypes = [vcf_data["genotypes"][s][variant_idx] for s in range(len(vcf_data["samples"]))]
    result = association_test_linear(genotypes, phenotypes, covariates)
```

### Population Structure Analysis

```python
from metainformant.gwas import estimate_population_structure

structure = estimate_population_structure(
    vcf_path="variants.vcf",
    config={"compute_pca": True, "n_components": 10},
    output_dir="output/structure"
)
```

### Visualization

```python
from metainformant.gwas import manhattan_plot, qq_plot

manhattan_plot(
    results="results/association_results.tsv",
    output_path="plots/manhattan.png",
    significance_threshold=5e-8
)

qq_plot(
    pvalues="results/association_results.tsv",
    output_path="plots/qq.png"
)
```

