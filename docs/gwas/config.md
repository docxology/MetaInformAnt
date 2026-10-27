# GWAS Configuration Reference

Complete reference for GWAS workflow configuration file format.

## Configuration File Format

GWAS configurations support YAML, TOML, or JSON formats. Environment variables can override settings using the `GWAS_` prefix (e.g., `GWAS_THREADS`, `GWAS_WORK_DIR`).

## Top-Level Configuration

### work_dir (required)

Working directory for GWAS workflow execution. All intermediate and output files will be written here or in subdirectories.

```yaml
work_dir: output/gwas/work
```

### log_dir (optional)

Directory for workflow execution logs.

```yaml
log_dir: output/gwas/logs
```

### threads (optional, default: 8)

Number of threads to use for parallel processing.

```yaml
threads: 8
```

## Genome Configuration

### genome (optional)

Reference genome configuration for variant calling or annotation.

```yaml
genome:
  accession: GCF_000001405.40    # NCBI assembly accession
  dest_dir: output/gwas/genome   # Destination directory
  include:                       # Files to download
    - genome
    - gff3
  ftp_url: <optional>           # Direct FTP URL
```

## Variants Configuration

### variants (optional)

Variant data source configuration. Three options available:

#### Option 1: Pre-existing VCF Files

```yaml
variants:
  vcf_files:
    - data/variants/cohort.vcf.gz
    - data/variants/additional.vcf.gz
  dest_dir: output/gwas/variants
```

#### Option 2: Download from Database

```yaml
variants:
  download:
    source: dbSNP              # or: 1000genomes, custom
    accession: GCF_000001405.40
    region: chr1:1000000-2000000  # Optional
  dest_dir: output/gwas/variants
```

#### Option 3: Variant Calling

```yaml
variants:
  calling:
    bam_files:
      - data/alignments/sample1.bam
      - data/alignments/sample2.bam
    reference: output/gwas/genome/genomic.fna
    method: bcftools          # or: gatk, freebayes
  dest_dir: output/gwas/variants
```

## Quality Control Configuration

### qc (optional)

Quality control filter thresholds.

```yaml
qc:
  min_maf: 0.01              # Minimum minor allele frequency
  max_missing: 0.05           # Maximum missing genotype rate
  min_call_rate: 0.95         # Minimum call rate per variant
  hwe_pval: 1e-6             # Hardy-Weinberg equilibrium p-value threshold
  exclude_indels: true        # Exclude insertion/deletion variants
  min_qual: 30.0             # Minimum variant quality (QUAL), inclusive; optional
```

`apply_qc_filters` applies the filters in this order: QUAL (`min_qual`), minor allele frequency (`min_maf`, default 0.01), missing data (`max_missing`, default 0.1) and Hardy-Weinberg equilibrium (`hwe_pval`, passed on as `min_hwe_p`, default 1e-6). A variant has to pass every filter that is enabled. `apply_qc_filters` does not read `min_call_rate` or `exclude_indels`.

`min_qual` is optional: leave it out (or set it to `null`) to switch the QUAL filter off. A variant passes when its QUAL is at least `min_qual`. A variant whose QUAL is missing (`.` in the VCF) cannot be shown to meet the threshold and is removed whatever the threshold value, so omit `min_qual` for a VCF that carries no QUAL values. The result of `apply_qc_filters` then holds a `qual_filter` entry with `min_qual`, `n_removed` and `n_missing_qual` (how many of the removed variants had no QUAL).

Numeric thresholds are read as numbers: `hwe_pval: 1e-6`, which PyYAML loads as a string, is converted to a float, and a value that is not a number raises `ValueError` (`quality_control.<key> must be numeric`) when the workflow normalises the configuration.

## Samples Configuration

### samples (optional)

Sample metadata and phenotype/covariate files.

```yaml
samples:
  sample_list: data/samples.txt        # Optional: list of sample IDs
  phenotype_file: data/phenotypes.tsv  # Required for association testing
  covariates_file: data/covariates.tsv  # Optional: covariates
```

## Population Structure Configuration

### structure (optional)

Population structure analysis settings.

```yaml
structure:
  compute_pca: true           # Compute principal components
  n_components: 10            # Number of PCs to compute
  compute_relatedness: true   # Compute kinship matrix
  kinship_method: vanraden    # or: astle, yang
```

## Association Testing Configuration

### association (required)

Association test model and trait configuration.

```yaml
association:
  model: linear               # or: logistic, mixed
  trait: height               # Column name in phenotype file
  covariates:                 # Covariates to include
    - age
    - sex
    - PC1
    - PC2
    - PC3
  min_sample_size: 50         # Minimum samples per variant
  relatedness_matrix: auto    # or: path to kinship matrix file
```

## Multiple Testing Correction Configuration

### correction (optional)

Multiple testing correction method.

```yaml
correction:
  method: bonferroni          # or: fdr, both
  alpha: 0.05                 # Significance level
```

## Output Configuration

### output (optional)

Output directories and format settings.

```yaml
output:
  results_dir: output/gwas/results  # Association results directory
  plots_dir: output/gwas/plots      # Visualization plots directory
  format: tsv                        # or: json, plink
```

## Environment Variable Overrides

The following environment variables can override configuration file settings:

- `GWAS_THREADS`: Number of threads
- `GWAS_WORK_DIR`: Working directory
- `GWAS_LOG_DIR`: Log directory

Example:

```bash
export GWAS_THREADS=16
export GWAS_WORK_DIR=/path/to/work
python -m metainformant gwas run --config config.yaml
```

## Full Example

See `config/gwas_template.yaml` for a complete example configuration file.

