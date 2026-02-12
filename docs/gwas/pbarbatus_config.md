# GWAS Configuration for Pogonomyrmex barbatus

This document describes the GWAS configuration setup for *Pogonomyrmex barbatus* (Red Harvester Ant).

## Configuration File

**Location**: `config/gwas/gwas_pbarbatus.yaml`

## Species Information

- **Species**: *Pogonomyrmex barbatus* (Red Harvester Ant)
- **NCBI Taxonomy ID**: 144034
- **Reference Assembly**: GCF_000187915.1 (Pbar_UMD_V03)
- **FTP URL**: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/187/915/GCF_000187915.1_Pbar_UMD_V03/

## Configuration Overview

The configuration file follows the GWAS template structure and includes:

### 1. Reference Genome

```yaml
genome:
  accession: GCF_000187915.1
  dest_dir: output/gwas/pbarbatus/genome
  include:
    - genome         # Genomic sequences (FASTA)
    - gff3          # Gene annotations
  ftp_url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/187/915/GCF_000187915.1_Pbar_UMD_V03/
```

### 2. Variant Data Sources

The configuration supports three options for variant data:
- **Pre-called VCF files**: Specify paths in `variants.vcf_files`
- **Download from database**: Configure `variants.download` section
- **Variant calling**: Use `variants.calling` section with BAM/CRAM files

### 3. Quality Control

Default QC parameters:
- MAF threshold: 0.01
- Max missing: 0.05
- Min call rate: 0.95
- HWE p-value: 1e-6
- Exclude indels: true
- Min quality: 30.0

### 4. Population Structure

- PCA computation enabled
- 10 principal components
- Kinship matrix computation (VanRaden method)

### 5. Association Testing

- Model: Linear regression
- Trait: `trait1` (specify in phenotype file)
- Min sample size: 50

### 6. Multiple Testing Correction

- Method: Bonferroni
- Alpha: 0.05

## Usage

### Configuration Validation

```bash
python -m metainformant gwas run --config config/gwas/gwas_pbarbatus.yaml --check
```

### Full Workflow Execution

```bash
python -m metainformant gwas run --config config/gwas/gwas_pbarbatus.yaml
```

### Python API

```python
from metainformant.gwas import load_gwas_config, execute_gwas_workflow

# Load configuration
cfg = load_gwas_config("config/gwas_pbarbatus.yaml")

# Validate configuration
result = execute_gwas_workflow(cfg, check=True)

# Execute workflow
result = execute_gwas_workflow(cfg, check=False)
```

## Required Data Files

Before running the full workflow, ensure you have:

1. **Phenotype file** (`data/phenotypes/pbarbatus/phenotypes.tsv`):
   - TSV format with `sample_id` and trait columns
   - Example:
     ```
     sample_id	trait1	trait2
     S1	10.5	15.2
     S2	11.0	16.0
     ```

2. **Variant data** (one of):
   - VCF file(s) specified in `variants.vcf_files`
   - Or BAM/CRAM files for variant calling in `variants.calling.bam_files`

3. **Covariates file** (optional):
   - TSV format with `sample_id` and covariate columns

## Workflow Steps

1. **Genome Download**: Downloads P. barbatus reference genome from NCBI
2. **Variant Acquisition**: Uses pre-existing VCF, downloads, or calls variants
3. **Quality Control**: Filters variants based on QC parameters
4. **Population Structure**: Computes PCA and kinship matrices
5. **Association Testing**: Performs GWAS with specified model
6. **Multiple Testing Correction**: Applies correction methods
7. **Visualization**: Generates Manhattan plots, Q-Q plots, etc.
8. **Results Export**: Writes results to output directories

## Output Directories

- **Work directory**: `output/gwas/pbarbatus/work`
- **Genome**: `output/gwas/pbarbatus/genome`
- **Variants**: `output/gwas/pbarbatus/variants`
- **Results**: `output/gwas/pbarbatus/results`
- **Plots**: `output/gwas/pbarbatus/plots`
- **Logs**: `output/gwas/pbarbatus/logs`

## Testing

Run configuration tests:

```bash
pytest tests/test_gwas_config_pbarbatus.py -v
```

## Related Documentation

- [GWAS Configuration Template](../../config/gwas/gwas_template.yaml)
- [GWAS Workflow Documentation](./workflow.md)
- [GWAS Configuration Guide](./config.md)

