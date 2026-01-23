# Agent Directives: config/gwas

## Role
GWAS pipeline configuration files.

## Contents
- `gwas_template.yaml` - Full template with all options
- `gwas_amellifera.yaml` - Apis mellifera GWAS config
- `gwas_pbarbatus.yaml` - Pogonomyrmex barbatus config
- `gwas_pbarbatus_synthetic.yaml` - Synthetic data for testing

## Configuration Structure
```yaml
work_dir: output/gwas/{species}
threads: 8
reference:
  accession: "GCF_..."
  fasta: null  # auto-downloaded
quality_control:
  min_maf: 0.01
  max_missing: 0.1
  min_hwe_p: 1e-6
association:
  method: linear  # or logistic
  covariates: []
```

## Environment Overrides
Use `GWAS_` prefix:
- `GWAS_THREADS=16`
- `GWAS_WORK_DIR=/path/to/output`
