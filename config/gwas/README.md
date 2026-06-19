# GWAS Configurations

Configuration files for GWAS workflows.

## Apis mellifera Profiles

- The real BeeWAS 2026 read-cohort profile lives in the nested project at
  `projects/apis_gwas/config/beewas_2026_full_guarded.yaml`. It points at
  `/Volumes/blue/data/beewas_2026`, the Amel_HAv3.1 reference, the clean-main
  FASTQ manifest, CRAM/VCF targets, curated phenotype workbook outputs,
  genomic-estimator readiness artifacts, and the reviewed phenotype gate.
- `gwas_amellifera.yaml` is the generated/demo Apis profile used by
  `scripts/gwas/run_amellifera_gwas.py`. It contains `data_generation` and
  should not be interpreted as real cohort evidence.

## Other Profiles

- `gwas_pbarbatus.yaml` - Pogonomyrmex barbatus profile
- `gwas_pbarbatus_synthetic.yaml` - synthetic P. barbatus test profile
- `gwas_template.yaml` - generic template
