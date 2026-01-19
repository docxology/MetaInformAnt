# SPEC: GWAS Pipelines

Technical specification for end-to-end GWAS pipeline orchestration.

> Inherits from [scripts/gwas/SPEC.md](../SPEC.md).

## Purpose

Scripts that orchestrate the full GWAS workflow from raw data to final results.

## Pipeline Stages

1. **Preparation**: Load genotypes, phenotypes, covariates.
2. **QC**: Filter variants and samples.
3. **Structure**: Calculate PCA and kinship.
4. **Association**: Run statistical tests.
5. **Visualization**: Generate plots.

## Related Documentation

- **Parent SPEC**: [scripts/gwas/SPEC.md](../SPEC.md)
- **Subdirectory SPECs**: [association/](./association/SPEC.md), [qc/](./qc/SPEC.md), [structure/](./structure/SPEC.md), [visualization/](./visualization/SPEC.md), [preparation/](./preparation/SPEC.md)
- **README**: [README.md](README.md)
- **AGENTS**: [AGENTS.md](AGENTS.md)
