# SPEC: GWAS Association

Technical specification for association testing orchestration.

> Inherits from [scripts/gwas/SPEC.md](../SPEC.md).

## Purpose

This subdirectory contains scripts that orchestrate the core statistical association testing phase of a GWAS workflow.

## Key Scripts

| Script | Description |
|--------|-------------|
| `run_association.py` | Runs linear/logistic regression across all variants |

## Workflow

1. Load filtered genotype matrix and phenotype vector.
2. Apply configured covariates (e.g., PCs from population structure).
3. Perform single-variant association testing.
4. Output summary statistics with p-values and effect sizes.

## Related Documentation

- **Parent SPEC**: [scripts/gwas/SPEC.md](../SPEC.md)
- **Core Module**: [src/metainformant/gwas/SPEC.md](../../../src/metainformant/gwas/SPEC.md)
- **README**: [README.md](README.md)
- **AGENTS**: [AGENTS.md](AGENTS.md)
