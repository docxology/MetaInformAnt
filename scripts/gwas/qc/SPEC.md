# SPEC: GWAS Quality Control

Technical specification for the GWAS quality control scripts.

> Inherits from [scripts/gwas/SPEC.md](../SPEC.md).

## Purpose

Scripts for filtering variants and samples before association testing.

## Key Filters

| Filter | Threshold |
|--------|-----------|
| Minor Allele Frequency (MAF) | ≥ 1% |
| Genotype Missingness | ≤ 10% |
| Hardy-Weinberg Equilibrium | p ≥ 1e-6 |

## Related Documentation

- **Parent SPEC**: [scripts/gwas/SPEC.md](../SPEC.md)
- **Core Module**: [src/metainformant/gwas/qc/](../../../src/metainformant/gwas/qc/)
- **README**: [README.md](README.md)
- **AGENTS**: [AGENTS.md](AGENTS.md)
