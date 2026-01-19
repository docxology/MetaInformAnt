# SPEC: GWAS Visualization

Technical specification for GWAS result visualization scripts.

> Inherits from [scripts/gwas/SPEC.md](../SPEC.md).

## Purpose

Scripts for generating publication-quality plots from GWAS summary statistics.

## Key Outputs

| Plot Type | Description |
|-----------|-------------|
| Manhattan Plot | Genome-wide view of p-values per locus |
| QQ Plot | Observed vs. expected p-value distribution |
| Regional Plot | Zoomed view of significant loci |

## Related Documentation

- **Parent SPEC**: [scripts/gwas/SPEC.md](../SPEC.md)
- **Visualization Module**: [src/metainformant/visualization/genomics/](../../../src/metainformant/visualization/genomics/)
- **README**: [README.md](README.md)
- **AGENTS**: [AGENTS.md](AGENTS.md)
