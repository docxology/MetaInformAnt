# SPEC: GWAS Scripts

End-to-end Genome-Wide Association Study (GWAS) orchestrators.

## Pipeline Overview

1. **Quality Control**: Filter variants by MAF, missingness, and HWE.
2. **Population Structure**: PCA and kinship matrix calculation.
3. **Association Testing**: Linear or logistic regression across all loci.
4. **Visualization**: Generation of Manhattan, QQ, and regional plots.

## Best Practices

- **Configuration**: Always use a YAML configuration file to define thresholds and output paths.
- **Performance**: Large-scale associations should leverage the `ParallelProcessor` in `core.parallel`.
