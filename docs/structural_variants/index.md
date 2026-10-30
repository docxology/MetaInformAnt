# Structural Variants Module

Structural variant (SV) and copy number variation (CNV) detection, annotation, and visualization.

## Submodules

### Detection (`detection/`)
- **cnv.py** - CNV detection from depth, coverage segmentation, CNV state calling, log2 ratio calculation
- **sv_calling.py** - SV calling, split-read/discordant-pair detection, SV typing, genotyping
- **breakpoints.py** - Breakpoint refinement, microhomology detection, clustering, confidence scoring

### Annotation (`annotation/`)
- **overlap.py** - Gene/regulatory overlap annotation, overlap fraction, nearest gene finding
- **functional_impact.py** - Functional impact prediction, dosage sensitivity, TAD disruption, pathogenicity scoring

### Filtering (`filtering/`)
- **quality_filter.py** - Quality/size/frequency filtering, blacklist application
- **merge.py** - Multi-caller merging, reciprocal overlap, SURVIVOR merge, deduplication

### Population (`population/`)
- **sv_population.py** - Cohort genotyping, allele frequency, association testing, PCA population structure, SV-SNP linkage disequilibrium, multi-sample callset merging

### Visualization (`visualization/`)
- **plots.py** - Circos plots, coverage tracks, SV size distributions, type summaries, breakpoint details, CNV profiles

## Usage

```python
from metainformant.structural_variants.detection import cnv, sv_calling, breakpoints
from metainformant.structural_variants.annotation import overlap, functional_impact
from metainformant.structural_variants.filtering import quality_filter, merge
from metainformant.structural_variants.population import sv_population
```

See [detection](detection.md), [annotation](annotation.md), [filtering](filtering.md) and [population](population.md) for the function reference of each submodule, and [SPEC.md](SPEC.md) for the data model.

Population PCA uses the same sample-score convention for either matrix shape.
Numerically undefined components have zero scores, loadings and eigenvalues;
repeated genotype rows do not create extra population dimensions.

Native pysam VCF input permits absent optional `CHR2`/`POS2` metadata when BND
ALT supplies the mate. ALT mates take precedence, so unused optional fields are not queried.
Errors in required native fallback metadata still propagate.
