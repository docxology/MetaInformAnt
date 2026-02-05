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

### Visualization (`visualization/`)
- **plots.py** - Circos plots, coverage tracks, SV size distributions, type summaries, breakpoint details, CNV profiles

## Usage

```python
from metainformant.structural_variants.detection import cnv, sv_calling, breakpoints
from metainformant.structural_variants.annotation import overlap, functional_impact
from metainformant.structural_variants.filtering import quality_filter, merge
```
