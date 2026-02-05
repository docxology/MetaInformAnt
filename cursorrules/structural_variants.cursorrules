# Structural Variants Module Rules

## Purpose
Structural variant (SV) and copy number variation (CNV) detection, annotation, filtering, merging, and visualization including circos plots.

## Source Structure
```
src/metainformant/structural_variants/
├── detection/
│   ├── cnv.py          # CNV detection from depth, segmentation, state calling
│   ├── sv_calling.py   # SV calling, split-read/discordant-pair, genotyping
│   └── breakpoints.py  # Breakpoint refinement, microhomology, clustering
├── annotation/
│   ├── overlap.py         # Gene/regulatory overlap, overlap fraction, nearest gene
│   └── functional_impact.py # Functional impact, dosage sensitivity, TAD disruption
├── filtering/
│   ├── quality_filter.py  # Quality/size/frequency filtering, blacklists
│   └── merge.py           # Multi-caller merging, reciprocal overlap, SURVIVOR
└── visualization/
    └── plots.py           # Circos plots, coverage tracks, SV size distributions
```

## Dependencies
- **Required**: numpy, pandas, pysam
- **Optional**: samtools, SURVIVOR (multi-caller merge), bedtools (interval operations)

## Import Patterns
```python
from metainformant.structural_variants.detection import cnv, sv_calling, breakpoints
from metainformant.structural_variants.annotation import overlap, functional_impact
from metainformant.structural_variants.filtering import quality_filter, merge
```

## Configuration
- Environment prefix: `SV_` (e.g., `SV_THREADS`, `SV_WORK_DIR`)
- Output path: `output/structural_variants/<analysis_type>/` (e.g., `detection/`, `annotation/`, `filtering/`)

## Integration
- **Structural Variants → DNA**: Genomic coordinates, variant calling
- **Structural Variants → Longread**: SV detection from long reads
- **Structural Variants → GWAS**: Structural variants in association studies
- **Structural Variants → Visualization**: Circos and coverage plots

## Testing
- Use `@pytest.mark.external_tool` for tests requiring samtools, SURVIVOR
- Generate synthetic BAM/VCF data programmatically
- All test outputs to `tmp_path`
