# Structural Variants

## Overview

Structural variant analysis module for METAINFORMANT. Covers CNV detection via circular binary segmentation, SV calling from split/discordant reads, annotation, filtering, population genotyping, and visualization.

## Contents

- **detection/** - CNV detection (CBS), SV calling (split/discordant reads), breakpoint refinement
- **annotation/** - Gene/regulatory overlap annotation, functional impact, TAD disruption
- **filtering/** - Quality filtering, size filtering, blacklist regions, multi-caller merging
- **population/** - Population genotyping, allele frequency, association testing, PCA, LD
- **visualization/** - Circos plots, coverage tracks, size distributions, CNV profiles

## Architecture

```mermaid
graph TD
    subgraph "Structural Variants Module"
        DT[detection/] --> |cnv.py| CN[CNV Detection (CBS)]
        DT --> |sv_calling.py| SC[SV Calling (Split/Discordant)]
        DT --> |breakpoints.py| BP[Breakpoint Refinement]

        AN[annotation/] --> |overlap.py| OV[Gene & Regulatory Overlap]
        AN --> |functional_impact.py| FI[Dosage Sensitivity / TAD]

        FL[filtering/] --> |quality_filter.py| QF[Quality & Size Filtering]
        FL --> |merge.py| MG[Multi-Caller Consensus]

        PO[population/] --> |sv_population.py| SP[Population Genotyping]

        VZ[visualization/] --> |plots.py| PL[Circos, Coverage, CNV Plots]
    end

    DT --> AN
    AN --> FL
    FL --> PO
```

## Usage

```python
from metainformant.structural_variants.detection import cnv, sv_calling, breakpoints
from metainformant.structural_variants.annotation import overlap, functional_impact
from metainformant.structural_variants.filtering import quality_filter, merge
from metainformant.structural_variants.population import sv_population
from metainformant.structural_variants.visualization import plots
```
