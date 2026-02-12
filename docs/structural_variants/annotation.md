# Structural Variants: Annotation

Gene overlap, regulatory annotation, and functional impact prediction.

## Concept Overview

Once structural variants are detected, annotation determines their biological
significance. This submodule provides two complementary layers:

- **Overlap Annotation** -- Determine which genes and regulatory elements a variant
  intersects, using interval indexing for efficient queries
- **Functional Impact** -- Predict pathogenicity by combining gene overlap with dosage
  sensitivity, TAD boundary disruption, and scoring models

Together these modules answer: "What does this variant affect, and how bad is it?"

## Overlap Annotation

### Core Types

**GenomicInterval** -- A genomic region with metadata.

| Field | Type | Description |
|-------|------|-------------|
| `chrom` | `str` | Chromosome |
| `start` | `int` | Start position |
| `end` | `int` | End position |
| `name` | `str` | Feature name (gene symbol, regulatory ID) |
| `feature_type` | `str` | `"gene"`, `"exon"`, `"enhancer"`, `"promoter"`, etc. |
| `strand` | `str` | `"+"`, `"-"`, or `"."` |
| `info` | `Dict` | Additional metadata |

**OverlapResult** -- A single variant-feature overlap.

| Field | Type | Description |
|-------|------|-------------|
| `variant_chrom` | `str` | Variant chromosome |
| `variant_start` | `int` | Variant start |
| `variant_end` | `int` | Variant end |
| `feature` | `GenomicInterval` | Overlapping feature |
| `overlap_bp` | `int` | Base pairs of overlap |
| `overlap_fraction_variant` | `float` | Fraction of variant overlapping feature |
| `overlap_fraction_feature` | `float` | Fraction of feature overlapping variant |
| `relationship` | `str` | `"contains"`, `"within"`, `"overlaps_start"`, `"overlaps_end"` |

**IntervalIndex** -- Sorted interval collection with binary search for fast queries.

### annotate_gene_overlap

Find genes overlapping structural variants.

```python
from metainformant.structural_variants.annotation.overlap import (
    annotate_gene_overlap,
    GenomicInterval,
)

gene_db = [
    GenomicInterval("chr1", 1000, 5000, "GENE_A", "gene", "+", {}),
    GenomicInterval("chr1", 8000, 12000, "GENE_B", "gene", "-", {}),
]

overlaps = annotate_gene_overlap(variants, gene_db)
for result in overlaps:
    print(f"Variant at {result.variant_chrom}:{result.variant_start}-{result.variant_end} "
          f"overlaps {result.feature.name} ({result.overlap_fraction_variant:.1%})")
```

### annotate_regulatory_overlap

Find regulatory elements (enhancers, promoters, CTCF sites) overlapping variants.

```python
from metainformant.structural_variants.annotation.overlap import annotate_regulatory_overlap

regulatory_db = [
    GenomicInterval("chr1", 900, 1100, "PROM_A", "promoter", ".", {}),
    GenomicInterval("chr1", 6000, 6500, "ENH_1", "enhancer", ".", {}),
]

reg_overlaps = annotate_regulatory_overlap(variants, regulatory_db)
```

### Supporting Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `calculate_overlap_fraction` | `(interval_a, interval_b)` | Reciprocal overlap between two intervals |
| `find_nearest_gene` | `(variant, gene_db, max_distance=100_000)` | Nearest gene when no direct overlap |

## Functional Impact Prediction

### Core Types

**FunctionalImpact** -- Comprehensive impact assessment for a variant.

| Field | Type | Description |
|-------|------|-------------|
| `variant_id` | `str` | Variant identifier |
| `impact_level` | `str` | `"HIGH"`, `"MODERATE"`, `"LOW"`, `"MODIFIER"` |
| `impact_type` | `str` | Specific mechanism (e.g., `"gene_deletion"`, `"enhancer_disruption"`) |
| `affected_genes` | `List[str]` | Genes affected by variant |
| `dosage_sensitive_genes` | `List[str]` | Haploinsufficient or triplosensitive genes |
| `tad_disrupted` | `bool` | Whether a TAD boundary is disrupted |
| `pathogenicity_score` | `float` | Combined pathogenicity score [0, 1] |
| `details` | `Dict` | Full annotation details |

**DosageSensitivity** -- Gene-level dosage sensitivity assessment.

| Field | Type | Description |
|-------|------|-------------|
| `gene_name` | `str` | Gene symbol |
| `haploinsufficiency_score` | `float` | HI score |
| `triplosensitivity_score` | `float` | TS score |
| `pli_score` | `float` | pLI (probability of LoF intolerance) |
| `loeuf` | `float` | LOEUF (LoF observed/expected upper bound) |
| `is_haploinsufficient` | `bool` | Exceeds HI threshold |
| `is_triplosensitive` | `bool` | Exceeds TS threshold |

**TADPrediction** -- TAD boundary disruption prediction.

| Field | Type | Description |
|-------|------|-------------|
| `disrupted_boundaries` | `List` | Affected TAD boundaries |
| `n_boundaries_disrupted` | `int` | Count of disrupted boundaries |
| `genes_in_affected_tads` | `List[str]` | Genes within affected TADs |
| `disruption_score` | `float` | Disruption severity [0, 1] |

### predict_functional_impact

Comprehensive impact prediction combining all evidence sources.

```python
from metainformant.structural_variants.annotation.functional_impact import (
    predict_functional_impact,
)

impact = predict_functional_impact(
    variant=sv,
    gene_annotations=gene_db,
    haploinsufficiency_db=hi_scores,
    tad_boundaries=tad_bounds
)

print(f"Impact: {impact.impact_level} - {impact.impact_type}")
print(f"Pathogenicity: {impact.pathogenicity_score:.3f}")
print(f"Affected genes: {', '.join(impact.affected_genes)}")
if impact.tad_disrupted:
    print("WARNING: TAD boundary disruption detected")
```

### assess_dosage_sensitivity

Evaluate a gene's sensitivity to copy number changes.

```python
from metainformant.structural_variants.annotation.functional_impact import (
    assess_dosage_sensitivity,
)

sensitivity = assess_dosage_sensitivity(
    gene="BRCA1",
    haploinsufficiency_db=hi_db,
    triplosensitivity_db=ts_db,
    pli_db=pli_db,
    loeuf_db=loeuf_db
)

if sensitivity.is_haploinsufficient:
    print(f"{sensitivity.gene_name}: haploinsufficient (pLI={sensitivity.pli_score:.3f})")
```

### predict_tad_disruption

Predict whether a variant disrupts topologically associating domain boundaries.

```python
from metainformant.structural_variants.annotation.functional_impact import (
    predict_tad_disruption,
)

tad = predict_tad_disruption(
    variant=sv,
    tad_boundaries=boundaries,
    boundary_window=10_000
)
print(f"Boundaries disrupted: {tad.n_boundaries_disrupted}")
print(f"Disruption score: {tad.disruption_score:.3f}")
```

### score_pathogenicity

Combined pathogenicity score from all annotation sources.

```python
from metainformant.structural_variants.annotation.functional_impact import (
    score_pathogenicity,
)

score = score_pathogenicity(variant=sv, annotations=all_annotations)
print(f"Pathogenicity: {score:.3f}")  # 0.0 (benign) to 1.0 (pathogenic)
```

## Configuration

Annotation parameters use the `SV_` environment prefix:

| Parameter | Default | Env Variable |
|-----------|---------|-------------|
| `max_distance` | `100000` | `SV_MAX_GENE_DISTANCE` |
| `boundary_window` | `10000` | `SV_TAD_WINDOW` |

## See Also

- [Detection](detection.md) -- Detecting SVs to annotate
- [Filtering](filtering.md) -- Filtering annotated variants
- [Population](population.md) -- Population frequency annotation
