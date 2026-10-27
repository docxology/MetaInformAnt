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

The annotation functions take variant dictionaries (`chrom`, `start`, `end`, `sv_type`, ...);
`dataclasses.asdict` turns a `StructuralVariant` into one. They read `chrom`, `start` and `end`
only. A translocation call has `end == start` (its first breakpoint) and carries the second
breakpoint in `chrom2` and `pos2`, which the annotation functions do not look at, so a
translocation is annotated at its first breakpoint.

## Overlap Annotation

### Core Types

**GenomicInterval** -- A genomic region with metadata.

| Field | Type | Description |
|-------|------|-------------|
| `chrom` | `str` | Chromosome |
| `start` | `int` | Start position (0-based, inclusive) |
| `end` | `int` | End position (0-based, exclusive) |
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
| `relationship` | `str` | `"contains"` (variant spans the feature), `"contained_in"` (variant lies inside the feature) or `"overlap"` |

**IntervalIndex** -- Interval collection for overlap and nearest-interval queries. It groups
intervals by chromosome, sorts them by `(start, end)` and keeps a segment tree of the maximum
interval end, so `query_overlap(chrom, start, end)` returns every interval with
`interval.start < end and interval.end > start`, sorted by `(start, end)`, in `O(k log n)` time
for `k` hits, however long the intervals are. `query_nearest(chrom, position, max_distance)`
considers every interval of the chromosome: an interval that spans the position is returned at
distance 0 (the one that starts last when several do); otherwise the closer of the interval with
the smallest start at or after the position and the interval with the largest end before it wins,
the former on a tie; `(None, -1)` means that none lies within `max_distance`.

### annotate_gene_overlap

Find genes overlapping structural variants. The input is a list of variant dictionaries and the
result is a list of the same dictionaries with the annotation keys added.

```python
from metainformant.structural_variants.annotation.overlap import (
    annotate_gene_overlap,
    GenomicInterval,
)

gene_db = [
    GenomicInterval("chr1", 1000, 5000, "GENE_A", "gene", "+", {}),
    GenomicInterval("chr1", 8000, 12000, "GENE_B", "gene", "-", {}),
]

annotated = annotate_gene_overlap(variants, gene_db)
for variant in annotated:
    for result in variant["gene_overlaps"]:
        print(f"Variant at {result.variant_chrom}:{result.variant_start}-{result.variant_end} "
              f"overlaps {result.feature.name} ({result.overlap_fraction_variant:.1%})")
```

Each annotated dictionary has `overlapping_genes` (names), `gene_overlaps` (`OverlapResult` objects) and `n_genes_affected`.

### annotate_regulatory_overlap

Find regulatory elements (enhancers, promoters, CTCF sites) overlapping variants.

```python
from metainformant.structural_variants.annotation.overlap import (
    annotate_regulatory_overlap,
    GenomicInterval,
)

regulatory_db = [
    GenomicInterval("chr1", 900, 1100, "PROM_A", "promoter", ".", {}),
    GenomicInterval("chr1", 6000, 6500, "ENH_1", "enhancer", ".", {}),
]

reg_overlaps = annotate_regulatory_overlap(variants, regulatory_db)
```

Each annotated dictionary has `overlapping_regulatory` (names), `regulatory_overlaps`
(`OverlapResult` objects), `regulatory_types` (a set of the feature types hit) and
`n_regulatory_affected`.

### Supporting Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `calculate_overlap_fraction` | `(interval_a, interval_b)` | Reciprocal overlap between two intervals |
| `find_nearest_gene` | `(variant, gene_db, max_distance=100_000)` | The gene with the largest overlap, else the nearest gene within `max_distance`; returns a dictionary with `nearest_gene`, `distance`, `direction` and `gene_interval` |

## Functional Impact Prediction

### Core Types

**FunctionalImpact** -- Comprehensive impact assessment for a variant.

| Field | Type | Description |
|-------|------|-------------|
| `variant_id` | `str` | Variant identifier (the variant's `id`, else `chrom:start-end`) |
| `impact_level` | `str` | `"HIGH"`, `"MODERATE"`, `"LOW"`, `"MODIFIER"` |
| `impact_type` | `str` | Specific mechanism: `dosage_loss`, `dosage_gain`, `gene_disruption`, `gene_fusion_candidate`, `gene_duplication`, `gene_inversion`, `gene_insertion`, `tad_disruption`, `large_intergenic`, `intergenic` or `benign` |
| `affected_genes` | `List[str]` | Genes affected by variant |
| `dosage_sensitive_genes` | `List[str]` | Haploinsufficient genes hit by a deletion, triplosensitive genes hit by a duplication |
| `tad_disrupted` | `bool` | Whether a TAD boundary is disrupted |
| `pathogenicity_score` | `float` | Combined pathogenicity score [0, 1] |
| `details` | `Dict` | `sv_type`, `size`, `n_overlapping_genes` and, with TAD boundaries, `tad_prediction` |

**DosageSensitivity** -- Gene-level dosage sensitivity assessment.

| Field | Type | Description |
|-------|------|-------------|
| `gene_name` | `str` | Gene symbol |
| `haploinsufficiency_score` | `float` | HI score |
| `triplosensitivity_score` | `float` | TS score |
| `pli_score` | `float` | pLI (probability of LoF intolerance) |
| `loeuf` | `float` | LOEUF (LoF observed/expected upper bound) |
| `is_haploinsufficient` | `bool` | HI score >= 0.5, pLI >= 0.9 or LOEUF < 0.35 |
| `is_triplosensitive` | `bool` | TS score >= 0.5 |

**TADPrediction** -- TAD boundary disruption prediction.

| Field | Type | Description |
|-------|------|-------------|
| `variant_chrom`, `variant_start`, `variant_end` | `str`, `int`, `int` | The variant that was tested |
| `disrupted_boundaries` | `List[int]` | Midpoints of the disrupted TAD boundaries |
| `n_boundaries_disrupted` | `int` | Count of disrupted boundaries |
| `genes_in_affected_tads` | `List[str]` | Genes within affected TADs (sorted) |
| `disruption_score` | `float` | Disruption severity [0, 1] |

### predict_functional_impact

Comprehensive impact prediction combining all evidence sources. The four databases are
optional dictionaries mapping gene names to scores; every one you pass takes part in the
dosage classification.

```python
from metainformant.structural_variants.annotation.functional_impact import (
    predict_functional_impact,
)

impact = predict_functional_impact(
    variant=sv,
    gene_annotations=gene_db,
    haploinsufficiency_db=hi_scores,
    tad_boundaries=tad_bounds,
    triplosensitivity_db=ts_scores,
    pli_db=pli_scores,
    loeuf_db=loeuf_values,
)

print(f"Impact: {impact.impact_level} - {impact.impact_type}")
print(f"Pathogenicity: {impact.pathogenicity_score:.3f}")
print(f"Affected genes: {', '.join(impact.affected_genes)}")
if impact.tad_disrupted:
    print("WARNING: TAD boundary disruption detected")
```

A deletion (DEL or HOMODEL) over a haploinsufficient gene is `HIGH` / `dosage_loss`, and a duplication over a triplosensitive gene (`triplosensitivity_db`) is `HIGH` / `dosage_gain`. Without a database no gene is marked dosage-sensitive, and a duplication over genes is `MODERATE` / `gene_duplication`. `gene_annotations` is a list of dictionaries with `chrom`, `start`, `end`, `name` and optionally `is_coding`.

### assess_dosage_sensitivity

Evaluate a gene's sensitivity to copy number changes. A gene is haploinsufficient when its HI
score is at least 0.5, its pLI at least 0.9, or its LOEUF below 0.35, and triplosensitive when
its TS score is at least 0.5; a gene missing from a dictionary scores 0 (1.0 for LOEUF).

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

Predict whether a variant disrupts topologically associating domain boundaries. It tests the
variant's `chrom`, `start` and `end`, so a translocation (`end == start`) is tested at its first
breakpoint only.

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

Combined pathogenicity score from the annotation results of one variant. `annotations` is a
dictionary with `overlapping_genes`, `dosage_sensitive_genes`, `tad_disrupted`, `impact_level`,
`n_coding_genes` and optionally `population_frequency`.

```python
from metainformant.structural_variants.annotation.functional_impact import (
    score_pathogenicity,
)

score = score_pathogenicity(variant=sv, annotations=all_annotations)
print(f"Pathogenicity: {score:.3f}")  # 0.0 (benign) to 1.0 (pathogenic)
```

## Configuration

Annotation parameters are function arguments; none has an environment variable:

| Parameter | Default | Function |
|-----------|---------|----------|
| `max_distance` | `100000` | `find_nearest_gene` |
| `boundary_window` | `10000` | `predict_tad_disruption` |

See [CONFIGURATION.md](CONFIGURATION.md) for every argument.

## See Also

- [Detection](detection.md) -- Detecting SVs to annotate
- [Filtering](filtering.md) -- Filtering annotated variants
- [Population](population.md) -- Population frequency annotation
