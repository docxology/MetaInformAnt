# Annotation

Structural variant annotation providing gene and regulatory element overlap detection, functional impact prediction, dosage sensitivity assessment, and TAD boundary disruption analysis.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports overlap and functional_impact submodules |
| `overlap.py` | Interval-based gene/regulatory overlap detection and nearest-gene finding |
| `functional_impact.py` | Functional impact prediction, dosage sensitivity, TAD disruption scoring |

## Key Functions

| Function | Description |
|----------|-------------|
| `overlap.annotate_gene_overlap()` | Find genes overlapping structural variants |
| `overlap.annotate_regulatory_overlap()` | Find regulatory elements overlapping SVs |
| `overlap.calculate_overlap_fraction()` | Compute reciprocal overlap fraction between intervals |
| `overlap.find_nearest_gene()` | Find nearest gene to a structural variant |
| `functional_impact.predict_functional_impact()` | Predict functional consequence of an SV |
| `functional_impact.assess_dosage_sensitivity()` | Assess gene dosage sensitivity for affected genes |
| `functional_impact.predict_tad_disruption()` | Predict TAD boundary disruption by SVs |
| `functional_impact.score_pathogenicity()` | Compute aggregate pathogenicity score |

## Usage

```python
from metainformant.structural_variants.annotation import overlap, functional_impact

overlaps = overlap.annotate_gene_overlap(variant, gene_annotations)
nearest = overlap.find_nearest_gene(variant, gene_annotations)
impact = functional_impact.predict_functional_impact(variant, gene_annotations)
score = functional_impact.score_pathogenicity(variant, gene_annotations, tad_boundaries)
```
