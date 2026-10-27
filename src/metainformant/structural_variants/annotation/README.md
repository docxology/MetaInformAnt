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
| `overlap.IntervalIndex` | Overlap and nearest-interval queries over a max-end segment tree; finds intervals of any length |
| `overlap.annotate_gene_overlap()` | Find genes overlapping structural variants |
| `overlap.annotate_regulatory_overlap()` | Find regulatory elements overlapping SVs |
| `overlap.calculate_overlap_fraction()` | Compute reciprocal overlap fraction between intervals |
| `overlap.find_nearest_gene()` | Find nearest gene to a structural variant |
| `functional_impact.predict_functional_impact()` | Predict functional consequence of an SV from the haploinsufficiency, triplosensitivity, pLI and LOEUF databases you supply |
| `functional_impact.assess_dosage_sensitivity()` | Assess gene dosage sensitivity for affected genes |
| `functional_impact.predict_tad_disruption()` | Predict TAD boundary disruption by SVs (a translocation is checked at its first breakpoint only) |
| `functional_impact.score_pathogenicity()` | Compute aggregate pathogenicity score |

## Usage

```python
from metainformant.structural_variants.annotation import overlap, functional_impact

variant = {"chrom": "chr1", "start": 1_000, "end": 90_000, "sv_type": "DEL"}
gene_annotations = [{"chrom": "chr1", "start": 5_000, "end": 20_000, "name": "GENE_A"}]

annotated = overlap.annotate_gene_overlap([variant], gene_annotations)  # list in, list out
nearest = overlap.find_nearest_gene(variant, gene_annotations)           # one variant dict
impact = functional_impact.predict_functional_impact(
    variant, gene_annotations, haploinsufficiency_db={"GENE_A": 0.9}
)
score = functional_impact.score_pathogenicity(
    variant,
    {"overlapping_genes": ["GENE_A"], "dosage_sensitive_genes": ["GENE_A"], "n_coding_genes": 1},
)
```
