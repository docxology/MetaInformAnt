# Specification: annotation

## 🎯 Scope
Structural variant annotation subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 3 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `functional_impact.py`
- `overlap.py`

### Function contracts (behavior pinned by tests)

- `overlap.IntervalIndex.query_overlap(chrom, start, end)` returns every interval with
  `interval.start < end and interval.end > start`, sorted by `(start, end)`, in `O(k log n)` for `k`
  hits (a segment tree of maximum ends), however long the intervals are.
- `overlap.IntervalIndex.query_nearest(chrom, position, max_distance)` considers every interval of
  the chromosome. An interval that spans the position is returned at distance 0 (the one that
  starts last when several do); otherwise the closer of the interval with the smallest start at
  or after the position and the interval with the largest end before it (the former on a tie);
  `(None, -1)` when none is within `max_distance`.
- `functional_impact.predict_functional_impact(variant, gene_annotations, haploinsufficiency_db,
  tad_boundaries, triplosensitivity_db, pli_db, loeuf_db)` classifies genes with every database it
  is given: a duplication over a triplosensitive gene (TS score >= 0.5) is `HIGH` /
  `dosage_gain`, and a deletion over a gene that is haploinsufficient (HI score >= 0.5,
  pLI >= 0.9 or LOEUF < 0.35) is `HIGH` / `dosage_loss`. The databases after `tad_boundaries` are
  trailing keyword parameters: positional calls with up to four arguments keep their meaning.
- `functional_impact.predict_tad_disruption(...)` tests the `chrom`/`start`/`end` of the variant
  only, so a translocation (`end == start`) is tested at its first breakpoint; its
  `genes_in_affected_tads` is sorted.
