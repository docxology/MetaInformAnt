# Filtering

Structural variant filtering and multi-caller consensus merging, providing quality-based filtering, size filtering, population frequency filtering, blacklist exclusion, SURVIVOR-style merging, and deduplication.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports quality_filter and merge submodules |
| `quality_filter.py` | Quality, size, frequency, and blacklist filtering of SV calls |
| `merge.py` | Multi-caller merging, reciprocal overlap, SURVIVOR-style consensus |

## Key Functions

| Function | Description |
|----------|-------------|
| `quality_filter.filter_by_quality()` | Filter SVs by quality score and read support |
| `quality_filter.filter_by_size()` | Filter SVs by minimum/maximum size; TRA/BND and unmeasured (size 0) INS calls are kept |
| `quality_filter.filter_by_frequency()` | Filter SVs by population allele frequency; raises `ValueError` when no database entry is readable |
| `quality_filter.apply_blacklist()` | Remove SVs overlapping blacklist regions |
| `merge.merge_callsets()` | Merge SV calls from multiple callers by reciprocal overlap |
| `merge.calculate_reciprocal_overlap()` | Compute reciprocal overlap between two SVs |
| `merge.survivor_merge()` | SURVIVOR-style distance-based merging; majority-vote genotype and agreement-weighted quality |
| `merge.deduplicate_variants()` | Remove duplicate SV calls |

## Usage

```python
from metainformant.structural_variants.filtering import quality_filter, merge

filtered, _ = quality_filter.filter_by_quality(variants, min_qual=20, min_support=3)
filtered, _ = quality_filter.filter_by_size(filtered, min_size=50, max_size=10_000_000)
merged, merge_stats = merge.merge_callsets({"caller1": caller1_svs, "caller2": caller2_svs}, min_overlap=0.5)
deduped, n_removed = merge.deduplicate_variants(filtered)
```

Native VCF parsing checks for optional `CHR2` and `POS2` fields before reading
them. A BND record can describe its mate entirely in ALT without declaring
these INFO fields. Errors in required fallback metadata propagate; they do not trigger
a text-parser fallback.
