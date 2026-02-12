# Structural Variants: Filtering & Merging

Quality-based filtering and multi-caller variant merging.

## Concept Overview

Raw SV call sets contain false positives and duplicate calls. This submodule
provides two complementary cleanup operations:

- **Quality Filtering** -- Remove low-confidence calls based on quality scores,
  read support, size, population frequency, and blacklist regions
- **Multi-Caller Merging** -- Combine calls from multiple SV callers into a
  consensus set using reciprocal overlap or distance-based matching

Each filtering step returns a `FilterStats` object tracking how many variants
passed, enabling pipeline quality monitoring.

## Quality Filtering

### FilterStats

Every filter function returns `FilterStats` alongside filtered variants.

| Field | Type | Description |
|-------|------|-------------|
| `input_count` | `int` | Variants before filtering |
| `output_count` | `int` | Variants after filtering |
| `filtered_count` | `int` | Variants removed |
| `filter_name` | `str` | Name of the filter applied |
| `parameters` | `Dict` | Filter parameters used |
| `pass_rate` | `float` | Property: `output_count / input_count` |

### filter_by_quality

Filter by quality score, read support, and mapping quality.

```python
from metainformant.structural_variants.filtering.quality_filter import filter_by_quality

filtered, stats = filter_by_quality(
    variants=raw_calls,
    min_qual=20.0,
    min_support=3,
    min_mapq=0.0,
    min_gq=0.0
)
print(f"Passed: {stats.output_count}/{stats.input_count} ({stats.pass_rate:.1%})")
```

### filter_by_size

Remove variants outside a size range.

```python
from metainformant.structural_variants.filtering.quality_filter import filter_by_size

filtered, stats = filter_by_size(
    variants=calls,
    min_size=50,        # Minimum 50 bp
    max_size=5_000_000  # Maximum 5 Mb (None for no upper limit)
)
```

### filter_by_frequency

Remove common variants based on population allele frequency.

```python
from metainformant.structural_variants.filtering.quality_filter import filter_by_frequency

filtered, stats = filter_by_frequency(
    variants=calls,
    population_db=gnomad_sv,        # Population frequency database
    max_af=0.01,                    # Keep variants with AF < 1%
    match_window=200,               # Position matching window
    min_reciprocal_overlap=0.5      # Overlap threshold for matching
)
```

### apply_blacklist

Remove variants overlapping known problematic regions.

```python
from metainformant.structural_variants.filtering.quality_filter import apply_blacklist

filtered, stats = apply_blacklist(
    variants=calls,
    blacklist_regions=encode_blacklist,  # List of GenomicInterval
    min_overlap_fraction=0.0            # Any overlap triggers removal
)
```

## Multi-Caller Merging

### Core Types

**MergedVariant** -- A consensus variant from multiple callers.

| Field | Type | Description |
|-------|------|-------------|
| `chrom` | `str` | Chromosome |
| `start` | `int` | Consensus start position |
| `end` | `int` | Consensus end position |
| `sv_type` | `str` | SV type (DEL, DUP, INV, etc.) |
| `size` | `int` | Variant size |
| `n_callers` | `int` | Number of callers supporting this variant |
| `caller_names` | `List[str]` | Which callers called it |
| `support_variants` | `List` | Original variants from each caller |
| `consensus_quality` | `float` | Combined quality score |
| `genotype` | `str` | Consensus genotype |

**MergeStats** -- Summary of the merging operation.

| Field | Type | Description |
|-------|------|-------------|
| `n_input_callsets` | `int` | Number of callers merged |
| `n_input_variants` | `int` | Total input variants |
| `n_output_variants` | `int` | Merged output variants |
| `n_multi_caller` | `int` | Variants supported by 2+ callers |
| `n_single_caller` | `int` | Variants from only one caller |
| `caller_counts` | `Dict[str, int]` | Variants per caller |

### merge_callsets

Merge variants from multiple callers using reciprocal overlap.

```python
from metainformant.structural_variants.filtering.merge import merge_callsets

callsets = {
    "delly": delly_calls,
    "manta": manta_calls,
    "lumpy": lumpy_calls,
}

merged, stats = merge_callsets(
    callsets=callsets,
    min_overlap=0.5,     # 50% reciprocal overlap required
    type_match=True      # SV types must match
)

print(f"Merged: {stats.n_input_variants} -> {stats.n_output_variants}")
print(f"Multi-caller: {stats.n_multi_caller}")
```

### survivor_merge

SURVIVOR-style merging based on breakpoint distance.

```python
from metainformant.structural_variants.filtering.merge import survivor_merge

merged = survivor_merge(
    vcf_files=["delly.vcf", "manta.vcf", "lumpy.vcf"],
    max_distance=1000,      # Maximum breakpoint distance
    min_callers=2,          # Minimum callers for inclusion
    type_match=True,
    strand_match=False
)
```

### Supporting Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `calculate_reciprocal_overlap` | `(sv1, sv2) -> float` | Min of two overlap fractions |
| `deduplicate_variants` | `(variants, max_distance=100, min_reciprocal_overlap=0.8, type_match=True)` | Remove duplicates within one callset |

## End-to-End Pipeline

```python
from metainformant.structural_variants.filtering.quality_filter import (
    filter_by_quality, filter_by_size, filter_by_frequency, apply_blacklist,
)
from metainformant.structural_variants.filtering.merge import merge_callsets

# Step 1: Merge multi-caller results
merged, merge_stats = merge_callsets(
    {"delly": delly, "manta": manta}, min_overlap=0.5
)

# Step 2: Quality filter cascade
variants = merged
variants, q_stats = filter_by_quality(variants, min_qual=20, min_support=3)
variants, s_stats = filter_by_size(variants, min_size=50, max_size=5_000_000)
variants, f_stats = filter_by_frequency(variants, population_db=gnomad, max_af=0.01)
variants, b_stats = filter_by_blacklist(variants, blacklist_regions=blacklist)

# Report
for stats in [q_stats, s_stats, f_stats, b_stats]:
    print(f"{stats.filter_name}: {stats.input_count} -> {stats.output_count} "
          f"({stats.pass_rate:.1%})")
```

## Configuration

Filtering parameters use the `SV_` environment prefix:

| Parameter | Default | Env Variable |
|-----------|---------|-------------|
| `min_qual` | `20.0` | `SV_MIN_QUAL` |
| `min_support` | `3` | `SV_MIN_SUPPORT` |
| `min_size` | `50` | `SV_MIN_SIZE` |
| `max_af` | `0.01` | `SV_MAX_AF` |
| `min_overlap` | `0.5` | `SV_MIN_OVERLAP` |

## See Also

- [Detection](detection.md) -- Generating raw SV calls
- [Annotation](annotation.md) -- Annotating filtered variants
- [Population](population.md) -- Population-scale analysis
