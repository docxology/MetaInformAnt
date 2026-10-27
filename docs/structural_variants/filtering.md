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

The filters and merge functions take variant dictionaries (`chrom`, `start`, `end`, `sv_type`, ...);
`dataclasses.asdict` turns a `StructuralVariant` into one. They read `chrom`, `start` and `end`
only. A translocation call has `end == start` (its first breakpoint) and carries the second
breakpoint in `chrom2` and `pos2`, which they do not look at (see [Zero-length calls](#zero-length-calls)).

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

Remove variants outside a size range. Translocations (`TRA`, `BND`) have no meaningful size
and are always kept, and so are insertions (`INS`) without a positive `size`: an insertion
called from split reads has `size` 0, meaning "not measured". An insertion with a positive
`size` is filtered on it like any other variant.

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
    max_af=0.01,                    # Remove variants with AF > 1%
    match_window=200,               # Position matching window
    min_reciprocal_overlap=0.5      # Overlap threshold for matching
)
```

`population_db` is a dictionary `{"chrom:start-end[-svtype]": af}` or a list of dictionaries
with `chrom`, `start`, `end`, `af` (or `allele_frequency`) and optionally `sv_type`; every
frequency must lie in [0, 1]. A non-empty database none of whose entries can be read raises
`ValueError` naming the supported layouts, because the filter would otherwise pass every
variant without comparing it to anything. Unreadable entries of an otherwise usable database
are skipped and counted in a logged warning. `None` skips the filter, unless the variants carry
their own `af` values.

### apply_blacklist

Remove variants overlapping known problematic regions.

```python
from metainformant.structural_variants.filtering.quality_filter import apply_blacklist

filtered, stats = apply_blacklist(
    variants=calls,
    blacklist_regions=encode_blacklist,  # list of {"chrom", "start", "end"} dicts or (chrom, start, end) tuples
    min_overlap_fraction=0.0            # Any overlap triggers removal
)
```

A zero-length variant (`start == end`, such as an insertion point) has no bases to overlap. It
counts as fully overlapping (fraction 1.0) when it lies strictly inside a region, and not at
all on a region's edge, the same as a variant that only abuts a region.

## Multi-Caller Merging

### Core Types

**MergedVariant** -- A consensus variant from multiple callers.

| Field | Type | Description |
|-------|------|-------------|
| `chrom` | `str` | Chromosome |
| `start` | `int` | Consensus start position (median) |
| `end` | `int` | Consensus end position (median) |
| `sv_type` | `str` | SV type (DEL, DUP, INV, etc.; majority vote) |
| `size` | `int` | Variant size (`end - start`) |
| `n_callers` | `int` | Number of callers supporting this variant |
| `caller_names` | `List[str]` | Which callers called it (sorted) |
| `support_variants` | `List` | Original variants from each caller |
| `consensus_quality` | `float` | Mean call quality plus an agreement bonus (the bonus differs per merge function) |
| `genotype` | `str` | Consensus genotype (`./.` when no merged call has one) |

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

Merge variants from multiple callers using reciprocal overlap. The consensus quality is the
mean `quality` (or `qual`) of the merged calls plus `100 * log10(n_callers + 1)`, and the genotype
is the majority vote among the non-missing genotypes (a tie goes to the genotype seen first).

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

SURVIVOR-style merging based on breakpoint distance. It returns the same `(merged, stats)` tuple.

```python
from metainformant.structural_variants.filtering.merge import survivor_merge

merged, stats = survivor_merge(
    vcf_files=["delly.vcf", "manta.vcf", "lumpy.vcf"],
    max_distance=1000,      # Maximum breakpoint distance
    min_callers=2,          # Minimum callers for inclusion
    type_match=True,
    strand_match=False
)
```

Each merged variant combines the calls of its group as follows:

- `genotype`: majority vote among the non-missing genotypes (allele order ignored, so `1/0` is `0/1`). A tie goes to the genotype whose calls have the larger summed quality, then to the one that sorts first. `./.` when no call has a genotype. Genotypes are read from each dictionary's `genotype` key and from the first sample of a VCF file (a phased `0|1` becomes `0/1`).
- `consensus_quality`: the mean `quality` (or `qual`; 0 when absent) of the merged calls plus `10 * n_callers * log10(n_callers + 1)`. Two callers with calls of quality 80 (genotype `0/1`) and 40 (genotype `1/1`) give genotype `0/1` and `consensus_quality` 69.54.

### Supporting Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `calculate_reciprocal_overlap` | `(sv1, sv2) -> float` | Min of two overlap fractions |
| `deduplicate_variants` | `(variants, max_distance=100, min_reciprocal_overlap=0.8, type_match=True)` | Remove duplicates within one callset; returns `(variants, n_removed)` |

### Zero-length calls

A call with `start == end`, such as a translocation, has a reciprocal overlap of 0 with every
other call (itself included). `merge_callsets` and `deduplicate_variants` therefore never join
or collapse such calls with a threshold above 0, and `filter_by_frequency` matches no database
entry to one unless `min_reciprocal_overlap` is 0. `survivor_merge` compares breakpoint
distances instead, so it joins translocations whose first breakpoints are close whatever their
partner chromosome.

## End-to-End Pipeline

```python
from metainformant.structural_variants.filtering.quality_filter import (
    filter_by_quality, filter_by_size, filter_by_frequency, apply_blacklist,
)
from metainformant.structural_variants.filtering.merge import merge_callsets

# Step 1: filter each caller's calls (dictionaries in, (variants, FilterStats) out)
delly, q_stats = filter_by_quality(delly, min_qual=20, min_support=3)
manta, _ = filter_by_quality(manta, min_qual=20, min_support=3)

# Step 2: merge the callers by reciprocal overlap
merged, merge_stats = merge_callsets({"delly": delly, "manta": manta}, min_overlap=0.5)

# Step 3: the remaining filters take dictionaries, so describe each consensus call as one
variants = [
    {"chrom": m.chrom, "start": m.start, "end": m.end, "sv_type": m.sv_type, "size": m.size,
     "quality": m.consensus_quality, "genotype": m.genotype, "n_callers": m.n_callers}
    for m in merged
]
variants, s_stats = filter_by_size(variants, min_size=50, max_size=5_000_000)
variants, f_stats = filter_by_frequency(variants, population_db=gnomad, max_af=0.01)
variants, b_stats = apply_blacklist(variants, blacklist_regions=blacklist)

# Report
for stats in [q_stats, s_stats, f_stats, b_stats]:
    print(f"{stats.filter_name}: {stats.input_count} -> {stats.output_count} "
          f"({stats.pass_rate:.1%})")
```

## Configuration

Filtering parameters are function arguments. Two have an `SV_` environment override that
replaces the matching argument when `filter_by_quality` runs:

| Parameter | Default | Env Variable |
|-----------|---------|-------------|
| `min_qual` | `20.0` | `SV_MIN_QUAL` |
| `min_support` | `3` | `SV_MIN_SUPPORT` |

`min_size`, `max_size`, `max_af`, `min_overlap` and the rest have no environment variable. See
[CONFIGURATION.md](CONFIGURATION.md) for every argument.

## See Also

- [Detection](detection.md) -- Generating raw SV calls
- [Annotation](annotation.md) -- Annotating filtered variants
- [Population](population.md) -- Population-scale analysis
