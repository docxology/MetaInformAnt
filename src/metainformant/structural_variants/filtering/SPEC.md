# Specification: filtering

## 🎯 Scope
Structural variant filtering subpackage.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 3 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `merge.py`
- `quality_filter.py`

### Function contracts (behavior pinned by tests)

- `quality_filter.filter_by_frequency(variants, population_db=...)` reads `population_db` as a dict
  `{"chrom:start-end[-svtype]": af}` or a list of dicts with `chrom`, `start`, `end`, `af` (or
  `allele_frequency`) and optionally `sv_type`, with `af` in [0, 1]. A non-empty database none of
  whose entries can be read raises `ValueError` (the filter would otherwise silently pass every
  variant); unreadable entries of an otherwise usable database are skipped and counted in a
  logged warning.
- `quality_filter.apply_blacklist(variants, regions)` treats a zero-length variant (`start == end`)
  strictly inside a region as fully overlapping it (fraction 1.0), and one on a region's edge or
  outside it as not overlapping, like an interval that only abuts the region.
- `quality_filter.filter_by_size(variants, ...)` always keeps TRA/BND variants and INS variants
  without a positive `size` (an insertion called from split reads has `size` 0, "not measured").
- `merge.survivor_merge(...)` combines the calls of a group as follows. Genotype: majority vote
  among the non-missing genotypes, allele order ignored; a tie goes to the genotype with the larger
  summed call quality, then to the one that sorts first; `./.` when none is called. Genotypes come
  from each dict's `genotype` and from the first sample of a VCF file (phased `0|1` becomes
  `0/1`). Quality: the mean call `quality` (or `qual`; 0 when absent) plus
  `10 * n_callers * log10(n_callers + 1)`.
- `MergedVariant.caller_names` (from `merge_callsets` and `survivor_merge`) is sorted.

### Known limits

- The filters and merge functions read `chrom`, `start` and `end` only and ignore `chrom2` and
  `pos2`, which carry the second breakpoint of a translocation (`SVType.TRA`, `end == start`).
  `survivor_merge` joins translocations whose first breakpoints are within `max_distance`, whatever
  their partner. `merge_callsets` and `deduplicate_variants` give a zero-length call a reciprocal
  overlap of 0, so with a threshold above 0 they never join or collapse translocations, and
  `filter_by_frequency` matches no database entry to a zero-length call unless
  `min_reciprocal_overlap` is 0.
