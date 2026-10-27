# DNA Methylation Analysis

Analysis of DNA methylation patterns, implemented in `metainformant.epigenome.assays.methylation`: CpG site modeling, beta value computation, differentially methylated region (DMR) detection, methylation-based CpG island screening, and methylation entropy.

## Key Concepts

**Beta values** represent the methylation level at each CpG site as the ratio of methylated reads to total reads: beta = M / (M + U), where M and U are the methylated and unmethylated read counts. Values range from 0.0 (unmethylated) to 1.0 (fully methylated); sites with no reads get 0.0. `MethylationSite.methylation_level` is the same ratio for site objects.

**Differentially Methylated Regions (DMRs)** are contiguous runs of CpG sites whose methylation differs between conditions. Sites present in both conditions are compared one by one (a minimum methylation difference plus a two-proportion z-test on the read counts), and consecutive significant sites are chained into regions of at least `min_sites` sites. There is no sliding window.

**CpG islands** are regions of high CG dinucleotide density, typically defined by observed/expected CpG ratio >= 0.6, GC content >= 50%, and length >= 200 bp. `identify_cpg_islands` does not read DNA sequence: it estimates GC content and CpG ratio from the methylation levels of the sites in each window, so treat its output as a coarse methylation-density screen rather than a sequence-based island call.

**Methylation entropy** quantifies how heterogeneous methylation is across the CpG sites of a genomic window: the Shannon entropy (in bits) of the distribution of site methylation levels over five equal-width bins. Windows whose sites share one level score near 0; windows that mix low, intermediate and high levels score higher.

## Data Model

### `MethylationSite`

Class representing a single CpG site measurement (a plain class, not a dataclass). The constructor requires a non-empty `chromosome`, `position >= 0`, `methylated_reads >= 0` and `total_reads >= methylated_reads`.

| Field | Type | Description |
|-------|------|-------------|
| `chromosome` | `str` | Chromosome name |
| `position` | `int` | Genomic coordinate |
| `methylated_reads` | `int` | Number of methylated reads |
| `total_reads` | `int` | Total reads covering the site (not the unmethylated count) |
| `strand` | `str` | Strand (default `"+"`) |

Read-only properties: `methylation_level` (`methylated_reads / total_reads`, 0.0 when there are no reads) and `coverage` (the same as `total_reads`). `to_dict()` returns all of these. There are no unmethylated-count, beta-value or sequence-context fields.

## Function Reference

### `load_cpg_table(path) -> pandas.DataFrame`

Read a tab-delimited table with the columns `chrom`, `pos`, `methylated` and `unmethylated` (read counts) into a DataFrame; raises `ValueError` if one of them is missing.

The functions below that take `methylation_data` work on a `Dict[str, List[MethylationSite]]` keyed by chromosome. Load one with `load_methylation_bedgraph(path, min_coverage=1)` (BEDgraph with a 0--1 level in column 4; it carries no read counts, so every site is given `max(min_coverage, 100)` estimated reads with the methylated count set from the level, and `min_coverage` never removes a site) or `load_methylation_cov(path, min_coverage=1)`, which reads a header-less Bismark `.cov` file (read as gzip when the name ends in `.gz`) with six tab-separated columns: chromosome, start, end, methylation percentage, methylated count and unmethylated count. The site position is `start`, its coverage is methylated plus unmethylated reads (`min_coverage` applies to that total) and its level is computed from the counts; the percentage column only has to parse as a number. A coverage file records no strand, so every site keeps the default `+`. Blank lines, `#` lines and malformed lines (for example an older layout with a strand symbol in column 3) are skipped, the malformed ones with a logged warning.

### `compute_beta_values(df) -> pandas.DataFrame`

Return a copy of `df` (which needs `methylated` and `unmethylated` columns) with a `beta` column: `beta = methylated / (methylated + unmethylated)`, set to 0.0 where both counts are zero. No offset or other regularization is applied.

### `summarize_beta_by_chromosome(df) -> pandas.DataFrame`

Per-chromosome summary of the `beta` column (`df` needs `chrom` and `beta`): a DataFrame indexed by chromosome with the columns `mean`, `std`, `min`, `max` and `count`.

### `find_differentially_methylated_regions(site_data1, site_data2, delta_threshold=0.2, p_value_threshold=0.05, min_sites=3, max_gap=None) -> List[Dict]`

Find DMRs between two conditions, each a `Dict[str, List[MethylationSite]]`. Only positions present in both conditions on a shared chromosome are compared. A position is significant when its absolute methylation difference is at least `delta_threshold` and its two-sided two-proportion z-test p-value (normal approximation from the read counts) is at most `p_value_threshold`. The standard error is built from each condition's own proportion, so when both levels are exactly 0 or 1 (for example 0% against 100% methylated) it is zero, the code returns p = 1 and the position is never significant. Runs of consecutive significant positions with at least `min_sites` positions become regions; `max_gap` (bp) optionally ends a run when neighboring significant positions are further apart (default `None`: no gap limit).

Returns dicts with `chromosome`, `start`, `end`, `sites` (the list of positions), `mean_delta`, `mean_p_value` (the smallest per-site p-value in the region) and `direction` (`"hyper"` or `"hypo"`: condition 1 relative to condition 2 at the region's first site).

### `identify_cpg_islands(methylation_data, window_size=200, step_size=50, min_gc=0.5, min_cpg_ratio=0.6) -> List[Dict]`

Slide a window along each chromosome that has at least 10 sites and keep windows with at least 5 sites whose estimated GC content and CpG ratio reach `min_gc` and `min_cpg_ratio`. Both values are derived from the window's mean methylation level (`gc_content = min(1, 0.3 + 0.4 * mean)` and `cpg_ratio = 1.2 * mean`); no sequence is read, so with the default thresholds a window qualifies when its mean methylation is at least 0.5. Overlapping windows are resolved by keeping the one with the highest CpG ratio.

Returns dicts with `chromosome`, `start`, `end`, `sites`, `avg_methylation`, `gc_content` and `cpg_ratio`.

### `calculate_methylation_entropy(methylation_data, window_size=1000) -> Dict[str, List[Tuple[int, float]]]`

Shannon entropy (bits) of the site methylation levels in sliding windows of `window_size` bp (step `window_size // 2`), per chromosome; chromosomes with fewer than 10 sites and windows with fewer than 5 sites are skipped. Returns `{chromosome: [(window_center, entropy), ...]}`. The value is 0 when every site falls in the same 0.2-wide level bin and reaches log2(5), about 2.32, when sites are spread evenly over the five bins. A site at exactly 1.0 falls in none of the bins (the last bin is [0.8, 1.0)) but still counts toward the window's site total, so such sites distort the value (a window of only 1.0 sites scores 0).

## Usage Examples

```python-snippet
from metainformant.epigenome.assays.methylation import (
    MethylationSite,
    calculate_methylation_entropy,
    compute_beta_values,
    find_differentially_methylated_regions,
    identify_cpg_islands,
    load_cpg_table,
)

# Beta values from a count table (columns: chrom, pos, methylated, unmethylated)
table = load_cpg_table("cpg_counts.tsv")
beta = compute_beta_values(table)
for row in beta.head(3).itertuples():
    print(f"{row.chrom}:{row.pos} beta={row.beta:.3f}")


def make_sites(levels, step=10):
    """One chromosome of sites with 100 reads each: (methylated reads, total reads)."""
    return {
        "chr1": [
            MethylationSite("chr1", 1000 + step * i, round(level * 100), 100)
            for i, level in enumerate(levels)
        ]
    }


control_sites = make_sites([0.1] * 30 + [0.8] * 30)
treated_sites = make_sites([0.7] * 30 + [0.8] * 30)

# Find DMRs between conditions
dmrs = find_differentially_methylated_regions(control_sites, treated_sites, delta_threshold=0.2)
for dmr in dmrs:
    print(
        f"{dmr['chromosome']}:{dmr['start']}-{dmr['end']} "
        f"mean_delta={dmr['mean_delta']:.3f} ({dmr['direction']})"
    )

# Methylation-based CpG island screen (estimated from methylation levels, not sequence)
islands = identify_cpg_islands(treated_sites)

# Methylation entropy in 200 bp windows: {chromosome: [(window_center, entropy), ...]}
entropy = calculate_methylation_entropy(control_sites, window_size=200)
```

## Configuration

No dedicated configuration file or environment prefix is used by the methylation module; call the functions above directly from Python.

## Optional Dependencies

- `pandas` -- used by `load_cpg_table`, `compute_beta_values` and `summarize_beta_by_chromosome`
- `scipy` is not used: the DMR test is a z-test implemented with the standard library

## Related Modules

- `metainformant.epigenome.assays.chipseq` -- histone mark analysis
- `metainformant.epigenome.assays.atacseq` -- chromatin accessibility
- `metainformant.epigenome.workflow` -- integrated epigenome pipelines
