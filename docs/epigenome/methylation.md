# DNA Methylation Analysis

Analysis of DNA methylation patterns, implemented in `metainformant.epigenome.assays.methylation`: CpG site modeling, beta value computation, differentially methylated region (DMR) detection, sequence-based CpG island calling, and methylation entropy.

## Key Concepts

**Beta values** represent the methylation level at each CpG site as the ratio of methylated reads to total reads: beta = M / (M + U), where M and U are the methylated and unmethylated read counts. Values range from 0.0 (unmethylated) to 1.0 (fully methylated); sites with no reads get 0.0. `MethylationSite.methylation_level` is the same ratio for site objects.

**Differentially Methylated Regions (DMRs)** are contiguous runs of CpG sites whose methylation differs between conditions. Sites present in both conditions are compared one by one: a site is significant when its methylation difference reaches a minimum and a two-proportion test on the read counts gives a small p-value (a pooled z-test with Yates' continuity correction when every expected count is at least 10, Fisher's exact test otherwise). Consecutive significant sites are chained into regions of at least `min_sites` sites. There is no sliding window.

**CpG islands** are regions of high CG dinucleotide density in GC-rich sequence. `find_cpg_islands` calls them from DNA sequence with the Gardiner-Garden & Frommer (1987) criteria: length >= 200 bp, GC content >= 50%, and observed/expected CpG ratio >= 0.6, where observed/expected CpG is `CpG count * length / (C count * G count)`. GC content and CpG depletion are properties of the DNA, so methylation levels cannot stand in for sequence: `identify_cpg_islands` also requires the sequence and only summarizes the methylation sites inside each island.

**Methylation entropy** quantifies how heterogeneous methylation is across the CpG sites of a genomic window: the Shannon entropy (in bits) of the distribution of site methylation levels over five equal-width bins, the last of which includes 1.0. Windows whose sites share one level score near 0; windows that mix low, intermediate and high levels score higher.

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

The functions below that take `methylation_data` work on a `Dict[str, List[MethylationSite]]` keyed by chromosome, which the two loaders that follow return.

### `load_methylation_bedgraph(path, min_coverage=1, level_scale="auto") -> Dict[str, List[MethylationSite]]`

Read a BEDgraph file (read as gzip when the name ends in `.gz`) whose column 4 is the methylation level; `#` and `track` lines are skipped. BEDgraph carries levels, not read counts, so every site is given `max(min_coverage, 100)` estimated reads with the methylated count set from the level, and `min_coverage` never removes a site. `level_scale` says how column 4 is written:

- `"fraction"`: levels run from 0 to 1.
- `"percent"`: levels run from 0 to 100.
- `"auto"` (default): the file is read as percentages when any level in it is above 1, and a warning is logged; otherwise it is read as fractions.

A line whose level is outside the range of the scale (negative, above 1 for fractions, above 100 for percentages, or `nan`) is skipped with a warning that names the line, so one bad value does not abort the load. A level that does not parse as a number is skipped as well. Pass `"fraction"` for a fraction file that may hold a stray value above 1, which `"auto"` would read as a percentage file. Any other `level_scale` raises `ValueError`.

### `load_methylation_cov(path, min_coverage=1) -> Dict[str, List[MethylationSite]]`

Read a header-less Bismark `.cov` file (read as gzip when the name ends in `.gz`) with six tab-separated columns: chromosome, start, end, methylation percentage, methylated count and unmethylated count. The site position is `start`, its coverage is methylated plus unmethylated reads (`min_coverage` applies to that total) and its level is computed from the counts; the percentage column only has to parse as a number. A coverage file records no strand, so every site keeps the default `+`. Blank lines, `#` lines and malformed lines (for example an older layout with a strand symbol in column 3) are skipped, the malformed ones with a logged warning.

### `compute_beta_values(df) -> pandas.DataFrame`

Return a copy of `df` (which needs `methylated` and `unmethylated` columns) with a `beta` column: `beta = methylated / (methylated + unmethylated)`, set to 0.0 where both counts are zero. No offset or other regularization is applied.

### `summarize_beta_by_chromosome(df) -> pandas.DataFrame`

Per-chromosome summary of the `beta` column (`df` needs `chrom` and `beta`): a DataFrame indexed by chromosome with the columns `mean`, `std`, `min`, `max` and `count`.

### `find_differentially_methylated_regions(site_data1, site_data2, delta_threshold=0.2, p_value_threshold=0.05, min_sites=3, max_gap=None) -> List[Dict]`

Find DMRs between two conditions, each a `Dict[str, List[MethylationSite]]`. Only positions present in both conditions on a shared chromosome are compared. A position is significant when its absolute methylation difference is at least `delta_threshold` and the two-sided p-value of the test that both conditions share one methylation proportion (from their methylated and total read counts) is at most `p_value_threshold`.

The test is a pooled two-proportion z-test with Yates' continuity correction when every expected count of the methylated/unmethylated by condition table is at least 10, and Fisher's exact test otherwise. The smallest expected count is `min(n1, n2) * min(M, U) / (n1 + n2)`, where `n1` and `n2` are the total reads of the two conditions and `M` and `U` the pooled methylated and unmethylated counts, so few reads or a pooled level close to 0% or 100% selects Fisher's exact test. A site that is fully methylated in one condition and unmethylated in the other is tested like any other site. When every read of both conditions agrees (no methylated reads at all, or no unmethylated ones), the conditions cannot differ and the p-value is 1.

Runs of consecutive significant positions with at least `min_sites` positions become regions; `max_gap` (bp) optionally ends a run when neighboring significant positions are further apart (default `None`: no gap limit).

Returns dicts, ordered by chromosome name and then position, with `chromosome`, `start`, `end`, `sites` (the list of positions), `mean_delta`, `mean_p_value` (the mean of the per-site p-values), `min_p_value` (the smallest of them) and `direction` (`"hyper"` or `"hypo"`: condition 1 relative to condition 2 at the region's first site).

### `find_cpg_islands(sequences, window_size=200, step_size=50, min_gc=0.5, min_cpg_ratio=0.6, min_length=200) -> List[Dict]`

Call CpG islands in DNA sequence. `sequences` maps a chromosome name to its sequence, which is upper-cased; bases other than C and G, including N, count as neither. A chromosome shorter than `window_size` is skipped; the others are scanned with windows of `window_size` bp every `step_size` bp, the last window aligned to the end of the sequence. Windows with a GC fraction of at least `min_gc` and an observed/expected CpG ratio of at least `min_cpg_ratio` are merged when they overlap or touch, and a merged stretch is reported when it is at least `min_length` bp long and still meets both criteria as a whole.

Island boundaries are therefore resolved to about `window_size` bp, because a window that is only partly island still qualifies. Use a smaller `window_size` with `step_size=1` for finer edges. `ValueError` if `window_size` or `step_size` is below 1.

Returns the islands in the order of `sequences` and then by start, as dicts with `chromosome`, `start` and `end` (0-based, end exclusive), `length`, `gc_content`, `cpg_ratio` and `cpg_count`.

### `identify_cpg_islands(methylation_data, window_size=200, step_size=50, min_gc=0.5, min_cpg_ratio=0.6, *, sequences=None, min_length=200) -> List[Dict]`

Call the islands with `find_cpg_islands` and summarize the methylation sites inside each one. `sequences` (a keyword-only mapping of chromosome name to DNA sequence) must be given: the default `None` raises `ValueError`, because GC content and CpG depletion cannot be derived from methylation levels. `methylation_data` may be empty or lack chromosomes, in which case the islands are reported with no sites. Site positions are compared with the 0-based sequence coordinates, so shift 1-based positions (such as Bismark's) first if the difference matters.

Returns the dicts of `find_cpg_islands` plus `sites` (the number of methylation sites with `start <= position < end`) and `avg_methylation` (their mean level, `None` without sites).

### `calculate_methylation_entropy(methylation_data, window_size=1000) -> Dict[str, List[Tuple[int, float]]]`

Shannon entropy (bits) of the site methylation levels in sliding windows of `window_size` bp (step `window_size // 2`), per chromosome; chromosomes with fewer than 10 sites and windows with fewer than 5 sites are skipped. Returns `{chromosome: [(window_center, entropy), ...]}`. Levels are binned into [0, 0.2), [0.2, 0.4), [0.4, 0.6), [0.6, 0.8) and [0.8, 1.0]; the last bin is closed, so fully methylated sites (1.0) are counted. The value is 0 when every site falls in the same bin and reaches log2(5), about 2.32, when sites are spread evenly over the five bins.

## Usage Examples

```python-snippet
from metainformant.epigenome.assays.methylation import (
    MethylationSite,
    calculate_methylation_entropy,
    compute_beta_values,
    find_cpg_islands,
    find_differentially_methylated_regions,
    identify_cpg_islands,
    load_cpg_table,
    load_methylation_bedgraph,
)

# Beta values from a count table (columns: chrom, pos, methylated, unmethylated)
table = load_cpg_table("cpg_counts.tsv")
beta = compute_beta_values(table)
for row in beta.head(3).itertuples():
    print(f"{row.chrom}:{row.pos} beta={row.beta:.3f}")

# A BEDgraph of levels; level_scale="auto" reads the file as percentages when a level exceeds 1
bedgraph_sites = load_methylation_bedgraph("sample.bedGraph", level_scale="auto")


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
        f"mean_delta={dmr['mean_delta']:.3f} min_p={dmr['min_p_value']:.2e} ({dmr['direction']})"
    )

# CpG islands are called from DNA sequence; here a CG repeat at 1000-1300 inside AT-only flanks.
# With window_size=50 and step_size=1 the island comes out as 975-1325.
sequences = {"chr1": "AT" * 500 + "CG" * 150 + "AT" * 500}
islands = find_cpg_islands(sequences, window_size=50, step_size=1)

# The same islands, with the methylation sites inside each one summarized
island_summary = identify_cpg_islands(
    treated_sites, window_size=50, step_size=1, sequences=sequences
)
print(island_summary[0]["sites"], island_summary[0]["avg_methylation"])

# Methylation entropy in 200 bp windows: {chromosome: [(window_center, entropy), ...]}
entropy = calculate_methylation_entropy(control_sites, window_size=200)
```

## Configuration

No dedicated configuration file or environment prefix is used by the methylation module; call the functions above directly from Python.

## Optional Dependencies

- `pandas` -- used by `load_cpg_table`, `compute_beta_values` and `summarize_beta_by_chromosome`
- `scipy` is not used: the DMR tests (a z-test and Fisher's exact test) are implemented with the standard library

## Related Modules

- `metainformant.epigenome.assays.chipseq` -- histone mark analysis
- `metainformant.epigenome.assays.atacseq` -- chromatin accessibility
- `metainformant.epigenome.workflow` -- integrated epigenome pipelines
