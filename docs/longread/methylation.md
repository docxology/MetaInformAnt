# Long-Read Methylation Analysis

The methylation module provides signal-level methylation calling, per-site aggregation, differentially methylated region (DMR) detection, and single-read methylation pattern analysis for nanopore-derived data. It works on plain per-read dictionaries: signal-level records that carry a log-likelihood ratio, and the per-read calls derived from them. The functions live in the `calling` submodule (`from metainformant.longread.methylation import calling`); they are not re-exported on `metainformant.longread.methylation`. For methylation that the basecaller already wrote into BAM `MM`/`ML` tags, see `longread.io.bam.extract_methylation_tags`.

## Key Concepts

### Nanopore Methylation Detection

Oxford Nanopore sequencing directly detects DNA base modifications by measuring changes in the electrical current as modified bases pass through the nanopore. The two most commonly detected modifications are:

- **5-methylcytosine (5mC)**: The most studied epigenetic mark, primarily at CpG dinucleotides. Associated with gene silencing and imprinting.
- **N6-methyladenine (6mA)**: Less common in mammals but important in prokaryotes and increasingly studied in eukaryotic gene regulation.

### Signal-Level Calling

Each per-read, per-site signal record carries a log-likelihood ratio (LLR) for methylation versus no methylation, as produced by a nanopore signal model. The `threshold` model converts the LLR to a probability with the logistic function and calls the site methylated when the probability reaches the threshold, one record at a time. The `bayesian` model pools the reads at one site: they share an unknown methylation rate with a Beta prior, and each read's probability is its posterior given the LLRs of all reads at the site, so a read that disagrees with the rest of its site is pulled toward them. To score raw current traces or per-site signal features directly, see `longread.analysis.modified_bases`.

### Methylation Aggregation

Per-read methylation calls are aggregated across multiple reads covering the same genomic position and strand to compute per-site methylation frequencies (fraction of reads showing the modification).

### Differentially Methylated Regions (DMRs)

DMRs are genomic regions where methylation levels differ significantly between two conditions (e.g., tumor vs. normal, treated vs. untreated). The module segments CpG sites shared by both samples into candidate regions by proximity, then tests each region with a Mann-Whitney U test (a Welch t-test approximation when scipy is unavailable).

### Pattern Analysis

Single-read methylation patterns reveal epigenetic heterogeneity. The module encodes the methylation state of each read across the CpG sites of a region as a string (for example `MMUUM`, with `?` for a CpG the read does not cover), counts the distinct epialleles and reports the Shannon entropy of their distribution.

## Function Reference

All functions live in `metainformant.longread.methylation.calling`.

### call_methylation_from_signal

```python-snippet
def call_methylation_from_signal(
    signal_data: list[dict[str, Any]],
    model: str = "threshold",
    threshold: float = 0.5,
    prior_alpha: float = 1.0,
    prior_beta: float = 1.0,
) -> list[dict[str, Any]]
```

Call methylation status for per-read, per-site signal records. Each record in `signal_data` carries `read_id`, `chrom`, `position`, `strand` and `log_likelihood_ratio` (positive favours methylation). Returns a list of the same records, each copied with two added fields: `methylated` (bool; true when the probability reaches `threshold`) and `confidence` (0-1: twice the distance of the probability from 0.5, so 0 at a probability of 0.5 and 1 at 0 or 1). `model` is `"threshold"` or `"bayesian"`; any other value raises `ValueError`. A missing `log_likelihood_ratio` raises `ValueError`; it is not treated as neutral evidence. A record whose LLR is NaN is called unmethylated with a NaN `confidence`, and under `"bayesian"` it does not contribute to its site.

- `"threshold"`: the probability of a record is the logistic function of its LLR, `1 / (1 + exp(-LLR))`, the posterior under a flat 50 % prior. `prior_alpha` and `prior_beta` are ignored.
- `"bayesian"`: a hierarchical Bayesian call. Records with the same `(chrom, position, strand)` share an unknown methylation rate `theta ~ Beta(prior_alpha, prior_beta)`; each read is methylated with probability `theta`, and its LLR is the likelihood ratio of its signal. The probability of a read is its posterior `P(methylated | LLRs of all reads at the site)`, the expectation of `logistic(LLR + logit(theta))` over the posterior of `theta`, computed by numerical integration (accurate to better than 1e-6). A read that disagrees with the rest of its site is pulled toward it, and the result does not depend on the order of the reads. A site with one read has probability `logistic(LLR + logit(prior mean))`, where the prior mean is `prior_alpha / (prior_alpha + prior_beta)`; with the default Beta(1, 1) prior that equals the threshold model. Records without both `chrom` and `position` are not pooled: each is its own site. `prior_alpha` and `prior_beta` must each be at least 0.25, otherwise `ValueError`. The model requires numpy and raises `ImportError` without it.

### aggregate_methylation

```python-snippet
def aggregate_methylation(
    calls: list[dict[str, Any]],
    min_coverage: int = 5,
) -> dict[str, Any]
```

Aggregate per-read calls (the output of `call_methylation_from_signal`, which need `chrom`, `position`, `strand` and `methylated`) into per-site methylation frequencies. Sites with fewer than `min_coverage` reads are dropped, and `min_coverage < 1` raises `ValueError`. Every call must supply non-null `position` and `methylated` fields; missing values raise `ValueError`. Returns a dictionary with `sites` (a list of dicts holding `chrom`, `position`, `strand`, `methylation_freq`, `coverage`, `n_methylated` and `n_unmethylated`), `n_sites` and `mean_coverage`.

### compute_methylation_stats

```python-snippet
def compute_methylation_stats(
    sites: list[dict[str, Any]],
) -> dict[str, Any]
```

Compute summary statistics for the `sites` list returned by `aggregate_methylation` (each site needs `methylation_freq`; `coverage` and `context` are optional). Returns `global_methylation`, `n_sites`, `mean_coverage`, `methylation_distribution` (`min`, `max`, `median`, `std`, `q25`, `q75`), `context_proportions`, `highly_methylated_fraction` (frequency of at least 0.8) and `lowly_methylated_fraction` (at most 0.2). Quartiles and the median use linearly interpolated percentiles. An empty list returns the same keys with zero or empty values.

### detect_dmrs

```python-snippet
def detect_dmrs(
    methylation_a: dict[str, Any],
    methylation_b: dict[str, Any],
    min_cpgs: int = 3,
    min_diff: float = 0.2,
    max_distance: int = 1000,
) -> list[dict[str, Any]]
```

Detect differentially methylated regions between two samples. Both arguments are `aggregate_methylation` results (the whole dictionary, with its `sites` list). Only positions present in both samples are used; consecutive shared CpGs no farther apart than `max_distance` form a candidate region, and regions with at least `min_cpgs` sites and an absolute mean difference of at least `min_diff` are tested. Returns a list of DMR dicts with `chrom`, `start`, `end`, `n_cpgs`, `mean_diff` (B minus A), `p_value` and `direction` (`"hyper"` if B is higher, otherwise `"hypo"`). An empty sample or no shared CpG positions returns an empty list with a warning; the p-value is not adjusted for multiple testing.

### methylation_pattern_analysis

```python-snippet
def methylation_pattern_analysis(
    calls: list[dict[str, Any]],
    region: dict[str, Any],
) -> dict[str, Any]
```

Analyze single-read methylation patterns inside a region. `calls` are per-read call dicts (`read_id`, `chrom`, `position`, `methylated`) and `region` is a dict with non-null `chrom`, `start` and `end`; an incomplete region raises `ValueError`. A read is counted only when it covers at least `max(1, n // 2)` of the `n` CpG positions seen in the region (so a single site is enough when `n` is below 4). Returns `patterns` (the distinct epiallele strings over `M`, `U` and `?` for an uncovered CpG, most frequent first), `entropy` (Shannon entropy in bits), `n_reads` and `epialleles` (pattern string to count).

## Usage Examples

```python
import random

from metainformant.longread.methylation import calling

random.seed(1)
positions = [1000 + 40 * i for i in range(8)]


def simulate_signal(mean_llr, n_reads=8):
    """Signal-level records: one per read and CpG site (log_likelihood_ratio > 0 favours methylation)."""
    return [
        {"read_id": f"read_{r:03d}", "chrom": "chr1", "position": pos, "strand": "+",
         "log_likelihood_ratio": random.gauss(mean_llr, 1.5)}
        for r in range(n_reads)
        for pos in positions
    ]


# Signal-level methylation calling (replace the simulated records with your own)
tumor_calls = calling.call_methylation_from_signal(simulate_signal(2.0), model="threshold", threshold=0.5)
normal_calls = calling.call_methylation_from_signal(simulate_signal(-2.0), model="threshold", threshold=0.5)

# Bayesian model (needs numpy): reads at one site share a methylation rate with a Beta prior, so a weak
# dissenting read is pulled toward its site. The prior defaults to Beta(1, 1).
site = [
    {"read_id": f"read_{i}", "chrom": "chr1", "position": 1000, "strand": "+", "log_likelihood_ratio": llr}
    for i, llr in enumerate([3.0, 3.0, 3.0, 3.0, 3.0, -1.0])
]
alone = calling.call_methylation_from_signal(site, model="threshold")
pooled = calling.call_methylation_from_signal(site, model="bayesian", prior_alpha=1.0, prior_beta=1.0)
print(alone[-1]["methylated"], pooled[-1]["methylated"])  # False True

# Aggregate across reads
tumor_methylation = calling.aggregate_methylation(tumor_calls, min_coverage=5)
normal_methylation = calling.aggregate_methylation(normal_calls, min_coverage=5)
print(f"{tumor_methylation['n_sites']} sites, mean coverage {tumor_methylation['mean_coverage']:.1f}")
stats = calling.compute_methylation_stats(tumor_methylation["sites"])
print(f"Global methylation: {stats['global_methylation']:.2%}")

# Detect DMRs between conditions (both arguments are aggregate_methylation results)
dmrs = calling.detect_dmrs(tumor_methylation, normal_methylation, min_cpgs=3, min_diff=0.2)
for dmr in dmrs:
    print(f"{dmr['chrom']}:{dmr['start']}-{dmr['end']}: "
          f"delta={dmr['mean_diff']:.2f} (p={dmr['p_value']:.4e}, {dmr['direction']})")

# Single-read pattern analysis
patterns = calling.methylation_pattern_analysis(
    tumor_calls, region={"chrom": "chr1", "start": 1000, "end": 1400},
)
print(f"Unique patterns: {len(patterns['patterns'])}, entropy: {patterns['entropy']:.2f} bits")
```

## Configuration

- **Environment prefix**: `LR_` is honoured only by `workflow.pipelines.load_pipeline_config` and `utils.batch.process_batch` (`LR_BATCH_WORKERS`); the functions on this page read no environment variables
- **Optional dependencies**: numpy (required by `model="bayesian"` in `call_methylation_from_signal`; the other functions do not use it), scipy (the Mann-Whitney test in `detect_dmrs` falls back to a Welch t-test approximation without scipy)
- Default methylation calling threshold: 0.5 probability
- Bayesian model prior: Beta(1, 1) by default (`prior_alpha` and `prior_beta`, each at least 0.25)
- Aggregation requires at least 5 reads per site by default
- DMR detection requires minimum 3 CpG sites per region and 20% methylation difference
- Pattern analysis counts a read only if it covers at least `max(1, n // 2)` of the region's `n` CpG sites

## Related Modules

- `longread.io` -- `fast5` reads signal files; `bam.extract_methylation_tags` decodes MM/ML tags
- `longread.analysis.modified_bases` -- Calling from raw current traces and per-site signal features (`detect_methylation`, `call_5mc`, `call_6ma`), region aggregation and differential methylation
- `longread.phasing` -- Allele-specific methylation analysis
- `longread.visualization` -- `plot_methylation_track` for visualization
