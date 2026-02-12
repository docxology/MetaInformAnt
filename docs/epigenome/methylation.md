# DNA Methylation Analysis

Analysis of DNA methylation patterns including CpG site modeling, beta value computation, differential methylation region (DMR) detection, CpG island identification, and methylation entropy.

## Key Concepts

**Beta values** represent the methylation level at each CpG site as the ratio of methylated signal to total signal: beta = M / (M + U + offset). Values range from 0.0 (unmethylated) to 1.0 (fully methylated).

**Differentially Methylated Regions (DMRs)** are contiguous genomic intervals where methylation differs significantly between conditions, identified by sliding-window analysis with statistical testing.

**CpG islands** are regions of high CG dinucleotide density, typically defined by observed/expected CpG ratio >= 0.6, GC content >= 50%, and length >= 200 bp.

**Methylation entropy** quantifies the heterogeneity of methylation patterns across reads at a locus, distinguishing uniform from stochastic methylation states.

## Data Model

### `MethylationSite`

Dataclass representing a single CpG site measurement.

| Field | Type | Description |
|-------|------|-------------|
| `chrom` | `str` | Chromosome name |
| `position` | `int` | Genomic coordinate |
| `methylated_count` | `int` | Number of methylated reads |
| `unmethylated_count` | `int` | Number of unmethylated reads |
| `beta_value` | `float` | Methylation level (0.0--1.0) |
| `context` | `str` | Sequence context (CG, CHG, CHH) |

## Function Reference

### `load_cpg_table(file_path, format="bismark") -> List[MethylationSite]`

Load CpG methylation data from tab-delimited files. Supports Bismark coverage format (chrom, start, end, percent, count_M, count_U) and custom formats.

### `compute_beta_values(sites, offset=100) -> List[MethylationSite]`

Calculate beta values for each site: `beta = M / (M + U + offset)`. The offset (default 100) regularizes estimates from low-coverage sites.

### `summarize_beta_by_chromosome(sites) -> Dict[str, Dict]`

Per-chromosome summary statistics: mean beta, median beta, standard deviation, and number of sites.

### `find_dmrs(group1_sites, group2_sites, window_size=1000, step_size=500, min_diff=0.2, min_sites=3, p_threshold=0.05) -> List[Dict]`

Identify differentially methylated regions between two groups using a sliding-window approach with Welch's t-test (scipy) or Mann-Whitney U fallback.

Returns dicts with `chrom`, `start`, `end`, `mean_diff`, `p_value`, `n_sites`.

### `find_cpg_islands(sequence, chrom="chr1", min_length=200, min_gc=0.5, min_obs_exp=0.6) -> List[Dict]`

Scan a DNA sequence for CpG islands using the Gardiner-Garden and Frommer (1987) criteria: minimum length, GC fraction, and observed/expected CpG ratio.

### `calculate_methylation_entropy(sites, window_size=4) -> List[Dict]`

Compute Shannon entropy of methylation patterns within sliding windows. High entropy indicates heterogeneous (stochastic) methylation; low entropy indicates uniform patterns.

## Usage Examples

```python
from metainformant.epigenome import (
    MethylationSite,
    compute_beta_values,
    find_dmrs,
    find_cpg_islands,
    calculate_methylation_entropy,
)

# Create methylation sites
sites = [
    MethylationSite("chr1", 1000, 80, 20),
    MethylationSite("chr1", 1050, 10, 90),
    MethylationSite("chr1", 1100, 45, 55),
]

# Compute beta values with regularization
sites = compute_beta_values(sites, offset=100)
for s in sites:
    print(f"{s.chrom}:{s.position} beta={s.beta_value:.3f}")

# Find DMRs between conditions
dmrs = find_dmrs(control_sites, treated_sites, min_diff=0.2)
for dmr in dmrs:
    print(f"{dmr['chrom']}:{dmr['start']}-{dmr['end']} diff={dmr['mean_diff']:.3f}")

# Identify CpG islands in a sequence
islands = find_cpg_islands("ACGTCGCGCGATCGCGCGATCG" * 50)

# Methylation entropy
entropy = calculate_methylation_entropy(sites, window_size=4)
```

## Configuration

Environment variable prefix: `EPI_`

## Optional Dependencies

- `scipy` -- statistical tests for DMR detection (Welch's t-test, Mann-Whitney U)
- Pure Python fallbacks are used when scipy is unavailable

## Related Modules

- `metainformant.epigenome.chipseq` -- histone mark analysis
- `metainformant.epigenome.atacseq` -- chromatin accessibility
- `metainformant.epigenome.workflow` -- integrated epigenome pipelines
