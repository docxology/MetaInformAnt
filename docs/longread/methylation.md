# Long-Read Methylation Analysis

The methylation module provides signal-level methylation calling, per-site aggregation, differentially methylated region (DMR) detection, and single-read methylation pattern analysis for nanopore-derived data. It works with both raw signal data and pre-computed methylation tags from BAM files.

## Key Concepts

### Nanopore Methylation Detection

Oxford Nanopore sequencing directly detects DNA base modifications by measuring changes in the electrical current as modified bases pass through the nanopore. The two most commonly detected modifications are:

- **5-methylcytosine (5mC)**: The most studied epigenetic mark, primarily at CpG dinucleotides. Associated with gene silencing and imprinting.
- **N6-methyladenine (6mA)**: Less common in mammals but important in prokaryotes and increasingly studied in eukaryotic gene regulation.

### Signal-Level Calling

Raw electrical signal is analyzed to detect modified bases. The module implements a threshold-based approach on signal features (current level deviations from expected unmodified values) to call modifications at individual positions.

### Methylation Aggregation

Per-read methylation calls are aggregated across multiple reads covering the same genomic position to compute per-site methylation levels (fraction of reads showing the modification).

### Differentially Methylated Regions (DMRs)

DMRs are genomic regions where methylation levels differ significantly between two conditions (e.g., tumor vs. normal, treated vs. untreated). The module implements sliding-window and site-level approaches with statistical testing.

### Pattern Analysis

Single-read methylation patterns reveal co-methylation and allele-specific methylation by examining the methylation state across multiple CpG sites on individual reads.

## Function Reference

### call_methylation_from_signal

```python
def call_methylation_from_signal(
    signal_data: Any,
    reference_positions: list[int],
    modification_type: str = "5mC",
    threshold: float = 0.5,
) -> dict[str, Any]
```

Call methylation from raw nanopore signal. Returns per-position methylation calls with probabilities, including `positions`, `probabilities`, `calls` (binary), and `quality_scores`.

### compute_methylation_stats

```python
def compute_methylation_stats(
    methylation_calls: list[dict[str, Any]],
    region: str | None = None,
) -> dict[str, Any]
```

Compute summary statistics for methylation data. Returns `total_sites`, `mean_methylation`, `median_methylation`, `coverage_distribution`, and per-site statistics.

### aggregate_methylation

```python
def aggregate_methylation(
    per_read_calls: list[dict],
    min_coverage: int = 5,
) -> dict[str, Any]
```

Aggregate per-read methylation calls into per-site methylation levels. Requires minimum coverage threshold. Returns per-site `methylation_level`, `coverage`, `n_methylated`, and `n_unmethylated`.

### detect_dmrs

```python
def detect_dmrs(
    group1_methylation: dict[str, Any],
    group2_methylation: dict[str, Any],
    min_difference: float = 0.2,
    min_sites: int = 3,
    p_value_threshold: float = 0.05,
) -> list[dict[str, Any]]
```

Detect differentially methylated regions between two groups. Returns a list of DMR records with `chrom`, `start`, `end`, `n_sites`, `mean_difference`, `p_value`, and `direction` ("hyper" or "hypo").

### methylation_pattern_analysis

```python
def methylation_pattern_analysis(
    per_read_calls: list[dict],
    region: str | None = None,
    min_cpg_per_read: int = 3,
) -> dict[str, Any]
```

Analyze single-read methylation patterns. Returns `patterns` (unique methylation patterns with frequencies), `co_methylation_matrix`, `allele_specific_signal`, and `entropy` per site.

## Usage Examples

```python
from metainformant.longread import methylation, io

# Signal-level methylation calling
signal = io.extract_signal("reads.fast5", read_id="read_001")
meth_calls = methylation.call_methylation_from_signal(
    signal, reference_positions=[100, 150, 200],
    modification_type="5mC", threshold=0.5,
)

# Aggregate across reads
per_site = methylation.aggregate_methylation(all_read_calls, min_coverage=5)
stats = methylation.compute_methylation_stats(all_read_calls)
print(f"Mean methylation: {stats['mean_methylation']:.2%}")

# Detect DMRs between conditions
dmrs = methylation.detect_dmrs(
    tumor_methylation, normal_methylation,
    min_difference=0.2, min_sites=3,
)
for dmr in dmrs:
    print(f"{dmr['chrom']}:{dmr['start']}-{dmr['end']}: "
          f"delta={dmr['mean_difference']:.2f} (p={dmr['p_value']:.4e})")

# Single-read pattern analysis
patterns = methylation.methylation_pattern_analysis(
    per_read_calls, region="chr1:1000-2000", min_cpg_per_read=3,
)
print(f"Unique patterns: {len(patterns['patterns'])}")
```

## Configuration

- **Environment prefix**: `LR_`
- **Optional dependencies**: numpy, scipy
- Default methylation calling threshold: 0.5 probability
- DMR detection requires minimum 3 CpG sites per region and 20% methylation difference
- Pattern analysis requires at least 3 CpGs per read

## Related Modules

- `longread.io` -- Read FAST5 files and extract signal data
- `longread.analysis.modified_bases` -- BAM-based methylation extraction (from MM/ML tags)
- `longread.phasing` -- Allele-specific methylation analysis
- `longread.visualization` -- `plot_methylation_track` for visualization
