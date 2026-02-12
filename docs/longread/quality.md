# Long-Read Quality Control

The quality module provides comprehensive quality metrics, read filtering, adapter detection, and chimeric read handling for long-read sequencing data. It computes N50/Nx statistics, quality score distributions, accuracy estimates, and throughput calculations.

## Key Concepts

### N50 and Nx Statistics

N50 is the read length such that 50% of the total sequenced bases are in reads at least this long. It is the standard metric for characterizing read length distributions. Nx generalizes this to any percentage x. L50 is the minimum number of reads needed to cover 50% of total bases.

### Quality Score Distribution

Phred quality scores (Q values) indicate per-base accuracy: Q10 = 90% accuracy, Q15 = 96.8%, Q20 = 99%, Q30 = 99.9%. The module tracks the fraction of reads meeting various quality thresholds (Q7, Q10, Q15, Q20).

### Adapter Detection and Trimming

ONT and PacBio reads may contain adapter sequences at their ends. The module detects known adapter sequences (LSK109, LSK110, rapid, RNA adapters for ONT; SMRTbell for PacBio) and trims them. Detection uses approximate string matching to handle sequencing errors.

### Chimeric Read Splitting

Chimeric reads (artifacts containing segments from different genomic locations) are detected by analyzing internal adapter sequences or split-alignment patterns, then split into their constituent segments.

## Data Structures

### ReadLengthStatistics

```python
@dataclass
class ReadLengthStatistics:
    count: int; total_bases: int; mean_length: float
    median_length: float; min_length: int; max_length: int
    std_dev: float; n50: int; n90: int; l50: int
    percentiles: dict[str, float]  # p5, p25, p75, p95
```

### QualityDistribution

```python
@dataclass
class QualityDistribution:
    mean_quality: float; median_quality: float
    min_quality: float; max_quality: float
    q7_fraction: float; q10_fraction: float
    q15_fraction: float; q20_fraction: float
    per_read_means: list[float]
    histogram: dict[int, int]
```

### ReadRecord

```python
@dataclass
class ReadRecord:
    read_id: str
    sequence: str
    quality_string: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)
```

## Function Reference

### Quality Metrics

```python
def calculate_n50(read_lengths: list[int]) -> int
def calculate_nx(read_lengths: list[int], x: float = 50.0) -> int
def read_length_stats(read_lengths: list[int]) -> ReadLengthStatistics
def quality_score_distribution(quality_strings: list[str]) -> QualityDistribution
def estimate_accuracy(quality_scores: list[float]) -> dict[str, float]
def calculate_throughput(
    read_lengths: list[int], run_duration_hours: float,
) -> dict[str, float]
```

`calculate_n50` computes the N50 statistic. `calculate_nx` generalizes to any Nx. `read_length_stats` provides comprehensive length statistics. `quality_score_distribution` analyzes Phred quality score distributions. `estimate_accuracy` converts quality scores to estimated accuracy. `calculate_throughput` computes bases per hour and reads per hour.

### Read Filtering

```python
def filter_by_length(
    reads: list[ReadRecord],
    min_length: int = 0, max_length: int | None = None,
) -> list[ReadRecord]

def filter_by_quality(
    reads: list[ReadRecord],
    min_mean_quality: float = 7.0,
) -> list[ReadRecord]
```

Filter reads by length range or minimum mean quality score.

### Adapter Handling

```python
def detect_adapters(
    reads: list[ReadRecord],
    adapter_sequences: dict[str, str] | None = None,
    max_mismatches: int = 3,
) -> list[dict[str, Any]]

def trim_adapters(
    reads: list[ReadRecord],
    adapter_sequences: dict[str, str] | None = None,
) -> list[ReadRecord]

def split_chimeric_reads(
    reads: list[ReadRecord],
    adapter_sequences: dict[str, str] | None = None,
) -> list[ReadRecord]
```

`detect_adapters` identifies adapter sequences at read ends (uses built-in ONT and PacBio adapter databases if none provided). `trim_adapters` removes detected adapters. `split_chimeric_reads` splits reads containing internal adapter sequences into separate segments.

## Usage Examples

```python
from metainformant.longread import quality, io

# Load reads
reads = io.read_fast5("reads.fast5")
lengths = [len(r.sequence) for r in reads]

# Compute N50 and comprehensive stats
n50 = quality.calculate_n50(lengths)
stats = quality.read_length_stats(lengths)
print(f"N50: {stats.n50:,} bp, Mean: {stats.mean_length:,.0f} bp")
print(f"Total bases: {stats.total_bases:,}, Reads: {stats.count:,}")

# Quality score analysis
qual_strings = [r.quality_string for r in reads if r.quality_string]
qual_dist = quality.quality_score_distribution(qual_strings)
print(f"Mean Q: {qual_dist.mean_quality:.1f}, Q20+: {qual_dist.q20_fraction:.1%}")

# Throughput
throughput = quality.calculate_throughput(lengths, run_duration_hours=24.0)
print(f"Throughput: {throughput['bases_per_hour']:,.0f} bases/hour")

# Filter reads
read_records = [quality.ReadRecord(r.read_id, r.sequence, r.quality_string) for r in reads]
filtered = quality.filter_by_length(read_records, min_length=1000)
filtered = quality.filter_by_quality(filtered, min_mean_quality=10.0)
print(f"After filtering: {len(filtered)} / {len(read_records)} reads")

# Detect and trim adapters
adapter_report = quality.detect_adapters(read_records)
trimmed = quality.trim_adapters(read_records)

# Split chimeric reads
split_reads = quality.split_chimeric_reads(read_records)
print(f"Split: {len(read_records)} -> {len(split_reads)} reads")
```

## Configuration

- **Environment prefix**: `LR_`
- **Optional dependencies**: numpy
- Built-in adapter databases include ONT (LSK109, LSK110, rapid, RNA) and PacBio (SMRTbell) adapters
- Default quality filter threshold: Q7 (ONT standard)
- Adapter detection allows up to 3 mismatches by default

## Related Modules

- `longread.io` -- Read FAST5/BAM files for quality analysis
- `longread.assembly` -- Quality filtering before assembly
- `longread.utils` -- Batch processing and QC summary generation
- `longread.visualization` -- `plot_read_length_histogram`, `plot_quality_vs_length`
