# Long-Read Quality Control

The quality module provides comprehensive quality metrics, read filtering, adapter detection, and chimeric read handling for long-read sequencing data. It computes N50/Nx statistics, quality score distributions, accuracy estimates, and throughput calculations.

The package has two submodules: `metrics` (N50/Nx, length and quality statistics, accuracy, throughput) and `filtering` (read filters, adapter detection and trimming, chimeric read splitting). Import them with `from metainformant.longread.quality import filtering, metrics`; the functions are not re-exported on `metainformant.longread.quality` itself.

## Key Concepts

### N50 and Nx Statistics

N50 is the read length such that 50% of the total sequenced bases are in reads at least this long. It is the standard metric for characterizing read length distributions. Nx generalizes this to any percentage x. L50 is the minimum number of reads needed to cover 50% of total bases.

### Quality Score Distribution

Phred quality scores (Q values) indicate per-base accuracy: Q10 = 90% accuracy, Q15 = 96.8%, Q20 = 99%, Q30 = 99.9%. The module tracks the fraction of reads meeting various quality thresholds (Q7, Q10, Q15, Q20).

### Adapter Detection and Trimming

ONT and PacBio reads may contain adapter sequences at their ends. The module detects known adapter sequences (LSK109, LSK110, rapid, RNA adapters for ONT; SMRTbell for PacBio) and trims them. Matching compares each adapter, and its reverse complement, with windows of the read and accepts windows above an identity threshold. That tolerates substitution errors but not insertions or deletions inside the adapter. `trim_adapters` tests every offset near the read ends, whereas `detect_adapters` (and therefore `split_chimeric_reads`) steps through the whole read in strides of half the adapter length, so an adapter at an unlucky offset can be missed.

### Chimeric Read Splitting

Chimeric reads (artifacts containing segments from different genomic locations) are detected by finding adapter sequences in the interior of the read, then split into their constituent segments.

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

### AdapterMatch

```python
@dataclass
class AdapterMatch:
    adapter_name: str           # Name in the adapter table; "_rc" is appended for a reverse-complement hit
    adapter_sequence: str       # Sequence that matched
    position: int               # 0-based start in the read
    end_position: int           # Exclusive end in the read
    score: float                # identity * adapter length
    identity: float             # Fraction of matching bases
    location: str               # "start", "end" or "internal"
```

## Function Reference

### Quality Metrics

These functions live in `metainformant.longread.quality.metrics`.

```python-snippet
def calculate_n50(read_lengths: Sequence[int]) -> int
def calculate_nx(read_lengths: Sequence[int], x: int = 50) -> int
def read_length_stats(reads: Sequence[dict | int | object]) -> ReadLengthStatistics
def quality_score_distribution(reads: Sequence[dict | str | object]) -> QualityDistribution
def estimate_accuracy(quality_scores: Sequence[int] | str) -> float
def calculate_throughput(
    reads: Sequence[dict | int | object], run_duration: float | None = None,
) -> dict[str, float]
```

`calculate_n50` computes the N50 statistic and `calculate_nx` generalizes it to any Nx (`x` between 0 and 100, otherwise `ValueError`); both return 0 for empty input. `read_length_stats` provides comprehensive length statistics; it accepts Python integer lengths (convert a numpy array with `.tolist()` first; numpy integers are not recognized), dicts with a `length` or `sequence` key, or objects with a `sequence` (or `query_length`) attribute. `quality_score_distribution` analyzes per-read mean Phred quality from Phred+33 strings, dicts with a `quality` or `quality_string` key, or objects with a `quality_string` attribute. `estimate_accuracy` takes integer Phred scores or a Phred+33 string and returns one number, the mean base accuracy as a fraction (`1 - mean(10^(-Q/10))`). `calculate_throughput` always returns `total_bases`, `total_reads` and `mean_read_length`; when `run_duration` (in hours) is given it also returns `bases_per_hour`, `reads_per_hour` and `gigabases_per_hour`.

### Read Filtering

These functions live in `metainformant.longread.quality.filtering`.

```python-snippet
def filter_by_length(
    reads: Sequence[T], min_length: int = 1000, max_length: int | None = None,
) -> list[T]

def filter_by_quality(
    reads: Sequence[T], min_q: float = 7.0,
) -> list[T]
```

Filter reads by length range or minimum mean quality score. Reads can be `ReadRecord` objects, dicts with a `sequence` (and `quality` or `quality_string`) key, or any object with those attributes; the same objects are returned. `filter_by_length` also drops reads that have no sequence, while `filter_by_quality` keeps reads that carry no quality string.

### Adapter Handling

```python-snippet
def detect_adapters(
    sequence: str,
    known_adapters: dict[str, str] | None = None,
    max_search_length: int = 300,
    min_identity: float = 0.70,
) -> list[AdapterMatch]

def trim_adapters(
    reads: Sequence[dict | ReadRecord],
    adapter_sequences: dict[str, str] | None = None,
    max_search_length: int = 200,
    min_identity: float = 0.75,
) -> list[ReadRecord]

def split_chimeric_reads(
    reads: Sequence[dict | ReadRecord],
    known_adapters: dict[str, str] | None = None,
    min_fragment_length: int = 500,
    min_adapter_identity: float = 0.75,
) -> list[ReadRecord]
```

`detect_adapters` scans one read sequence (a string, not a list of reads) for adapters anywhere along its length and returns `AdapterMatch` objects labelled `start`, `end` or `internal`; the built-in ONT and PacBio adapter tables (`ALL_KNOWN_ADAPTERS`, also available as `ONT_ADAPTERS` and `PACBIO_ADAPTERS`) are used when no dictionary is given. `trim_adapters` searches the first `max_search_length` bases for adapters in the forward orientation and the last `max_search_length` bases in both orientations, trims what it finds, and returns new `ReadRecord`s whose `metadata` records `original_length`, `trimmed_start` and `trimmed_end`; reads that nothing remains of are dropped. `split_chimeric_reads` splits reads containing internal adapter sequences into fragments named `<read_id>_fragment_<n>` (with `chimeric`, `parent_read` and fragment coordinates in `metadata`), drops fragments shorter than `min_fragment_length`, and returns every other read as a `ReadRecord` with `metadata["chimeric"] = False`.

**Caution:** the built-in table contains two 8 bp LSK110 sequences. At the default identity thresholds (6 of 8 bases is enough) they match random sequence, so with the default table `trim_adapters` can trim real bases from adapter-free reads and `detect_adapters` / `split_chimeric_reads` can report spurious `internal` hits. For real data pass an explicit adapter dictionary, for example only the entries of `ALL_KNOWN_ADAPTERS` that are at least 20 bp long, or raise `min_identity`.

## Usage Examples

```python
from metainformant.longread.io import fast5
from metainformant.longread.quality import filtering, metrics

# Load reads (skip any without basecalls)
reads = [r for r in fast5.read_fast5("reads.fast5") if r.sequence]
lengths = [len(r.sequence) for r in reads]

# Compute N50 and comprehensive stats
n50 = metrics.calculate_n50(lengths)
stats = metrics.read_length_stats(lengths)
print(f"N50: {stats.n50:,} bp, Mean: {stats.mean_length:,.0f} bp")
print(f"Total bases: {stats.total_bases:,}, Reads: {stats.count:,}")

# Quality score analysis
qual_strings = [r.quality_string for r in reads if r.quality_string]
qual_dist = metrics.quality_score_distribution(qual_strings)
print(f"Mean Q: {qual_dist.mean_quality:.1f}, Q20+: {qual_dist.q20_fraction:.1%}")
print(f"Estimated accuracy: {metrics.estimate_accuracy(qual_strings[0]):.2%}")

# Throughput
throughput = metrics.calculate_throughput(lengths, run_duration=24.0)
print(f"Throughput: {throughput['bases_per_hour']:,.0f} bases/hour")

# Filter reads
read_records = [filtering.ReadRecord(r.read_id, r.sequence, r.quality_string) for r in reads]
filtered = filtering.filter_by_length(read_records, min_length=1000)
filtered = filtering.filter_by_quality(filtered, min_q=10.0)
print(f"After filtering: {len(filtered)} / {len(read_records)} reads")

# Detect and trim adapters (leave out the 8 bp LSK110 entries, see the caution above)
adapters = {name: seq for name, seq in filtering.ALL_KNOWN_ADAPTERS.items() if len(seq) >= 20}
for match in filtering.detect_adapters(read_records[0].sequence, adapters):
    print(f"{match.adapter_name}: {match.location} at {match.position}-{match.end_position}")
trimmed = filtering.trim_adapters(read_records, adapters)

# Split chimeric reads
split_reads = filtering.split_chimeric_reads(read_records, adapters)
print(f"Split: {len(read_records)} -> {len(split_reads)} reads")
```

## Configuration

- **Environment prefix**: `LR_` is honoured only by `workflow.pipelines.load_pipeline_config` and `utils.batch.process_batch` (`LR_BATCH_WORKERS`); the functions on this page read no environment variables
- **Optional dependencies**: none in use (the modules try to import numpy but never use it)
- Built-in adapter databases include ONT (LSK109, LSK110, rapid, RNA) and PacBio (SMRTbell) adapters
- Default quality filter threshold: Q7 (ONT standard); default minimum read length: 1000 bp
- Adapter matching is identity-based, not mismatch-count based: 0.70 minimum identity for `detect_adapters`, 0.75 for `trim_adapters` and `split_chimeric_reads`

## Related Modules

- `longread.io` -- Read FAST5/BAM files for quality analysis
- `longread.assembly` -- Assembly from the filtered reads
- `longread.utils` -- Batch processing and QC summary generation
- `longread.visualization` -- `plot_read_length_histogram`, `plot_quality_vs_length`
