# Long-Read Quality Control

The quality module provides comprehensive quality metrics, read filtering, adapter detection, and chimeric read handling for long-read sequencing data. It computes N50/Nx statistics, quality score distributions, accuracy estimates, and throughput calculations.

The package has two submodules: `metrics` (N50/Nx, length and quality statistics, accuracy, throughput) and `filtering` (read filters, adapter detection and trimming, chimeric read splitting). Import them with `from metainformant.longread.quality import filtering, metrics`; the functions are not re-exported on `metainformant.longread.quality` itself.

## Key Concepts

### N50 and Nx Statistics

N50 is the read length such that 50% of the total sequenced bases are in reads at least this long. It is the standard metric for characterizing read length distributions. Nx generalizes this to any percentage x. L50 is the minimum number of reads needed to cover 50% of total bases.

### Quality Score Distribution

Phred quality scores (Q values) indicate per-base accuracy: Q10 = 90% accuracy, Q15 = 96.8%, Q20 = 99%, Q30 = 99.9%. The module tracks the fraction of reads meeting various quality thresholds (Q7, Q10, Q15, Q20).

### Adapter Detection and Trimming

ONT and PacBio reads may contain adapter sequences at their ends, and chimeric reads carry them in the interior. The module detects known adapters (LSK109, rapid and RNA adapters for ONT; SMRTbell and barcoded adapters for PacBio) in both orientations, then trims or splits at them.

An adapter is located by semi-global alignment: the whole adapter must align to a stretch of the read, the stretch may start and end anywhere, and each substitution, insertion or deletion costs one edit. A placement is accepted when it needs at most `max_errors` edits, the smallest of three limits:

- the divergence the caller tolerates, `floor(len(adapter) * (1 - min_identity))`;
- the largest budget whose false-positive bound stays within `max_false_positive_rate` (1e-9 per read offset by default);
- `len(adapter) // 7 - 1`, which keeps the exact seeds that propose candidate placements at least 7 bp long.

An adapter with budget `k` is cut into `k + 1` exact seeds, and every placement within the budget contains one of them unchanged (pigeonhole). The seeds propose candidate regions and the alignment confirms them, so the search loses no placement the model accepts and finds an adapter at any offset in the read.

### False-Positive Bound

For i.i.d. uniform bases, the probability that a given read offset ends a spurious placement of an `L` bp adapter within `k` edits is at most `adapter_false_positive_bound(L, k)`. It is a union bound over edit scripts, so it is conservative. A read of `n` adapter-free bases has at most `2 * n * bound` expected spurious hits per adapter, counting both orientations; at the default rate that is at most 2e-5 hits for a 10 kb read. The bound assumes uniform base composition: GC- or AT-rich sequence raises the true rate, so lower `max_false_positive_rate` for such genomes.

The shortest adapter whose exact match stays within a given rate is `minimum_adapter_length(rate)`, 15 bp at the default. An adapter shorter than that is searched by exact match only, and a warning names it. The edit budget grows with adapter length: at the default rate and `min_identity` 0.75, an adapter of 15 to 18 bp is matched exactly, 19 to 21 bp tolerates one edit, 22 to 27 bp two, and 28 to 34 bp three.

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
    adapter_sequence: str       # Sequence that matched (the reverse complement for an "_rc" hit)
    position: int               # 0-based start in the read
    end_position: int           # Exclusive end in the read
    score: float                # Adapter length minus edit_distance
    identity: float             # 1 - edit_distance / adapter length
    location: str               # "start", "end" or "internal"
    edit_distance: int = 0      # Substitutions, insertions and deletions in the placement
    false_positive_bound: float = 0.0  # Per-offset bound for this adapter's length and edit budget
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

`calculate_n50` computes the N50 statistic and `calculate_nx` generalizes it to any Nx (`x` between 0 and 100, otherwise `ValueError`); both return 0 for empty input. `read_length_stats` provides comprehensive length statistics; it accepts integer lengths (Python integers or numpy integer scalars, so an integer numpy array can be passed directly), dicts with a `length` or `sequence` key, or objects with a `sequence` (or `query_length`) attribute. Whole-number floats such as `1000.0` (including float numpy arrays) are accepted; fractional, NaN and infinite lengths raise `ValueError` instead of being skipped. `quality_score_distribution` analyzes per-read mean Phred quality from Phred+33 strings, dicts with a `quality` or `quality_string` key, or objects with a `quality_string` attribute. `estimate_accuracy` takes integer Phred scores or a Phred+33 string and returns one number, the mean base accuracy as a fraction (`1 - mean(10^(-Q/10))`). `calculate_throughput` takes the same read forms as `read_length_stats` (integer lengths, numpy integers included, dicts or objects) and applies the same whole-number validation to float lengths. It always returns `total_bases`, `total_reads` and `mean_read_length`; when `run_duration` (in hours) is given it also returns `bases_per_hour`, `reads_per_hour` and `gigabases_per_hour`. `plot_read_length_histogram` in `longread.visualization` also takes integer lengths, numpy integers included.

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
    max_false_positive_rate: float = 1e-9,
) -> list[AdapterMatch]

def trim_adapters(
    reads: Sequence[dict | ReadRecord],
    adapter_sequences: dict[str, str] | None = None,
    max_search_length: int = 200,
    min_identity: float = 0.75,
    max_false_positive_rate: float = 1e-9,
) -> list[ReadRecord]

def split_chimeric_reads(
    reads: Sequence[dict | ReadRecord],
    known_adapters: dict[str, str] | None = None,
    min_fragment_length: int = 500,
    min_adapter_identity: float = 0.75,
    max_false_positive_rate: float = 1e-9,
) -> list[ReadRecord]

def adapter_false_positive_bound(adapter_length: int, max_errors: int) -> float
def adapter_error_budget(
    adapter_length: int, min_identity: float = 0.75, max_false_positive_rate: float = 1e-9,
) -> int
def minimum_adapter_length(max_false_positive_rate: float = 1e-9) -> int
```

`detect_adapters` scans one read sequence (a string, not a list of reads) for adapters at every offset and returns `AdapterMatch` objects sorted by position, with overlapping matches resolved by the lowest edit rate per adapter base, then by the longer adapter at equal rates. A match is labelled `start` when it begins before position `max_search_length`, `end` when it begins after `len(sequence) - max_search_length`, and `internal` otherwise. The built-in ONT, PacBio and Illumina adapter tables (`ALL_KNOWN_ADAPTERS`, also available as `ONT_ADAPTERS`, `PACBIO_ADAPTERS` and `ILLUMINA_ADAPTERS`) are used when no dictionary is given. `trim_adapters` searches the first `max_search_length` bases for adapters in the forward orientation and the last `max_search_length` bases in both orientations, trims up to the end of an adapter found at the start and from the start of one found at the end, and returns new `ReadRecord`s whose `metadata` records `original_length` and the number of bases removed at each end (`trimmed_start`, `trimmed_end`); reads that nothing remains of are dropped. `split_chimeric_reads` splits reads containing internal adapter sequences into fragments named `<read_id>_fragment_<n>` (with `chimeric`, `parent_read` and fragment coordinates in `metadata`), drops fragments shorter than `min_fragment_length`, and returns every other read as a `ReadRecord` with `metadata["chimeric"] = False`.

`min_identity` (`min_adapter_identity` in `split_chimeric_reads`) is the minimum fraction of the adapter that must match, from 0 to 1. It caps the tolerated edits at `floor(len(adapter) * (1 - min_identity))`; the false-positive bound and the seed limit can lower the cap but never raise it, so identity 1.0 allows no edits. `max_false_positive_rate` (default 1e-9, also available as `DEFAULT_MAX_FALSE_POSITIVE_RATE`) is the upper bound on the per-offset probability of a spurious placement in adapter-free i.i.d. uniform sequence (see False-Positive Bound above). All three functions raise `ValueError` when the identity is outside 0 to 1, the rate is outside the open interval (0, 1), or an adapter contains characters that are not IUPAC DNA letters.

`adapter_false_positive_bound(adapter_length, max_errors)` returns the bound for one adapter length and edit budget (`4 ** -adapter_length` when `max_errors` is 0) and raises `ValueError` for a length below 1 or a negative budget. `adapter_error_budget(adapter_length, min_identity, max_false_positive_rate)` returns the edit budget the searches use, or -1 when the adapter is shorter than `minimum_adapter_length`; it raises `ValueError` for a length below 1 and for the identity and rate values the searches reject. `minimum_adapter_length(max_false_positive_rate)` returns the shortest adapter whose exact match stays within the rate and raises `ValueError` for a rate outside (0, 1).

The built-in adapters, with the edits each tolerates at the default rate (`min_identity` 0.70 and 0.75 give the same budgets):

| Adapter | Length (bp) | Edits tolerated |
|---------|-------------|-----------------|
| `ONT_LSK109_top` | 28 | 3 |
| `ONT_LSK109_bottom` | 22 | 2 |
| `ONT_RNA_adapter` | 24 | 2 |
| `ONT_rapid_top` | 50 | 6 |
| `ONT_rapid_barcode` | 36 | 4 |
| `PacBio_SMRTbell` | 45 | 5 |
| `Illumina_Nextera_read1` | 33 | 3 |

`PacBio_barcoded` is a deprecated lookup alias for `Illumina_Nextera_read1`; the sequence is an Illumina adapter, not a PacBio adapter.

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

# Detect adapters with the built-in table: any offset, both orientations
for match in filtering.detect_adapters(read_records[0].sequence):
    print(f"{match.adapter_name}: {match.location} at {match.position}-{match.end_position}, "
          f"{match.edit_distance} edits, false-positive bound {match.false_positive_bound:.1e}")

# Trim adapters from the read ends; a lower rate is stricter (useful for GC- or AT-rich genomes)
trimmed = filtering.trim_adapters(read_records)
strict = filtering.trim_adapters(read_records, max_false_positive_rate=1e-12)

# Search only the ONT adapters
ont_trimmed = filtering.trim_adapters(read_records, filtering.ONT_ADAPTERS)

# Split chimeric reads at internal adapters
split_reads = filtering.split_chimeric_reads(read_records)
print(f"Split: {len(read_records)} -> {len(split_reads)} reads")

# Edit budget and false-positive bound of one adapter
adapter_length = len(filtering.ONT_ADAPTERS["ONT_LSK109_bottom"])
budget = filtering.adapter_error_budget(adapter_length)
bound = filtering.adapter_false_positive_bound(adapter_length, budget)
print(f"{adapter_length} bp tolerates {budget} edits, bound {bound:.1e}")
print(f"Shortest adapter at the default rate: {filtering.minimum_adapter_length()} bp")
```

## Configuration

- **Environment prefix**: `LR_` is honoured only by `workflow.pipelines.load_pipeline_config` and `utils.batch.process_batch` (`LR_BATCH_WORKERS`); the functions on this page read no environment variables
- **Optional dependencies**: none in use (the modules try to import numpy but never use it)
- Built-in adapter databases include ONT (LSK109, rapid, RNA) and PacBio (SMRTbell, barcoded) adapters; the shortest is 22 bp
- Default quality filter threshold: Q7 (ONT standard); default minimum read length: 1000 bp
- Adapter matching tolerates substitutions, insertions and deletions. `min_identity` caps the edits (0.70 for `detect_adapters`, 0.75 for `trim_adapters` and `split_chimeric_reads`), and `max_false_positive_rate` (1e-9 per read offset) and the seed limit `len(adapter) // 7 - 1` can lower the cap
- Minimum adapter length at the default rate: 15 bp; shorter custom adapters are matched exactly, with a warning

## Related Modules

- `longread.io` -- Read FAST5/BAM files for quality analysis
- `longread.assembly` -- Assembly from the filtered reads
- `longread.utils` -- Batch processing and QC summary generation
- `longread.visualization` -- `plot_read_length_histogram`, `plot_quality_vs_length`
