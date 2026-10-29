# Quality

Long-read quality assessment and filtering, providing N50/Nx statistics, read length distributions, quality score analysis, adapter trimming, and chimeric read detection.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports metrics and filtering submodules |
| `metrics.py` | N50, Nx, read length distributions, accuracy estimation, throughput |
| `filtering.py` | Length/quality filtering, adapter trimming, chimeric read splitting |

## Key Functions

| Function | Description |
|----------|-------------|
| `metrics.calculate_n50()` | Calculate N50 read length statistic |
| `metrics.calculate_nx()` | Calculate arbitrary Nx percentile statistic |
| `metrics.read_length_stats()` | Comprehensive read length distribution statistics; takes lengths (integers, numpy integers, whole-number floats such as `1000.0`; a fractional or NaN length raises `ValueError`), dicts or objects |
| `metrics.quality_score_distribution()` | Analyze quality score distribution across reads |
| `metrics.estimate_accuracy()` | Estimate read accuracy from Phred quality scores |
| `metrics.calculate_throughput()` | Calculate sequencing throughput over time |
| `filtering.filter_by_length()` | Filter reads by minimum/maximum length |
| `filtering.filter_by_quality()` | Filter reads by minimum quality score |
| `filtering.trim_adapters()` | Detect and trim ONT/PacBio adapter sequences |
| `filtering.detect_adapters()` | Detect adapter sequences without trimming |
| `filtering.split_chimeric_reads()` | Detect and split chimeric reads |
| `filtering.adapter_false_positive_bound()` | Bound on the chance that a random read offset ends a spurious adapter placement, for an adapter length and edit budget |
| `filtering.adapter_error_budget()` | Edits tolerated for an adapter length, `min_identity` and false-positive rate (-1 when even an exact match exceeds the rate; the searches then match exactly, with a warning) |
| `filtering.minimum_adapter_length()` | Shortest adapter whose exact match stays within a false-positive rate (15 bp at the default) |

The default adapter table is `ALL_KNOWN_ADAPTERS` (`ONT_ADAPTERS`, `PACBIO_ADAPTERS`, `ILLUMINA_ADAPTERS`); the Illumina Nextera read-1 sequence is `Illumina_Nextera_read1` (it used to be listed as `PacBio_barcoded`, which still resolves with a `DeprecationWarning`). Where adapter placements overlap, the lower edit rate per adapter base wins and, at an equal rate, the longer adapter, so a longer adapter cannot take a nested shorter adapter's exact site by spending edits.

The adapter functions locate adapters by semi-global alignment, so substitutions, insertions and deletions count as edits, and they find an adapter at any offset in a read. `trim_adapters()`, `detect_adapters()` and `split_chimeric_reads()` take `max_false_positive_rate` (default 1e-9 per read offset); matches report `edit_distance` and `false_positive_bound`. An adapter shorter than `minimum_adapter_length()` is searched by exact match only, with a warning. See [docs/longread/quality.md](../../../../docs/longread/quality.md) for the matching model and the bound.

## Usage

```python
from metainformant.longread.quality import metrics, filtering

stats = metrics.read_length_stats(reads)
n50 = metrics.calculate_n50(read_lengths)
filtered = filtering.filter_by_length(reads, min_length=1000)
trimmed = filtering.trim_adapters(filtered)  # ReadRecord objects
matches = filtering.detect_adapters(trimmed[0].sequence)  # AdapterMatch objects for one read sequence
```
