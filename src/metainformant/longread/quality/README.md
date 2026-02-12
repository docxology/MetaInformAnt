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
| `metrics.read_length_stats()` | Comprehensive read length distribution statistics |
| `metrics.quality_score_distribution()` | Analyze quality score distribution across reads |
| `metrics.estimate_accuracy()` | Estimate read accuracy from Phred quality scores |
| `metrics.calculate_throughput()` | Calculate sequencing throughput over time |
| `filtering.filter_by_length()` | Filter reads by minimum/maximum length |
| `filtering.filter_by_quality()` | Filter reads by minimum quality score |
| `filtering.trim_adapters()` | Detect and trim ONT/PacBio adapter sequences |
| `filtering.detect_adapters()` | Detect adapter sequences without trimming |
| `filtering.split_chimeric_reads()` | Detect and split chimeric reads |

## Usage

```python
from metainformant.longread.quality import metrics, filtering

stats = metrics.read_length_stats(reads)
n50 = metrics.calculate_n50(read_lengths)
filtered = filtering.filter_by_length(reads, min_length=1000)
trimmed = filtering.trim_adapters(filtered)
```
