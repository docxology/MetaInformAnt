# Methylation

Signal-level methylation calling from nanopore long-read data, including per-site aggregation, differentially methylated region detection, and single-read epiallele pattern analysis.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports calling submodule |
| `calling.py` | Methylation calling pipeline: signal-level, aggregation, DMR, patterns |

## Key Functions

| Function | Description |
|----------|-------------|
| `calling.call_methylation_from_signal()` | Call 5mC/6mA from signal data using threshold or likelihood models |
| `calling.aggregate_methylation()` | Aggregate per-read calls to per-site methylation frequencies |
| `calling.detect_dmrs()` | Detect differentially methylated regions via segmentation |
| `calling.methylation_pattern_analysis()` | Analyze single-read epiallele methylation patterns |
| `calling.compute_methylation_stats()` | Compute summary statistics for methylation data |

## Usage

```python
from metainformant.longread.methylation import calling

calls = calling.call_methylation_from_signal(signal_data, model="threshold")
sites = calling.aggregate_methylation(calls, min_coverage=5)
dmrs = calling.detect_dmrs(group_a_sites, group_b_sites)
stats = calling.compute_methylation_stats(sites)
```
