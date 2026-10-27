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
| `calling.call_methylation_from_signal()` | Call methylation from per-read log-likelihood ratios: the `threshold` model (logistic) or the `bayesian` model, which pools the reads at each site under a Beta prior (`prior_alpha`, `prior_beta`; needs numpy) |
| `calling.aggregate_methylation()` | Aggregate per-read calls to per-site methylation frequencies |
| `calling.detect_dmrs()` | Detect differentially methylated regions via segmentation |
| `calling.methylation_pattern_analysis()` | Analyze single-read epiallele methylation patterns |
| `calling.compute_methylation_stats()` | Compute summary statistics for methylation data |

## Usage

```python
from metainformant.longread.methylation import calling

calls = calling.call_methylation_from_signal(signal_data, model="bayesian", prior_alpha=1.0, prior_beta=1.0)
methylation = calling.aggregate_methylation(calls, min_coverage=5)
dmrs = calling.detect_dmrs(methylation_a, methylation_b)  # both are aggregate_methylation results
stats = calling.compute_methylation_stats(methylation["sites"])
```

See [docs/longread/methylation.md](../../../../docs/longread/methylation.md) for the Bayesian model.
