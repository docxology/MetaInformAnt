# Quality Reporting

QC report generation, threshold checking, multi-sample aggregation, and trend analysis for sequencing metrics. Provides configurable pass/warn/fail thresholds and batch-level quality summaries.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `multiqc_integration` module |
| `multiqc_integration.py` | QC thresholds, reporting, aggregation, and trend analysis |

## Key Functions

| Function | Description |
|----------|-------------|
| `default_qc_thresholds()` | Return default thresholds for common sequencing metrics |
| `check_qc_thresholds()` | Check sample metrics against pass/warn/fail thresholds |
| `aggregate_sample_qc()` | Aggregate QC metrics across multiple samples |
| `generate_qc_report()` | Generate a structured QC report from collected metrics |
| `qc_trend_analysis()` | Analyze QC metric trends over time or across batches |

## Usage

```python
from metainformant.quality.reporting import multiqc_integration

thresholds = multiqc_integration.default_qc_thresholds()
status = multiqc_integration.check_qc_thresholds(metrics, thresholds)
report = multiqc_integration.generate_qc_report(all_metrics)
```
