# Utils

Utility functions for long-read sequencing analysis, providing batch processing of multiple samples with parallelism and run summary generation for reporting.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports batch and summary submodules |
| `batch.py` | Concurrent batch processing with progress tracking and error handling |
| `summary.py` | Run summary generation, aggregation, and export for analysis reports |

## Key Functions

| Function | Description |
|----------|-------------|
| `batch.process_batch()` | Process multiple items concurrently with progress tracking |
| `batch.batch_filter_reads()` | Batch-apply read filtering across multiple files |
| `batch.batch_compute_metrics()` | Batch-compute QC metrics for multiple samples |
| `summary.generate_qc_summary()` | Generate QC summary from analysis results |
| `summary.generate_assembly_summary()` | Generate assembly quality summary |
| `summary.generate_methylation_summary()` | Generate methylation analysis summary |
| `summary.generate_sv_summary()` | Generate structural variant summary |
| `summary.build_run_summary()` | Build comprehensive RunSummary from pipeline results |
| `summary.export_run_summary()` | Export summary to JSON/text file |
| `summary.compare_run_summaries()` | Compare multiple run summaries side-by-side |

## Usage

```python
from metainformant.longread.utils import batch, summary

result = batch.process_batch(items, process_func, max_workers=4)
run = summary.build_run_summary(pipeline_results, sample_name="sample_1")
summary.export_run_summary(run, "output/run_summary.json")
```
