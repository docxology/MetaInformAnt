# Information Workflow

Batch information-theoretic analysis workflows: entropy profiling across sequence collections, cross-dataset comparison, and report generation.

## Contents

| File | Purpose |
|------|---------|
| `workflows.py` | Batch entropy analysis, dataset comparison, markdown/text report generation |

## Key Functions

| Function | Description |
|----------|-------------|
| `batch_entropy_analysis()` | Run entropy analysis across a collection of sequences |
| `information_workflow()` | Complete information-theoretic analysis pipeline for sequence data |
| `compare_datasets()` | Compare entropy distributions and sequence diversity between two datasets |
| `information_report()` | Generate formatted report (markdown or text) from analysis results |

## Workflow Steps

1. Provide sequence data to `information_workflow()` or `batch_entropy_analysis()`
2. Each sequence is profiled for entropy, complexity, and information content
3. Results are aggregated with summary statistics
4. Optionally compare two datasets via `compare_datasets()`
5. Generate publication-ready reports with `information_report()`

## Usage

```python
from metainformant.information.workflow.workflows import (
    information_workflow, compare_datasets, information_report
)

results = information_workflow(sequences, method="shannon")
comparison = compare_datasets(dataset_a, dataset_b)
report = information_report(results, format="markdown")
```
