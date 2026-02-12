# Epigenome Workflow

End-to-end workflow orchestration for epigenomic analyses, coordinating methylation, ChIP-seq, and ATAC-seq pipelines with cross-assay integration.

## Contents

| File | Purpose |
|------|---------|
| `workflow.py` | EpigenomeConfig, per-assay workflows, multi-assay integration |

## Key Functions

| Function | Description |
|----------|-------------|
| `load_epigenome_config()` | Load and validate YAML workflow configuration |
| `run_methylation_workflow()` | Full methylation pipeline: load, QC, DMR detection, reporting |
| `run_chipseq_workflow()` | Full ChIP-seq pipeline: peak loading, filtering, enrichment, motifs |
| `run_atacseq_workflow()` | Full ATAC-seq pipeline: peaks, TSS enrichment, accessibility |
| `integrate_epigenome_results()` | Cross-assay integration: methylation-ChIP, methylation-ATAC associations |

## Workflow Steps

1. Load configuration via `load_epigenome_config()`
2. Run individual assay workflows (methylation, ChIP-seq, ATAC-seq)
3. Integrate results across assays to find correlated epigenetic signals
4. Generate integration report with association statistics

## Usage

```python
from metainformant.epigenome.workflow.workflow import (
    load_epigenome_config, run_methylation_workflow, integrate_epigenome_results
)

config = load_epigenome_config("config/epigenome.yaml")
meth_results = run_methylation_workflow(config)
chip_results = run_chipseq_workflow(config)
integrated = integrate_epigenome_results(config)
```
