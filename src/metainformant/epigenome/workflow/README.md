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
| `run_methylation_workflow()` | Methylation pipeline: load BEDgraph and Bismark `.cov` files, per-sample statistics, text report; returns the per-site `sites` list |
| `run_chipseq_workflow()` | ChIP-seq pipeline: detect the peak layout (narrowPeak, broadPeak or BED), load peaks, filter on q-value (an ordinary FDR limit) and length, per-sample statistics, filtered narrowPeak output; returns the kept `peaks` |
| `run_atacseq_workflow()` | ATAC-seq pipeline: the same steps for ATAC-seq peaks; returns the kept `peaks` |
| `integrate_epigenome_results()` | Cross-assay integration: methylation-ChIP and methylation-ATAC site/peak proximity, and ChIP-ATAC peak overlap with the base-pair Jaccard index |

## Workflow Steps

1. Load configuration via `load_epigenome_config()`
2. Run individual assay workflows (methylation, ChIP-seq, ATAC-seq); each returns the `sites` or `peaks` list that the integration step reads
3. Integrate results across assays with `integrate_epigenome_results()`
4. Generate integration report with association statistics

## Usage

```python
from metainformant.epigenome.workflow.workflow import (
    integrate_epigenome_results,
    load_epigenome_config,
    run_atacseq_workflow,
    run_chipseq_workflow,
    run_methylation_workflow,
)

config = load_epigenome_config("config/epigenome.yaml")
meth_results = run_methylation_workflow("data/methylation/", "output/methylation/", config)
chip_results = run_chipseq_workflow("data/chipseq/", "output/chipseq/", config)
atac_results = run_atacseq_workflow("data/atacseq/", "output/atacseq/", config)
integrated = integrate_epigenome_results(meth_results, chip_results, atac_results, output_dir="output/integrated/")
```
