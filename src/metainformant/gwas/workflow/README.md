# GWAS Workflow

End-to-end GWAS workflow orchestration: configuration management, VCF loading, QC filtering, association testing, multi-trait analysis, and results export.

## Contents

| File | Purpose |
|------|---------|
| `workflow.py` | Re-export module aggregating config and execution APIs |
| `workflow_config.py` | GWASWorkflowConfig dataclass, YAML config loading, phenotype loading |
| `workflow_execution.py` | Pipeline execution: single-trait GWAS, multi-trait GWAS |

## Key Classes and Functions

| Name | Description |
|------|-------------|
| `GWASWorkflowConfig` | Dataclass holding VCF path, phenotypes, QC thresholds, methods |
| `load_gwas_config()` | Load YAML config into GWASWorkflowConfig with env overrides |
| `validate_gwas_config()` | Validate config paths exist and parameters are in range |
| `execute_gwas_workflow()` | Full pipeline: load VCF, QC, LD pruning, structure, association, correction, heritability, summary tables and plots; returns `status`, `steps_completed` and `outputs` (with `n_tests_failed`) |
| `run_gwas()` | Single-trait association with automatic method selection; returns `results["results"]` with `n_tests_failed` |
| `run_multi_trait_gwas()` | Run association across multiple phenotype columns; each `trait_results[trait]` has `n_tests` and `n_tests_failed` |

## Workflow Steps

1. Load config via `load_gwas_config()`
2. Parse VCF and apply QC filters (QUAL when `min_qual` is set, MAF, missingness, HWE)
3. Compute population structure (PCA, kinship)
4. Run association test (linear, logistic, or mixed model)
5. Apply multiple testing correction
6. Export summary statistics and significant hits

## Usage

```python
from metainformant.gwas.workflow.workflow import load_gwas_config, execute_gwas_workflow

config = load_gwas_config("config/gwas.yaml")
results = execute_gwas_workflow(config)
```

## Failed Tests, MAF and Genotype Layout

- A test that cannot be fitted keeps its row with NaN statistics (see the [analysis README](../analysis/README.md#failed-tests)), so rows line up with the variants and the summary tables print it as `NA`. It is left out of the corrections, the plots, the top hit and the fine-mapping, and the workflows report how many there were as `n_tests_failed`.
- `execute_gwas_workflow()` and `run_gwas()` give every association result a `maf`, the frequency of the genotype row the test ran on, and all three workflows pass the genotypes to the writers, so the MAF column of `summary_statistics.tsv` is filled. `execute_gwas_workflow()` also draws `maf_spectrum_plot.png` from these values; the MAF panel of `amellifera_pipeline.generate_visualizations()` uses the same `maf` and is skipped, with a warning, when no result carries one.
- The genotypes of `parse_vcf_full()` and `apply_qc_filters()` are sample-major. `execute_gwas_workflow()` declares that layout (`layout="sample_major"`) when it transposes them to variants x samples; without a declared layout, `_extract_genotype_matrix()` infers it from the shape and assumes a square matrix is variant-major, logging a warning.
- `execute_gwas_workflow()` passes the association results that have a p-value to `generate_all_plots()`, which takes records, a DataFrame or a `.json`, `.tsv`, `.txt` or `.csv` file.
