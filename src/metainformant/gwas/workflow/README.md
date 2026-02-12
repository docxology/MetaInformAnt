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
| `execute_gwas_workflow()` | Full pipeline: load VCF, QC, structure, association, output |
| `run_gwas()` | Single-trait association with automatic method selection |
| `run_multi_trait_gwas()` | Run association across multiple phenotype columns |

## Workflow Steps

1. Load config via `load_gwas_config()`
2. Parse VCF and apply QC filters (MAF, missingness, HWE)
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
