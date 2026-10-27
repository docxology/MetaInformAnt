# Specification: workflow

## 🎯 Scope
GWAS workflow adapters for shared RNA-seq inputs.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 5 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `rna_quantification.py`
- `workflow.py`
- `workflow_config.py`
- `workflow_execution.py`

## Behavioral Contracts
- `execute_gwas_workflow(config, *, check=False)` returns a dictionary with `status` (`"completed"` or `"failed"`; `"validated"` or `"invalid"` with `check=True`), `success`, `steps_completed`, `errors` and `outputs`. `outputs` includes `association_results`, `n_tests_failed`, the correction results, `heritability`, `kinship_metadata` and `plots`.
- `run_gwas(...)` puts `n_tests_failed` in `results["results"]`; `run_multi_trait_gwas(...)` puts `n_tests` and `n_tests_failed` in each `results["trait_results"][trait]`.
- Association rows of a failed test keep their place with NaN statistics; corrections, plots, top hit and fine-mapping take only the rows with a p-value.
- `_extract_genotype_matrix(vcf_data, *, layout=None)` returns variants x samples. `layout` is `"sample_major"`, `"variant_major"` or `None`; with a layout the shape must fit it (`ValueError` otherwise), without one a square matrix is assumed variant-major and a warning is logged. `execute_gwas_workflow` passes `layout="sample_major"`.
