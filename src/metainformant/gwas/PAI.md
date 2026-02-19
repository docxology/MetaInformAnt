# Personal AI Infrastructure (PAI) - gwas

## 🧠 Context & Intent

- **Path**: `src/metainformant/gwas`
- **Purpose**: Genome-Wide Association Studies (GWAS) module for METAINFORMANT.
- **Domain**: metainformant

## 🏗️ Virtual Hierarchy

- **Type**: Source Code
- **Parent**: `metainformant`
- **Sub-packages**: `analysis`, `data`, `finemapping`, `heritability`, `visualization`, `workflow`

## 📝 Maintenance Notes

- **System**: Part of the METAINFORMANT Domain layer.
- **Style**: Strict type hinting, no mocks in tests, pure Python fallbacks for numpy.
- **Stability**: API boundaries should be respected.
- **Benchmarking**: `analysis/benchmarking.py` provides pilot-run → full-genome time extrapolation using known computational complexity models (O(n·m), O(n²·m), etc.).
- **Expression**: `data/expression.py` provides `ExpressionLoader` for integrating Amalgkit kallisto quantification with GWAS/eQTL pipelines.

## 🔄 AI Workflows

- **Modification**: Run functional tests in `tests/test_gwas_*.py` before committing.
- **Documentation**: Update `SPEC.md` if architectural patterns change.
- **Benchmarking**: Use `benchmark_subset_run()` on small data before running full-genome analysis to estimate compute costs.
- **Expression Integration**: Use `ExpressionLoader` to bridge RNA-seq TPM data into eQTL analysis.
