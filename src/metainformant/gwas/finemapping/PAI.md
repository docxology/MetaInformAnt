# Personal AI Infrastructure (PAI) - finemapping

## 🧠 Context & Intent
- **Path**: `src/metainformant/gwas/finemapping`
- **Purpose**: GWAS statistical fine-mapping submodule: Wakefield approximate Bayes factors, credible sets, SuSiE regression, colocalization, conditional analysis and eQTL scans.
- **Domain**: gwas

## 🏗️ Virtual Hierarchy
- **Type**: Source Code
- **Parent**: `gwas`

## 📝 Maintenance Notes
- **System**: Part of the METAINFORMANT Domain layer.
- **Style**: Strict type hinting, real implementations in tests.
- **Stability**: API boundaries should be respected.
- **Numerics**: Bayes factors, PIPs and colocalization posteriors are accumulated in log space. The scale of each effect estimate is an input (`standard_errors=` or `sample_size=` with `maf=` and `case_fraction=`); without one the z-score scale is used and `variance_model` reports `"z_scale"` (per trait under `bayes_factor_model` for colocalization).

## 🔄 AI Workflows
- **Modification**: Run functional tests in `tests/` before committing; the fine-mapping suites are `tests/gwas/test_gwas_finemapping_abf.py` and `tests/gwas/test_gwas_finemapping_statistics.py`.
- **Documentation**: Update `SPEC.md` if architectural patterns change, and keep the `SPEC.md` signatures and the `README.md` function tables in step with the code.
