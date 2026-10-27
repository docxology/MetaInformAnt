# Agent Directives: eqtl

## 🤖 Role
Specialized agent context for the `eqtl` component.

## 🛠️ Tools & Capabilities
- **Context**: Zero-mocks test suite for the `metainformant.eqtl` pipeline — parameter resolution, synthetic cohorts, variant calling, and variant statistics — and for the eQTL scan contracts of `metainformant.gwas.finemapping.eqtl` and `scripts/eqtl/run_eqtl_real.py`.
- **Pattern**: Test Suite Pattern
- **Files** (verified 2026-10-27): __init__.py, test_eqtl_failed_fits.py, test_eqtl_real_script_labels.py, test_eqtl_unique_variant_ids.py, test_pipeline.py, test_support.py, test_synthetic.py, test_variant_calling.py, test_variant_stats.py (8 test files).

## Contract tests

- `test_eqtl_failed_fits.py`: a fit that cannot be made (fewer than 3 complete observations, a constant genotype, a zero standard error) has NaN statistics, `status == "error"` and a reason in `cis_eqtl_scan`, `conditional_eqtl` and `eqtl_effect_sizes` (NaN `r_squared`); `eqtl_summary_stats` leaves it out of its Benjamini-Hochberg family and reports `n_tests_failed`; `write_eqtl_results` writes `NA`.
- `test_eqtl_unique_variant_ids.py`: synthetic variant IDs are `var_<full gene id>_<i>`, unique across transcripts that share a chromosome; `cis_eqtl_scan` raises `ValueError` for a repeated ID in the gene positions, the variant positions or the genotype index, and tests each pair once on its own genotypes.
- `test_eqtl_real_script_labels.py`: `scripts/eqtl/run_eqtl_real.py` runs as a separate process on a quantification tree the test writes itself (no network) and labels what it produced: `gene_positions_synthetic.tsv` and the `*_synthetic_genotypes.*` outputs, the provenance keys of the two JSON summaries, the log warning that the associations are not biological findings, and no output named `*_real`.

## ⚠️ Rules & Constraints
- **Imports**: Prefer absolute imports from `metainformant`.
- **I/O**: Use `metainformant.core.io` for all file operations.
- **Logging**: Use `metainformant.core.utils.logging`.
- **Run**: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/eqtl`
