# AGENTS.md — `MetaInformAnt/tests/structural_variants/`

tests for the structural-variants module.
Files (verified 2026-10-27): __init__.py, test_structural_variants.py, test_structural_variants_detection_advanced.py, test_structural_variants_population.py, test_quality_filter_depth.py, test_functional_impact_depth.py, test_merge_depth.py, test_detection_depth.py, test_numpy_fallback_parity.py, test_sv_classification_direction.py, test_sv_consensus_median.py, test_sv_calling_fixes.py, test_breakpoint_refinement_fixes.py, test_interval_index_correctness.py, test_dosage_gain_reachable.py, test_quality_filter_fixes.py, test_survivor_merge_fixes.py, test_cnv_fixes.py, test_population_variance_explained.py, test_reproducibility.py.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/structural_variants -p no:cacheprovider` (verified 2026-10-27: 481 passed).
- Defect regressions (`test_*_fixes.py`, `test_dosage_gain_reachable.py`, `test_interval_index_correctness.py`,
  `test_population_variance_explained.py`, `test_reproducibility.py`) check hand-derived values, not the
  implementation's own output. Cost claims are ratios between two sizes timed in the same process, never absolute
  limits tuned to one machine. `test_cnv_fixes.py` keeps a scalar CBS as an oracle for the vectorised one;
  `test_reproducibility.py` runs fresh interpreters with different `PYTHONHASHSEED` values against the package copy
  the test process imported.
Repo-wide policy: see the repository-root `AGENTS.md`.
