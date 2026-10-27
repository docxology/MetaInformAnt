# AGENTS.md — `MetaInformAnt/tests/metabolomics/`

tests for the metabolomics module.
Files (verified 2026-10-27): __init__.py, test_metabolomics_differential_abundance.py, test_metabolomics_identification.py, test_metabolomics_identification_validation.py, test_metabolomics_io.py, test_metabolomics_io_format_claims.py, test_metabolomics_io_strict_parsing.py, test_metabolomics_missing_value_handling.py, test_metabolomics_pathways.py, test_metabolomics_pathways_pvalue.py, test_metabolomics_visualization.py, test_metabolomics_visualization_numpy_compat.py.

- `test_metabolomics_differential_abundance.py` — Welch t-test statistics against an independent scipy reference.
- `test_metabolomics_identification_validation.py` — `ppm_tolerance` must be positive and finite; `ion_mode` is
  case-insensitive and anything but positive/negative raises.
- `test_metabolomics_io_format_claims.py` — pins the readers and writers that exist (CSV and MGF; no mzML/mzXML) and
  that the module docstring says the same.
- `test_metabolomics_io_strict_parsing.py` — `read_csv` returns NaN for missing cells, raises on other non-numeric
  cells and ragged rows, and fills `metadata`; `read_mgf` reads `.5 99`/`+200.0 20` peak lines, skips comments and
  raises on a bad header or peak line.
- `test_metabolomics_missing_value_handling.py` — TIC/median normalization refuse NaN, log2 keeps it, `fold_change`
  returns NaN for a metabolite with a missing value, and the read, impute, normalize order is correct.
- `test_metabolomics_pathways_pvalue.py` — the hypergeometric p-value is exact against scipy oracles, and an invalid
  `background_size` raises instead of returning 0.0.
- `test_metabolomics_visualization_numpy_compat.py` — peak areas do not depend on a NumPy 2-only function.

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/metabolomics` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.
