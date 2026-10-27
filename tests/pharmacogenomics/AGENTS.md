# AGENTS.md — `MetaInformAnt/tests/pharmacogenomics/`

tests for the pharmacogenomics module.
Files (verified 2026-10-27): __init__.py, test_pharmacogenomics_alleles.py, test_pharmacogenomics_annotations.py, test_pharmacogenomics_clinical.py, test_pharmacogenomics_metabolism.py, test_pharmacogenomics_visualization.py, test_drug_interactions_depth.py, test_reporting_depth.py, test_cpic_phenotype_definitions.py, test_allele_definition_loading.py, test_package_determinism.py.

- `test_cpic_phenotype_definitions.py` — phenotype bands and allele-function rules checked against the CPIC
  diplotype tables, plus the single source of truth shared by `alleles` and `metabolism`.
- `test_allele_definition_loading.py` — the `load_allele_definitions` source/filepath contract and the
  genes (UGT1A1, CYP1A2) that have activity scores but no calling definitions.
- `test_package_determinism.py` — randomness audit: no random-number use, results independent of `PYTHONHASHSEED`.
- `test_reporting_depth.py` — `generate_clinical_report`: `drugs=None` analyses every CPIC drug of the genotyped genes,
  `drugs=[]` none, recommendations follow the classified phenotype (a diplotype alone suffices, a phenotype in any
  form is accepted), plus the exports and `add_disclaimer`.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/pharmacogenomics`
Repo-wide policy: see the repository-root `AGENTS.md`.