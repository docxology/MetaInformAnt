# AGENTS.md — `MetaInformAnt/tests/pharmacogenomics/`

tests for the pharmacogenomics module.
Files (verified 2026-10-28): __init__.py, test_pharmacogenomics_alleles.py, test_pharmacogenomics_annotations.py, test_pharmacogenomics_clinical.py, test_pharmacogenomics_metabolism.py, test_pharmacogenomics_visualization.py, test_drug_interactions_depth.py, test_reporting_depth.py, test_cpic_phenotype_definitions.py, test_allele_definition_loading.py, test_package_determinism.py, test_unknown_genes_alleles_diplotypes.py.

- `test_cpic_phenotype_definitions.py` — phenotype bands and allele-function rules checked against the CPIC
  diplotype tables, plus the single source of truth shared by `alleles` and `metabolism`.
- `test_unknown_genes_alleles_diplotypes.py` — what CPIC has no answer for is never given a normal one: a gene
  without a CPIC phenotype (CYP1A2 has activity values only) raises `ValueError` naming the supported genes in every
  phenotype function; an allele with no activity value makes the diplotype Indeterminate (activity score `None`, the
  allele named) in `alleles`, `metabolism`, population frequencies, `dose_adjustment` and the report; an unscored
  `Diplotype` (`activity_score` `None`) is scored from its alleles, never read as a poor metabolizer.
- `test_allele_definition_loading.py` — the `load_allele_definitions` source/filepath contract, the
  genes (UGT1A1, CYP1A2) that have activity scores but no calling definitions, and that a definition file needs an
  `activity_value` for every allele and a gene-keyed JSON file must contain the requested gene.
- `test_package_determinism.py` — randomness audit: no random-number use, results independent of `PYTHONHASHSEED`.
- `test_reporting_depth.py` — `generate_clinical_report`: `drugs=None` analyses every CPIC drug of the genotyped genes,
  `drugs=[]` none, recommendations follow the classified phenotype (a diplotype alone suffices, a phenotype in any
  form is accepted), Indeterminate genes (unknown allele, gene without a CPIC phenotype) and their exports, plus the
  exports and `add_disclaimer`.
- `test_pharmacogenomics_clinical.py` also covers phenotype labels in any form in `drug_interaction`
  (`"poor"` is a poor metabolizer), ACMG frequency criteria that need a frequency, contradictory ACMG evidence
  (VUS), and `check_gnomad_frequency` inputs; `test_pharmacogenomics_annotations.py` covers drug-label
  classification from text, PGx term matching, `list_actionable_genes` levels and strict CPIC allele-definition parsing.


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV .venv/bin/python -m pytest -q tests/pharmacogenomics`
Repo-wide policy: see the repository-root `AGENTS.md`.