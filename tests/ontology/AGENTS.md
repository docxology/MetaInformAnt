# AGENTS.md — `MetaInformAnt/tests/ontology/`

tests for the ontology module.
Files (verified 2026-08-29): __init__.py, test_ontology_api.py, test_ontology_comprehensive.py, test_ontology_enrichment.py, test_ontology_go_basic.py, test_ontology_obo_parser.py, test_ontology_query.py, test_ontology_serialization.py … (11 files).


## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/ontology` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.