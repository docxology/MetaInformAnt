# AGENTS.md — `MetaInformAnt/tests/ontology/`

tests for the ontology module.
Files (verified 2026-10-26): __init__.py, test_ontology_annotation_regression.py, test_ontology_api.py, test_ontology_comprehensive.py, test_ontology_enrichment.py, test_ontology_go_basic.py, test_ontology_obo_definitions.py, test_ontology_obo_parser.py, test_ontology_query.py, test_ontology_query_exports.py, test_ontology_query_semantics.py, test_ontology_serialize.py, test_ontology_serialization.py, test_ontology_types.py, test_ontology_visualization.py (15 files).

## Coverage notes
- `test_ontology_obo_definitions.py`: `def:` parsing (quoted text, dbxref list, escapes, trailing modifiers, malformed values), `metadata["definition_xrefs"]`, `[Typedef]` definitions, and OBO save/load round trips of definitions and dbxrefs.
- `test_ontology_query_exports.py`: every public function of `query.py` and `serialize.py` is re-exported by the `metainformant.ontology.query` package and listed in `__all__`; named imports, star imports and the submodules stay importable.
- `test_ontology_query_semantics.py`: `ancestors`, `descendants` and `common_ancestors` exclude the queried term itself, for `is_a` and other relation types; the docstring examples run and give the documented results.

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/ontology` (unverified — not run in this pass).
Repo-wide policy: see the repository-root `AGENTS.md`.