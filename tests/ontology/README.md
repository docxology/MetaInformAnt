# ontology tests

pytest suite for the `ontology` domain of METAINFORMANT. Tests import from `src/metainformant/ontology` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_ontology_api.py`
- `test_ontology_comprehensive.py`
- `test_ontology_enrichment.py`
- `test_ontology_go_basic.py`
- `test_ontology_obo_parser.py`
- `test_ontology_query.py`
- `test_ontology_serialization.py`
- `test_ontology_serialize.py`
- `test_ontology_types.py`
- `test_ontology_visualization.py`

Run from the repo root: `pytest tests/ontology/ -v`. Tests follow the real-implementation policy (no mocks).
