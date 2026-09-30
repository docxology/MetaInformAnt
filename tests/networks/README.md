# networks tests

pytest suite for the `networks` domain of METAINFORMANT. Tests import from `src/metainformant/networks` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_networks_community.py`
- `test_networks_comprehensive.py`
- `test_networks_graph.py`
- `test_networks_pathway.py`
- `test_networks_ppi.py`
- `test_networks_regulatory.py`

Run from the repo root: `pytest tests/networks/ -v`. Tests follow the real-implementation policy (no mocks).
