# simulation tests

pytest suite for the `simulation` domain of METAINFORMANT. Tests import from `src/metainformant/simulation` and follow the repo's real-implementation policy.

## Files

- `__init__.py`
- `test_simulation.py`
- `test_simulation_agents.py`
- `test_simulation_popgen.py`
- `test_simulation_rna_advanced.py`
- `test_simulation_workflow.py`

Run from the repo root: `pytest tests/simulation/ -v`. Tests follow the real-implementation policy (no mocks).
