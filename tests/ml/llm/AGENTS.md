# AGENTS.md — `MetaInformAnt/tests/ml/llm`

Verified against disk 2026-08-30 (doc-realization fleet pass). Ollama LLM client/chain tests (`test_ollama_chains.py`, `test_ollama_client.py`) — require a local Ollama endpoint; skip otherwise (unverified skip mechanism).

## Layout

- `__init__.py`
- `test_ollama_chains.py`
- `test_ollama_client.py`

## Invariants & gotchas

- Real-implementation policy: tests exercise real I/O and real computation — no mock frameworks (see `docs/REAL_IMPLEMENTATION_POLICY.md`).
- Run: `pytest tests/ml/llm/ -v` (single dir) from the MetaInformAnt repo root; fast/full modes via `bash scripts/package/test.sh`.
- Do not add new top-level test dirs without updating `tests/AGENTS.md`.
