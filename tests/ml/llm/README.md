# llm tests

Ollama LLM client/chain tests (`test_ollama_chains.py`, `test_ollama_client.py`) — require a local Ollama endpoint; skip otherwise (unverified skip mechanism).

## Files

- `__init__.py`
- `test_ollama_chains.py`
- `test_ollama_client.py`

Run from the repo root: `pytest tests/ml/llm/ -v`. Tests follow the real-implementation policy (no mocks).
