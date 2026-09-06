# AGENTS.md — `MetaInformAnt/tests/ml/llm/`

Nested test subpackage under `tests/ml/` (verified 2026-09-05).
Contents: __init__.py, test_ollama_chains.py, test_ollama_client.py, test_ollama_local_server.py, test_ollama_prompts.py.
`test_ollama_client.py` and `test_ollama_chains.py` require a live Ollama
server (opt in via `METAINFORMANT_RUN_OLLAMA_INTEGRATION=1`); the other
files are fully offline.
Follow the repo real-implementation test policy — no mocks.
Repo-wide policy: see the repository-root `AGENTS.md`.
