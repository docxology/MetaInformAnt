# Specification: ollama

## 🎯 Scope
Ollama backend for local LLM inference.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 5 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `chains.py`
- `client.py`
- `config.py`
- `prompts.py`

### Behavior Contracts

#### `config.py`
- `OllamaConfig.max_retries` (default 3) is the maximum number of attempts per request, the first one included, so 1
  means no retry. `__post_init__` raises `ValueError` for a value below 1 (like `timeout <= 0`, `top_k < 1` and a
  `temperature` outside 0.0-2.0).

#### `client.py`
- `OllamaClient._request` raises `ValueError` before any network call when `config.max_retries` was lowered below 1 after
  the config was created. A connection failure (`URLError`) is retried with pauses of `2 ** attempt` seconds until the
  attempts are used up, then raises `ConnectionError`; an `HTTPError` raises `RuntimeError` at once.
