# AGENTS.md - LLM Module

## Module Purpose

Provides local LLM inference capabilities for bioinformatics workflows using Ollama and other local backends.

## Key Components

| Component | Path | Purpose |
|-----------|------|---------|
| `OllamaClient` | `ollama/client.py` | Core client for Ollama API |
| `OllamaConfig` | `ollama/config.py` | Configuration dataclass |
| `PromptTemplate` | `ollama/prompts.py` | Template with variable substitution |
| `Chain` | `ollama/chains.py` | Base class for composable workflows |

## Development Guidelines

### Adding New Backends

1. Create new subpackage (e.g., `llm/llamacpp/`)
2. Implement same interface as `OllamaClient`
3. Export from `llm/__init__.py`
4. Add tests in `tests/ml/llm/`

### Testing Strategy

- **Real functional tests only** - no mocks
- Use smallest available model (`smollm2:135m`)
- Test with short prompts to minimize latency
- Skip tests gracefully if Ollama unavailable

### Common Patterns

```python
# Check if Ollama is available before tests
client = OllamaClient()
if not client.is_available():
    pytest.skip("Ollama not available")

# Use streaming for long responses
response = client.generate(prompt, stream=True, stream_callback=print)

# Build chains for complex workflows
chain = SequentialChain([step1, step2, step3])
```

## Dependencies

- No external Python dependencies (uses urllib)
- Requires Ollama server running locally
- Uses standard library only for portability
