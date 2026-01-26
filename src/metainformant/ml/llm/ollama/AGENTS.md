# AGENTS.md - Ollama Backend

## Component Purpose

Provides the core Ollama integration for local LLM inference.

## File Roles

| File | Purpose |
|------|---------|
| `config.py` | Configuration dataclass with validation |
| `client.py` | HTTP client for Ollama REST API |
| `prompts.py` | Template system for prompt construction |
| `chains.py` | Composable workflow patterns |

## Key Design Decisions

### Pure urllib (No Dependencies)

Uses only Python standard library to avoid dependency conflicts in bioinformatics environments.

### Dataclass-Based Responses

All API responses wrapped in dataclasses for type safety and IDE support.

### Composable Chains

Chain pattern enables building complex workflows from simple operations.

## Testing Notes

```python
# Always check availability first
client = OllamaClient()
if not client.is_available():
    pytest.skip("Ollama not running")

# Use smollm2 for fast tests
config = OllamaConfig(model="smollm2:135m-instruct-q4_K_S")
client = OllamaClient(config)

# Keep prompts short in tests
response = client.generate("Say hello.", temperature=0.1)
```

## Extension Points

1. **New chain types**: Extend `Chain` base class
2. **Custom prompts**: Create `PromptTemplate` instances
3. **System prompts**: Use `SystemPrompt` builder
4. **Response processing**: Add methods to response dataclasses
