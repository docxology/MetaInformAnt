# SPEC.md - LLM Module Technical Specification

## Architecture

```
ml/llm/
├── __init__.py          # Package exports
├── README.md            # User documentation
├── AGENTS.md            # Agent guidance
├── SPEC.md              # This file
└── ollama/              # Ollama backend
    ├── __init__.py      # Subpackage exports
    ├── config.py        # OllamaConfig dataclass
    ├── client.py        # OllamaClient class
    ├── prompts.py       # Prompt templates
    ├── chains.py        # Composable chains
    └── README.md        # Ollama-specific docs
```

## Core Interfaces

### OllamaClient

```python
class OllamaClient:
    def __init__(self, config: Optional[OllamaConfig] = None): ...
    def is_available(self) -> bool: ...
    def list_models(self) -> list[ModelInfo]: ...
    def generate(self, prompt: str, **kwargs) -> GenerateResponse: ...
    def chat(self, messages: list[ChatMessage], **kwargs) -> ChatResponse: ...
    def conversation(self, message: str, **kwargs) -> ChatResponse: ...
```

### Chain Interface

```python
class Chain(ABC, Generic[T, R]):
    def run(self, input_data: T) -> ChainResult[R]: ...
    def then(self, next_chain: Chain) -> SequentialChain: ...
```

## API Endpoints

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/api/generate` | POST | Text completion |
| `/api/chat` | POST | Chat completion |
| `/api/tags` | GET | List models |

## Response Types

### GenerateResponse

- `text`: Generated content
- `model`: Model used
- `eval_count`: Tokens generated
- `total_duration`: Time in nanoseconds

### ChatResponse

- `message`: ChatMessage with role/content
- `model`: Model used
- `eval_count`: Tokens generated

## Error Handling

- `ConnectionError`: Server unreachable
- `RuntimeError`: API errors (4xx/5xx)
- Automatic retry with exponential backoff
