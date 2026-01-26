# Ollama Backend

Local LLM inference using the Ollama server.

## Overview

This subpackage provides a complete interface to Ollama for local large language model inference. It includes:

- **`OllamaClient`** - Core client for API communication
- **`OllamaConfig`** - Configuration with validation
- **`PromptTemplate`** - Composable prompt templates
- **Chains** - Sequential, MapReduce, Router patterns

## Installation

Ollama must be installed and running:

```bash
# Install Ollama
brew install ollama

# Start the service
brew services start ollama

# Pull a model
ollama pull smollm2:135m-instruct-q4_K_S
```

## Usage Examples

### Basic Generation

```python
from metainformant.ml.llm.ollama import OllamaClient, OllamaConfig

# Default config uses smollm2 for fast responses
client = OllamaClient()

# Simple generation
response = client.generate("What is DNA?")
print(response.text)
print(f"Generated {response.eval_count} tokens in {response.duration_seconds:.2f}s")
```

### Custom Configuration

```python
config = OllamaConfig(
    model="llama3:latest",
    temperature=0.3,      # Lower = more deterministic
    top_p=0.9,
    timeout=120.0,
)
client = OllamaClient(config)
```

### Chat Conversations

```python
from metainformant.ml.llm.ollama import ChatMessage

messages = [
    ChatMessage("system", "You are a molecular biology expert."),
    ChatMessage("user", "Explain transcription."),
]

response = client.chat(messages)
print(response.text)

# Continue the conversation
messages.append(response.message)
messages.append(ChatMessage("user", "What enzymes are involved?"))
response = client.chat(messages)
```

### Streaming Responses

```python
def on_token(text: str):
    print(text, end="", flush=True)

response = client.generate(
    "Write a haiku about genetics.",
    stream=True,
    stream_callback=on_token,
)
print()  # Newline after streaming
```

### Prompt Templates

```python
from metainformant.ml.llm.ollama import PromptTemplate

template = PromptTemplate(
    "Analyze the {organism} gene {gene_name}. "
    "Focus on: {focus_area}",
    {"organism": "human", "focus_area": "function"}
)

prompt = template.format(gene_name="BRCA1")
response = client.generate(prompt)
```

### Chain Workflows

```python
from metainformant.ml.llm.ollama import (
    SequentialChain,
    PromptChain,
    PromptTemplate,
)

# Multi-step analysis pipeline
extract = PromptChain(PromptTemplate("Extract key points from: {item}"))
analyze = PromptChain(PromptTemplate("Analyze these points: {item}"))
summarize = PromptChain(PromptTemplate("Summarize: {item}"))

pipeline = SequentialChain([extract, analyze, summarize])
result = pipeline.run({"item": document_text})

if result.success:
    print(result.value)
else:
    print(f"Error: {result.error}")
```

### Bioinformatics Helpers

```python
from metainformant.ml.llm.ollama import (
    build_bioinformatics_prompt,
    GENE_ANALYSIS_TEMPLATE,
    SystemPrompt,
)

# Pre-built gene analysis template
prompt = GENE_ANALYSIS_TEMPLATE.format(
    gene_name="TP53",
    organism="Homo sapiens",
)

# Build custom system prompt
system = (SystemPrompt()
    .role("computational biologist")
    .context("analyzing RNA-seq differential expression")
    .instruction("cite specific genes and pathways")
    .build())
```

## Available Models

Check models on your system:

```python
client = OllamaClient()
for model in client.list_models():
    print(f"{model.name}: {model.size_gb:.1f}GB")
```

Common models:

| Model | Size | Best For |
|-------|------|----------|
| `smollm2:135m-instruct-q4_K_S` | 102MB | Quick tests, simple tasks |
| `gemma3:4b` | 3.3GB | Balanced performance |
| `llama3:latest` | 4.7GB | Complex reasoning |
| `codellama:7b` | 3.8GB | Code generation |

## Error Handling

```python
try:
    response = client.generate(prompt)
except ConnectionError as e:
    print(f"Ollama server not available: {e}")
except RuntimeError as e:
    print(f"API error: {e}")
```

## Performance Tips

1. **Use smallest model that works** - `smollm2` for simple tasks
2. **Enable streaming** for long responses to get early output
3. **Reuse client instances** - they maintain connection state
4. **Set appropriate timeouts** - larger models need more time
5. **Use low temperature** (0.1-0.3) for deterministic outputs
