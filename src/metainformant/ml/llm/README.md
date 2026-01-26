# LLM Module

Local Large Language Model integration for METAINFORMANT bioinformatics workflows.

## Overview

This module provides composable interfaces to local LLM backends, with Ollama as the primary supported provider. It enables natural language processing capabilities for:

- **Gene/protein annotation summarization**
- **Literature analysis and extraction**
- **Data interpretation and reporting**
- **Interactive analysis assistance**

## Subpackages

| Package | Description |
|---------|-------------|
| `ollama` | Ollama backend for local LLM inference |

## Quick Start

```python
from metainformant.ml.llm import OllamaClient, OllamaConfig

# Initialize with default config (uses smollm2 for fast responses)
client = OllamaClient()

# Simple generation
response = client.generate("Explain the central dogma of molecular biology.")
print(response.text)

# Chat conversation
from metainformant.ml.llm import ChatMessage
messages = [
    ChatMessage("system", "You are a bioinformatics expert."),
    ChatMessage("user", "What is a FASTA file?"),
]
response = client.chat(messages)
print(response.text)
```

## Available Models

The module works with any Ollama-compatible model. Recommended models:

| Model | Size | Use Case |
|-------|------|----------|
| `smollm2:135m-instruct-q4_K_S` | 102MB | Fast, lightweight tasks |
| `gemma3:4b` | 3.3GB | Balanced quality/speed |
| `llama3:latest` | 4.7GB | High quality responses |

## Chain Composition

Build complex workflows with composable chains:

```python
from metainformant.ml.llm import SequentialChain, PromptTemplate
from metainformant.ml.llm.ollama import PromptChain

# Create a multi-step analysis chain
extract_chain = PromptChain(PromptTemplate("Extract key findings: {text}"))
summarize_chain = PromptChain(PromptTemplate("Summarize: {text}"))

pipeline = SequentialChain([extract_chain, summarize_chain])
result = pipeline.run({"text": document})
```

## Requirements

- Ollama server running locally (`brew services start ollama`)
- At least one model pulled (`ollama pull smollm2:135m-instruct-q4_K_S`)

## See Also

- [Ollama Documentation](https://ollama.ai/docs)
- [ollama/README.md](ollama/README.md) - Detailed Ollama backend docs
