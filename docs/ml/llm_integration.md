# LLM Integration (Ollama)

Integration with local large language models via Ollama for text generation, multi-turn chat and prompt-driven analysis of biological data.

## Overview

The LLM module is a small client for a local Ollama server (plain `urllib`, no extra dependencies), enabling:
- Text generation and multi-turn chat, with optional streaming
- Prompt templates and builders for bioinformatics tasks
- Composable chains (prompt, sequential, map-reduce, router, conversation)
- Listing the models installed on the server

There are no dedicated helpers for biological queries, code assistance or documentation generation: these are ordinary prompts sent with `generate` or `chat`, as shown below. Chain examples are in the [Ollama package README](../../src/metainformant/ml/llm/ollama/README.md).

## Quick Start

```python
from metainformant.ml.llm.ollama.client import OllamaClient
from metainformant.ml.llm.ollama.config import OllamaConfig

# Initialize client (the model must already be pulled with `ollama pull`)
client = OllamaClient(OllamaConfig(model="llama3:latest"))

# Ask about biological data by putting the data in the prompt
expression_summary = "gene,log2_fold_change\nMYC,2.4\nTP53,-0.3\nCDK4,1.9"
response = client.generate(
    "Which genes are upregulated in the cancer sample?\n\n" + expression_summary
)

print(response.text)
print(f"{response.eval_count} tokens in {response.duration_seconds:.2f}s")
```

`generate` returns a `GenerateResponse` (`text`, `model`, `done`, `total_duration`, `prompt_eval_count`, `eval_count`, `context`); `chat` returns a `ChatResponse` whose `.text` is the assistant message. Extra keyword arguments to either call (for example `temperature=0.2`) are passed to Ollama as sampling options.

## Configuration

### Basic Setup

```python-snippet
from metainformant.ml.llm.ollama.client import OllamaClient
from metainformant.ml.llm.ollama.config import OllamaConfig

# Configure Ollama
config = OllamaConfig(
    host="http://localhost:11434",
    model="llama2",
    temperature=0.7,
    num_predict=2048,  # maximum generated tokens
)

client = OllamaClient(config)
```

`OllamaConfig` fields and defaults: `host="http://localhost:11434"`, `model="smollm2:135m-instruct-q4_K_S"`, `temperature=0.7`, `top_p=0.9`, `top_k=40`, `timeout=120.0` (seconds), `max_retries=3` (the maximum number of attempts per request, the first one included, so `1` means "do not retry"; it must be at least 1), `stream=False`, `context_length=4096` (sent as `num_ctx`), `num_predict=None` (model default), `system_prompt=None` and `stop_sequences=[]`. Values are validated on construction (for example `temperature` must be between 0 and 2 and `max_retries` at least 1, otherwise `ValueError` is raised) and a `host` without a scheme gets `http://` prepended.

### Environment Variables

`OllamaConfig` and `OllamaClient` do not read any environment variables, including `OLLAMA_*` ones: the configuration comes only from the `OllamaConfig` you pass in. To drive the settings from the environment, read the variables in your own script (the names below are just a convention):

```python
import os
from metainformant.ml.llm.ollama.client import OllamaClient
from metainformant.ml.llm.ollama.config import OllamaConfig

config = OllamaConfig(
    host=os.environ.get("OLLAMA_HOST", "http://localhost:11434"),
    model=os.environ.get("OLLAMA_MODEL", "llama3:latest"),
    timeout=float(os.environ.get("OLLAMA_TIMEOUT", "300")),
)
client = OllamaClient(config)
```

## Features

### Biological Data Querying

```python
import pandas as pd
from metainformant.ml.llm.ollama.prompts import build_bioinformatics_prompt

# Put the data you want to ask about into the prompt as text
expression = pd.DataFrame({
    "gene": ["MYC", "TP53", "CDK4"],
    "log2_fold_change": [2.4, -0.3, 1.9],
    "adjusted_p_value": [0.001, 0.4, 0.02],
})

prompt = build_bioinformatics_prompt(
    task="Find genes with fold change > 2 and adjusted p-value < 0.05",
    data=expression.to_string(index=False),
    organism="Homo sapiens",
    output_format="structured",
)
response = client.generate(prompt)
print(response.text)
```

`build_bioinformatics_prompt(task, data, organism="Homo sapiens", output_format="structured")` returns a plain prompt string; the model, not the library, does the analysis, so check the answer against the data.

### Code Assistance

```python
from metainformant.ml.llm.ollama.prompts import SystemPrompt

system = (
    SystemPrompt()
    .role("Python programmer who writes bioinformatics code")
    .instruction("Return only runnable code with short comments")
    .build()
)

# Generate analysis code
response = client.generate(
    "Create a Manhattan plot for GWAS results with matplotlib.",
    system=system,
    temperature=0.2,
)
print(response.text)

# Explain existing code
code = "def gc_content(seq):\n    return (seq.count('G') + seq.count('C')) / len(seq)"
response = client.generate(f"Explain what this function does:\n\n{code}", system=system)
print(response.text)
```

### Documentation Generation

```python
from metainformant.ml.llm.ollama.prompts import PromptTemplate

report_template = PromptTemplate(
    "Write a short {format} report for the analysis '{analysis_name}'.\n\nResults:\n{results}",
    {"format": "markdown"},
)

# Generate an analysis report from a results table
report = client.generate(
    report_template.format(
        analysis_name="differential_expression",
        results=expression.to_string(index=False),
    )
)
print(report.text)
```

## Available Models

Ollama supports various models. Install models with `ollama pull <model>`, and list what is installed with `client.list_models()`:

| Model | Description | Use Case |
|-------|-------------|----------|
| `smollm2:135m-instruct-q4_K_S` | Tiny (library default) | Quick tests |
| `llama2` | General purpose | Code generation |
| `codellama` | Code-focused | Programming assistance |
| `mistral` | Fast inference | Quick queries |
| `llama3` | Latest Llama | General tasks |

**Note**: Biomedical-specific models require custom fine-tuning or third-party models. Check [Ollama library](https://ollama.ai/library) for available models.

## Prompts and Templates

The module includes pre-built prompts for common biological tasks in `metainformant.ml.llm.ollama.prompts`:

```python
from metainformant.ml.llm.ollama.prompts import (
    GENE_ANALYSIS_TEMPLATE,
    SEQUENCE_ANALYSIS_TEMPLATE,
    LITERATURE_SUMMARY_TEMPLATE,
    PromptTemplate,
    build_conversation_messages,
)

# Use pre-built prompts (organism defaults to "Homo sapiens")
prompt = GENE_ANALYSIS_TEMPLATE.format(gene_name="TP53")

# SEQUENCE_ANALYSIS_TEMPLATE needs `sequence` (and `sequence_type`, default "DNA");
# LITERATURE_SUMMARY_TEMPLATE needs `title` and `abstract` (and `style`, default "concise")
prompt = SEQUENCE_ANALYSIS_TEMPLATE.format(sequence="ATGGCGTACGTTAGC")

# Define your own template, for example for GWAS interpretation
gwas_template = PromptTemplate(
    "Interpret these GWAS association results for {trait}.\n\n"
    "Top variants:\n{hits}\n\n"
    "Summarize likely biological relevance and next validation steps.",
    {"trait": "the phenotype"},
)
prompt = gwas_template.format(hits="variant_3 chr1:1327 p=3.0e-09\nvariant_20 chr2:5110 p=8.0e-08")

# Build a chat message list for client.chat()
messages = build_conversation_messages(
    "You are a bioinformatics assistant.",
    ["What does a GC-rich promoter suggest?"],
)
response = client.chat(messages)
print(response.text)
```

There is no GWAS-specific prompt in the module; templates are plain `{placeholder}` strings.

## Integration with METAINFORMANT

### With GWAS Module

```python
from metainformant.gwas import run_gwas
from metainformant.ml.llm.ollama.client import OllamaClient

# Run GWAS analysis (gwas_config is a dictionary of GWAS settings)
results = run_gwas(vcf_path, phenotype_path, gwas_config)

# Pick the strongest association rows (dicts with "variant_id", "chrom", "pos", "p_value", ...)
rows = results["results"]["association_results"]
top_hits = sorted(rows, key=lambda row: row.get("p_value", 1.0))[:10]
hits = "\n".join(
    f"{row['variant_id']} chr={row.get('chrom', '?')} pos={row.get('pos', '?')} p={row.get('p_value', 1.0):.2e}"
    for row in top_hits
)

# Interpret results with the LLM through the template defined above
client = OllamaClient()
interpretation = client.generate(gwas_template.format(trait="varroa resistance", hits=hits))
print(interpretation.text)
```

### With RNA Module

```python
from metainformant.rna.analysis.expression_analysis import differential_expression
from metainformant.ml.llm.ollama.client import OllamaClient

# Run differential expression (counts: genes x samples, conditions: one label per sample)
de_results = differential_expression(counts, conditions)

# Generate biological interpretation of the strongest genes
client = OllamaClient()
top_genes = de_results.nsmallest(10, "adjusted_p_value")
biological_story = client.generate(
    "Explain the likely biology behind these differentially expressed genes:\n\n"
    + top_genes.to_string(index=False)
)
print(biological_story.text)
```

## Error Handling

```python
from metainformant.ml.llm.ollama.client import OllamaClient

try:
    client = OllamaClient()
    response = client.generate("Analyze this data")
except ConnectionError:
    print("Ollama server not running. Start with: ollama serve")
except RuntimeError as e:
    # Errors returned by the server, for example when the model has not been pulled
    print(f"Ollama API error: {e}")
```

`ConnectionError` is raised when the server cannot be reached (`generate` and `chat` make up to `max_retries` attempts in total, with pauses of 1 s, 2 s, 4 s, ... between them) and `RuntimeError` carries the HTTP status and message of an API error, which is not retried. A request needs at least one attempt: `OllamaConfig` rejects `max_retries` below 1, and a client whose `config.max_retries` was lowered below 1 after construction raises `ValueError` before it contacts the server. `client.is_available()` returns `True` when the server answers `list_models()`.

## Troubleshooting

### Ollama Not Running

```bash
# Start Ollama server
ollama serve

# Or start in background
nohup ollama serve > ollama.log 2>&1 &
```

### Model Not Found

```bash
# List available models
ollama list

# Pull a model
ollama pull llama2
```

### Connection Timeout

```python
# Increase timeout
from metainformant.ml.llm.ollama.client import OllamaClient
from metainformant.ml.llm.ollama.config import OllamaConfig

config = OllamaConfig(timeout=600)  # 10 minutes
client = OllamaClient(config)
```

## Performance Considerations

- **Local models**: Use GPU for best performance
- **Streaming**: `client.generate(prompt, stream=True, stream_callback=print)` hands you each chunk as it arrives
- **Token limits**: Set `num_predict` to cap generated tokens
- **Context window**: `context_length` (default 4096) sets the context size; raise it for large prompts
- **Retries and timeouts**: `timeout` (default 120 s) applies to every request; `max_retries` (default 3, at least 1) is the total number of attempts and applies to connection failures in `generate` and `chat`

## Related Documentation

- [ML Index](./index.md) - Machine learning overview
- [Deep Learning](./deep_learning.md) - Neural networks
- [Ollama Documentation](https://github.com/ollama/ollama) - Ollama project
