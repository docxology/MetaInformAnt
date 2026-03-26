# LLM Integration (Ollama)

Integration with local large language models via Ollama for biological data analysis, natural language queries, and automated documentation generation.

## Overview

The LLM module provides integration with Ollama, enabling:
- Natural language queries about biological data
- Automated documentation generation
- Code assistance and explanation
- Data analysis guidance

## Quick Start

```python
from metainformant.ml.llm import OllamaClient

# Initialize client
client = OllamaClient(model="llama2")

# Query about biological data
response = client.query(
    "What genes are upregulated in the cancer sample?",
    context=expression_data
)

print(response)
```

## Configuration

### Basic Setup

```python
from metainformant.ml.llm import OllamaClient
from metainformant.ml.llm.config import LLMConfig

# Configure Ollama
config = LLMConfig(
    base_url="http://localhost:11434",
    model="llama2",
    temperature=0.7,
    max_tokens=2048,
)

client = OllamaClient(config)
```

### Environment Variables

```bash
# Set Ollama configuration
export OLLAMA_BASE_URL=http://localhost:11434
export OLLAMA_MODEL=llama2
export OLLAMA_TIMEOUT=300
```

## Features

### Biological Data Querying

```python
# Query gene expression data
response = client.query_biological(
    "Find genes with fold change > 2 and p-value < 0.05",
    data=expression_matrix,
    annotations=gene_annotations
)

# Query variant data
response = client.query_variants(
    "What are the missense variants in the BRCA1 gene?",
    vcf_data=vcf_data
)
```

### Code Generation

```python
# Generate analysis code
code = client.generate_code(
    "Create a Manhattan plot for GWAS results",
    library="metainformant"
)

# Explain existing code
explanation = client.explain_code(
    "Explain this GWAS association test function"
)
```

### Documentation Generation

```python
# Generate module documentation
doc = client.generate_docs(
    module="metainformant.gwas.association",
    format="markdown"
)

# Generate analysis report
report = client.generate_report(
    analysis_name="differential_expression",
    results=deseq_results
)
```

## Available Models

Ollama supports various models. Install models with `ollama pull <model>`:

| Model | Description | Use Case |
|-------|-------------|----------|
| `llama2` | General purpose | Code generation |
| `codellama` | Code-focused | Programming assistance |
| `mistral` | Fast inference | Quick queries |
| `llama3` | Latest Llama | General tasks |

**Note**: Biomedical-specific models require custom fine-tuning or third-party models. Check [Ollama library](https://ollama.ai/library) for available models.

## Prompts and Templates

The module includes pre-built prompts for common biological tasks:

```python
from metainformant.ml.llm import prompts

# Use pre-built prompts
prompt = prompts.gwas_interpretation(
    pvalues=results['pvalue'],
    effect_sizes=results['beta'],
    annotations=gene_db
)
```

## Integration with METAINFORMANT

### With GWAS Module

```python
from metainformant.gwas import run_gwas
from metainformant.ml.llm import OllamaClient

# Run GWAS analysis
results = run_gwas(vcf_path, phenotype_path, config)

# Interpret results with LLM
client = OllamaClient()
interpretation = client.interpret_gwas(results)
```

### With RNA Module

```python
from metainformant.rna import differential_expression
from metainformant.ml.llm import OllamaClient

# Run differential expression
de_results = differential_expression(counts, conditions)

# Generate biological interpretation
client = OllamaClient()
biological_story = client.explain_de_results(de_results)
```

## Error Handling

```python
from metainformant.ml.llm import OllamaClient

try:
    client = OllamaClient()
    response = client.query("Analyze this data")
except ConnectionError:
    print("Ollama server not running. Start with: ollama serve")
except FileNotFoundError:
    print("Model not available. Pull with: ollama pull llama2")
except Exception as e:
    print(f"Error: {e}")
```

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
config = LLMConfig(timeout=600)  # 10 minutes
client = OllamaClient(config)
```

## Performance Considerations

- **Local models**: Use GPU for best performance
- **Batch queries**: Group queries for efficiency
- **Caching**: Enable response caching for repeated queries
- **Context window**: Use appropriate context size for your data

## Related Documentation

- [ML Index](./index.md) - Machine learning overview
- [Deep Learning](./deep_learning.md) - Neural networks
- [Ollama Documentation](https://github.com/ollama/ollama) - Ollama project