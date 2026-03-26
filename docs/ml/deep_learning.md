# Deep Learning for Biological Sequences

Deep learning module for biological sequence analysis, including sequence embeddings, transformers, and neural network architectures tailored for genomics and transcriptomics data.

## Overview

The deep learning module provides neural network architectures specifically designed for biological sequence analysis, with support for:
- Sequence embedding generation
- Transformer-based models for DNA/RNA/Protein
- Convolutional networks for motif detection
- Recurrent networks for sequence modeling

## Quick Start

```python
from metainformant.ml.deep_learning import sequences

# Generate embeddings for DNA sequences
embeddings = sequences.embed_dna_sequences(
    sequences=["ATGCGT...", "GCTAT..."],
    model_name="dna-transformer",
    embedding_dim=128
)

# Train a classifier on embeddings
from metainformant.ml import classification
model = classification.train_classifier(embeddings, labels, method="mlp")
```

## Available Model Architectures

The module supports various neural network architectures for biological sequence analysis. Model availability depends on installed checkpoints and dependencies.

### Supported Architectures

| Architecture | Description | Use Case |
|--------------|-------------|----------|
| Transformer | Attention-based models | Sequence classification, embedding
| CNN | Convolutional networks | Motif detection
| RNN/LSTM | Recurrent networks | Sequence modeling
| K-mer Embedder | K-mer feature extraction | Feature representation

**Note**: Specific model checkpoints (e.g., pre-trained DNA transformers, ESM2 weights) must be installed separately. Check the source code in `src/metainformant/ml/deep_learning/` for currently available models.

## API Reference

### Sequence Embedding

```python
from metainformant.ml.deep_learning import sequences

# Generate embeddings
embeddings = sequences.embed_sequences(
    sequences: list[str],
    model_name: str,
    batch_size: int = 32,
    device: str = "cuda",
) -> np.ndarray
```

### Fine-tuning

```python
# Fine-tune a model on your biological data
trained_model = sequences.fine_tune(
    base_model="dna-transformer",
    train_sequences=train_seqs,
    train_labels=train_labels,
    val_sequences=val_seqs,
    epochs=10,
    learning_rate=1e-4,
)
```

## Configuration

Deep learning models can be configured via YAML:

```yaml
models:
  dna_transformer:
    architecture: transformer
    num_layers: 6
    num_heads: 8
    embedding_dim: 256
    dropout: 0.1

training:
  batch_size: 32
  epochs: 50
  optimizer: adam
  learning_rate: 0.001
  scheduler: cosine
```

## Integration with METAINFORMANT

### With DNA Module

```python
from metainformant.dna import sequences as dna_seqs
from metainformant.ml.deep_learning import sequences as dl_seqs

# Read DNA sequences
dna_data = dna_seqs.read_fasta("genome.fasta")

# Generate embeddings
embeddings = dl_seqs.embed_sequences(
    list(dna_data.values()),
    model_name="dna-transformer"
)
```

### With ML Classification

```python
from metainformant.ml import classification

# Use deep learning embeddings with ML classifiers
model = classification.train_classifier(
    embeddings,
    labels,
    method="ensemble",
    hyperparameter_tuning=True
)
```

## Performance Notes

- **GPU**: Use CUDA for optimal performance
- **Batch Size**: Adjust based on available memory
- **Mixed Precision**: Enable for faster training on modern GPUs
- **Distributed**: Support for multi-GPU training

## Related Documentation

- [ML Index](./index.md) - Machine learning overview
- [Features](./features.md) - Feature engineering
- [Models](./models.md) - Classical ML models
- [Evaluation](./evaluation.md) - Model validation