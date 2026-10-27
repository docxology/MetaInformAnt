# Deep Learning for Biological Sequences

Deep learning module for biological sequence analysis: one-hot encoding of DNA/RNA sequences and a pure-NumPy 1D convolutional neural network for sequence classification or regression.

## Overview

The deep learning module (`metainformant.ml.deep_learning.sequences`) uses only NumPy (no PyTorch or TensorFlow) and provides:
- One-hot encoding of single sequences and of padded batches
- A 1D convolution forward pass with ReLU and global max/average pooling
- A configurable sequence CNN: `SequenceCNNConfig`, He-initialized `SequenceCNNWeights` and `predict_sequences` for inference

The module has no training loop, no pre-trained checkpoints, no transformer or recurrent architectures and no GPU support. `init_sequence_cnn_weights` returns random weights; to use trained weights, fill a `SequenceCNNWeights` yourself.

## Quick Start

```python
from metainformant.ml.deep_learning.sequences import (
    SequenceCNNConfig,
    batch_encode,
    init_sequence_cnn_weights,
    predict_sequences,
)

# One-hot encode DNA sequences (columns: A, C, G, T, N), padded to the longest sequence
encoded = batch_encode(["ATGCGTACGT", "GCTATTAGC", "ATGCGTACGTAA"])
print(encoded.shape)  # (3, 12, 5)

# Configure a small CNN and initialize random weights
config = SequenceCNNConfig(n_filters=[16, 8], kernel_sizes=[4, 3], n_classes=2)
weights = init_sequence_cnn_weights(config, random_state=0)

# Run inference: a (n_sequences, n_classes) array of scores
scores = predict_sequences(["ATGCGTACGT", "GCTATTAGC", "ATGCGTACGTAA"], weights, config)
print(scores.shape)  # (3, 2)
```

## Available Model Architectures

The module implements one architecture family, a 1D convolutional network, together with its building blocks.

### Supported Architectures

| Architecture | Description | Use Case |
|--------------|-------------|----------|
| CNN | Stacked 1D convolutions with ReLU, global max or average pooling and a dense output layer | Motif-like pattern detection, sequence classification or regression |

**Note**: The weights are not trained by this module. Pre-trained DNA transformers, ESM2 weights and recurrent models are not part of it.

## API Reference

### Sequence Encoding

```python-snippet
def one_hot_encode(
    sequence: str,
    alphabet: str = "DNA",
) -> np.ndarray: ...

def batch_encode(
    sequences: list[str],
    max_length: int | None = None,
    alphabet: str = "DNA",
) -> np.ndarray: ...
```

`one_hot_encode` returns a `(sequence_length, 5)` float32 array whose columns are A, C, G, T (U for `alphabet="RNA"`) and N; any other character is encoded as N and lowercase is accepted. `batch_encode` returns a `(n_sequences, max_length, 5)` array, padding shorter sequences with zero rows and truncating longer ones; `max_length` defaults to the longest sequence.

### Convolution Building Blocks

```python
import numpy as np
from metainformant.ml.deep_learning.sequences import (
    Conv1DParams,
    batch_encode,
    conv1d_forward,
    global_max_pool,
    relu,
)

x = batch_encode(["ATGCGTACGT", "GCTATTAGC"])        # (2, 10, 5)
rng = np.random.default_rng(0)
layer = Conv1DParams(
    weights=rng.normal(size=(8, 5, 4)).astype(np.float32),  # (out_channels, in_channels, kernel_size)
    bias=np.zeros(8, dtype=np.float32),
)

hidden = relu(conv1d_forward(x, layer, stride=1))      # (2, 7, 8)
pooled = global_max_pool(hidden)                       # (2, 8)
print(hidden.shape, pooled.shape)
```

`conv1d_forward(x, params, stride=1)` takes a `(batch, length, in_channels)` array and returns `(batch, out_length, out_channels)`; it uses explicit Python loops, so it suits small models and inference rather than large batches. A sequence shorter than the kernel makes the output length non-positive and raises `ValueError`, so pad with `max_length`. `global_avg_pool` is the average-pooling counterpart of `global_max_pool`.

### Sequence CNN

```python
from metainformant.ml.deep_learning.sequences import (
    SequenceCNNConfig,
    init_sequence_cnn_weights,
    predict_sequences,
)

config = SequenceCNNConfig(
    n_filters=[32, 64],     # filters per convolution layer
    kernel_sizes=[8, 4],    # kernel size per layer
    n_classes=2,            # 1 for regression
    alphabet="DNA",         # or "RNA"
    pool="max",             # or "avg"
)
weights = init_sequence_cnn_weights(config, random_state=42)

scores = predict_sequences(
    ["ACGT" * 10, "TTGA" * 10],
    weights,
    config,
    max_length=40,          # pad or truncate to this length
)
print(scores.shape)  # (2, 2)
```

`SequenceCNNWeights` holds `conv_layers` (a list of `Conv1DParams`), `dense_weights` of shape `(n_filters[-1], n_classes)` and `dense_bias`; `predict_sequences` returns the raw scores (logits) of the dense layer.

## Configuration

The network is configured through the `SequenceCNNConfig` dataclass shown above (defaults: `n_filters=[32, 64]`, `kernel_sizes=[8, 4]`, `n_classes=2`, `alphabet="DNA"`, `pool="max"`). There is no YAML model or training configuration in this module.

## Integration with METAINFORMANT

### With DNA Module

```python-snippet
from metainformant.dna.sequence import core as dna_seqs
from metainformant.ml.deep_learning.sequences import (
    SequenceCNNConfig,
    init_sequence_cnn_weights,
    predict_sequences,
)

# Read DNA sequences (a dict of sequence ID -> sequence string)
dna_data = dna_seqs.read_fasta("genome.fasta")

# Score every sequence with a CNN (weights are random here; load trained ones for real use)
config = SequenceCNNConfig()
weights = init_sequence_cnn_weights(config)
scores = predict_sequences(list(dna_data.values()), weights, config, max_length=200)
```

### With ML Classification

```python
import numpy as np
from metainformant.ml.deep_learning.sequences import (
    SequenceCNNConfig,
    batch_encode,
    conv1d_forward,
    global_max_pool,
    init_sequence_cnn_weights,
    relu,
)
from metainformant.ml.evaluation.validation import train_test_split_biological
from metainformant.ml.models.classification import BiologicalClassifier

# Toy sequences: 40 GC-rich (label 1) and 40 AT-rich (label 0)
rng = np.random.default_rng(0)

def random_sequence(gc_fraction, length=60):
    probabilities = [(1 - gc_fraction) / 2, gc_fraction / 2, gc_fraction / 2, (1 - gc_fraction) / 2]
    return "".join(rng.choice(list("ACGT"), size=length, p=probabilities))

sequences = [random_sequence(0.7) for _ in range(40)] + [random_sequence(0.3) for _ in range(40)]
labels = np.array([1] * 40 + [0] * 40)

# Convolutional features: encode, convolve, pool over the sequence
config = SequenceCNNConfig(n_filters=[16], kernel_sizes=[6], n_classes=2)
weights = init_sequence_cnn_weights(config, random_state=0)
hidden = batch_encode(sequences, alphabet=config.alphabet)
for layer in weights.conv_layers:
    hidden = relu(conv1d_forward(hidden, layer))
features = global_max_pool(hidden)  # (n_sequences, n_filters)

# Use the features with a classical classifier
X_train, X_test, y_train, y_test = train_test_split_biological(features, labels, test_size=0.25, random_state=42)
model = BiologicalClassifier(algorithm="logistic_regression", random_state=42)
model.fit(X_train, y_train)
print(f"Test accuracy: {(model.predict(X_test) == y_test).mean():.2f}")
```

## Performance Notes

- **Pure NumPy**: the forward pass runs on the CPU with Python loops; there is no GPU, mixed-precision or multi-device support
- **Batch size**: memory grows with `n_sequences x max_length x channels`; encode and score in chunks for large inputs
- **Padding**: with stride 1 every encoded sequence must be at least `sum(kernel_sizes) - len(kernel_sizes) + 1` long, so set `max_length` at least that large

## Related Documentation

- [ML Index](./index.md) - Machine learning overview
- [Features](./features.md) - Feature engineering
- [Models](./models.md) - Classical ML models
- [Evaluation](./evaluation.md) - Model validation
