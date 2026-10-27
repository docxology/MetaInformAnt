# deep_learning

Sub-package of `metainformant.ml`: a pure-NumPy 1D convolutional network for DNA/RNA sequences (inference only; there is
no training loop). See the module's `AGENTS.md`/`README.md` for the domain overview and
[`docs/ml/deep_learning.md`](../../../../docs/ml/deep_learning.md) for the guide.

## Contents

| File | Purpose |
|------|---------|
| `sequences.py` | One-hot encoding, `conv1d_forward`, pooling, `SequenceCNNConfig`, `SequenceCNNWeights`, `init_sequence_cnn_weights`, `predict_sequences` |

## Key Functions

| Function | Description |
|----------|-------------|
| `one_hot_encode()` / `batch_encode()` | One-hot encode a sequence, or a padded batch (columns A, C, G, T or U, N) |
| `conv1d_forward()` | Stride-aware 1D convolution of a `(batch, length, in_channels)` array; raises `ValueError` for a sequence shorter than the kernel |
| `init_sequence_cnn_weights()` | He-initialized random weights drawn from a local `numpy.random.RandomState` (default seed 42) |
| `predict_sequences()` | Scores for a list of sequences; raises `ValueError` when the padded or truncated sequences are too short for the convolution layers |

## Short sequences

With stride 1 every convolution layer shortens a sequence by `kernel_size - 1`.

- `conv1d_forward` raises a `ValueError` naming the sequence length and the kernel size when the input is shorter than the
  kernel.
- `predict_sequences` checks the whole stack before running any layer: the (padded or truncated) sequences need at least
  `1 + sum(kernel_size - 1)` positions (11 for the default kernel sizes 8 and 4), otherwise it raises a `ValueError`
  that names the length, the kernel sizes and the minimum length. Use longer sequences, a larger `max_length`, or
  smaller kernels.

## Usage

```python
from metainformant.ml.deep_learning.sequences import SequenceCNNConfig, init_sequence_cnn_weights, predict_sequences

config = SequenceCNNConfig(n_filters=[16, 8], kernel_sizes=[4, 3], n_classes=2)
weights = init_sequence_cnn_weights(config, random_state=0)
scores = predict_sequences(["ATGCGTACGT", "GCTATTAGC"], weights, config)  # shape (2, 2)
```
