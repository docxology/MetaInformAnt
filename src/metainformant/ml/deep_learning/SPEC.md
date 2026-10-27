# Specification: deep_learning

## 🎯 Scope
Deep learning sub-package for biological sequence and image analysis.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 2 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `sequences.py`

### Behavior Contracts

#### `sequences.py`
- `conv1d_forward(x, params, stride=1)` raises `ValueError` when `x.shape[1]` is smaller than the kernel size of
  `params`; the message names the sequence length and the kernel size.
- `predict_sequences(sequences, weights, config, max_length=None)` encodes the batch, then requires
  `x.shape[1] >= 1 + sum(kernel_size - 1)` over `weights.conv_layers` (stride 1) and raises `ValueError` otherwise; the
  message names the length, the kernel sizes and the minimum. The check happens before any layer runs.
- `init_sequence_cnn_weights(config, random_state=42)` draws the weights from `numpy.random.RandomState(random_state)`;
  the global NumPy random state is not touched.
