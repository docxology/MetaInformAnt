"""Deep learning models for biological sequence analysis.

Implements 1D convolutional neural networks for DNA/RNA sequence
classification and regression tasks. Uses pure NumPy for inference
to avoid hard PyTorch/TensorFlow dependencies; training utilities
are available when torch is importable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


# ──────────────────────────────────────────────────────────────
# Encoding
# ──────────────────────────────────────────────────────────────

_DNA_MAP = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
_RNA_MAP = {"A": 0, "C": 1, "G": 2, "U": 3, "N": 4}


def one_hot_encode(
    sequence: str,
    alphabet: str = "DNA",
) -> np.ndarray:
    """One-hot encode a nucleotide sequence.

    Args:
        sequence: DNA or RNA string.
        alphabet: 'DNA' (ACGT) or 'RNA' (ACGU).

    Returns:
        2D array (seq_len × 4/5) of one-hot encoded values.
    """
    mapping = _DNA_MAP if alphabet.upper() == "DNA" else _RNA_MAP
    n_classes = len(mapping)
    seq_upper = sequence.upper()
    encoded = np.zeros((len(seq_upper), n_classes), dtype=np.float32)
    for i, base in enumerate(seq_upper):
        idx = mapping.get(base, mapping.get("N", n_classes - 1))
        encoded[i, idx] = 1.0
    return encoded


def batch_encode(
    sequences: list[str],
    max_length: int | None = None,
    alphabet: str = "DNA",
) -> np.ndarray:
    """One-hot encode a batch of sequences with padding.

    Args:
        sequences: List of nucleotide strings.
        max_length: Pad/truncate to this length. Defaults to max in batch.
        alphabet: 'DNA' or 'RNA'.

    Returns:
        3D array (batch × seq_len × channels).
    """
    if max_length is None:
        max_length = max(len(s) for s in sequences)

    mapping = _DNA_MAP if alphabet.upper() == "DNA" else _RNA_MAP
    n_classes = len(mapping)
    batch = np.zeros((len(sequences), max_length, n_classes), dtype=np.float32)

    for b, seq in enumerate(sequences):
        enc = one_hot_encode(seq[:max_length], alphabet)
        batch[b, : len(enc)] = enc

    return batch


# ──────────────────────────────────────────────────────────────
# Conv1D Layer (NumPy inference)
# ──────────────────────────────────────────────────────────────


@dataclass
class Conv1DParams:
    """Parameters for a 1D convolutional layer.

    Attributes:
        weights: 3D array (out_channels × in_channels × kernel_size).
        bias: 1D array (out_channels,).
    """

    weights: np.ndarray
    bias: np.ndarray


def conv1d_forward(
    x: np.ndarray,
    params: Conv1DParams,
    stride: int = 1,
) -> np.ndarray:
    """Forward pass for 1D convolution.

    Args:
        x: 3D input (batch × seq_len × in_channels).
        params: Conv1D weights and bias.
        stride: Convolution stride.

    Returns:
        3D output (batch × out_len × out_channels).
    """
    batch_size, seq_len, in_ch = x.shape
    out_ch, _, kernel_size = params.weights.shape
    out_len = (seq_len - kernel_size) // stride + 1

    output = np.zeros((batch_size, out_len, out_ch), dtype=np.float32)
    for b in range(batch_size):
        for i in range(out_len):
            start = i * stride
            patch = x[b, start : start + kernel_size, :]  # (kernel × in_ch)
            for c in range(out_ch):
                output[b, i, c] = np.sum(patch * params.weights[c].T) + params.bias[c]

    return output


def relu(x: np.ndarray) -> np.ndarray:
    """ReLU activation."""
    return np.maximum(x, 0)


def global_max_pool(x: np.ndarray) -> np.ndarray:
    """Global max pooling over sequence dimension.

    Args:
        x: 3D array (batch × seq_len × channels).

    Returns:
        2D array (batch × channels).
    """
    return x.max(axis=1)


def global_avg_pool(x: np.ndarray) -> np.ndarray:
    """Global average pooling over sequence dimension.

    Args:
        x: 3D array (batch × seq_len × channels).

    Returns:
        2D array (batch × channels).
    """
    return x.mean(axis=1)


# ──────────────────────────────────────────────────────────────
# Simple Sequence CNN
# ──────────────────────────────────────────────────────────────


@dataclass
class SequenceCNNConfig:
    """Configuration for a simple sequence CNN.

    Attributes:
        n_filters: List of filter counts per conv layer.
        kernel_sizes: List of kernel sizes per conv layer.
        n_classes: Number of output classes (1 for regression).
        alphabet: 'DNA' or 'RNA'.
        pool: Pooling strategy ('max' or 'avg').
    """

    n_filters: list[int] = field(default_factory=lambda: [32, 64])
    kernel_sizes: list[int] = field(default_factory=lambda: [8, 4])
    n_classes: int = 2
    alphabet: str = "DNA"
    pool: str = "max"


@dataclass
class SequenceCNNWeights:
    """Stored weights for a SequenceCNN.

    Attributes:
        conv_layers: List of Conv1DParams per layer.
        dense_weights: 2D array (features × n_classes).
        dense_bias: 1D array (n_classes,).
    """

    conv_layers: list[Conv1DParams]
    dense_weights: np.ndarray
    dense_bias: np.ndarray


def init_sequence_cnn_weights(
    config: SequenceCNNConfig,
    random_state: int = 42,
) -> SequenceCNNWeights:
    """Initialize random weights for a SequenceCNN.

    Uses He initialization for convolutional layers.

    Args:
        config: CNN configuration.
        random_state: Random seed.

    Returns:
        Initialized SequenceCNNWeights.
    """
    rng = np.random.RandomState(random_state)
    mapping = _DNA_MAP if config.alphabet.upper() == "DNA" else _RNA_MAP
    in_ch = len(mapping)

    conv_layers = []
    for i, (n_f, k_s) in enumerate(zip(config.n_filters, config.kernel_sizes)):
        scale = np.sqrt(2.0 / (in_ch * k_s))
        w = rng.randn(n_f, in_ch, k_s).astype(np.float32) * scale
        b = np.zeros(n_f, dtype=np.float32)
        conv_layers.append(Conv1DParams(weights=w, bias=b))
        in_ch = n_f

    # Dense layer
    n_features = config.n_filters[-1]
    scale = np.sqrt(2.0 / n_features)
    dw = rng.randn(n_features, config.n_classes).astype(np.float32) * scale
    db = np.zeros(config.n_classes, dtype=np.float32)

    return SequenceCNNWeights(
        conv_layers=conv_layers,
        dense_weights=dw,
        dense_bias=db,
    )


def predict_sequences(
    sequences: list[str],
    weights: SequenceCNNWeights,
    config: SequenceCNNConfig,
    max_length: int | None = None,
) -> np.ndarray:
    """Run inference on a batch of sequences.

    Args:
        sequences: List of nucleotide strings.
        weights: Trained model weights.
        config: Model configuration.
        max_length: Pad/truncate length.

    Returns:
        2D array (batch × n_classes) of predictions (logits or scores).
    """
    x = batch_encode(sequences, max_length=max_length, alphabet=config.alphabet)

    for conv_params in weights.conv_layers:
        x = conv1d_forward(x, conv_params)
        x = relu(x)

    if config.pool == "avg":
        pooled = global_avg_pool(x)
    else:
        pooled = global_max_pool(x)

    logits = pooled @ weights.dense_weights + weights.dense_bias
    return logits
