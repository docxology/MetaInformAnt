"""Tests for deep learning sequence encoding and NumPy CNN inference.

Pure NumPy implementations - real implementations, no torch required.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.ml.deep_learning.sequences import (
    Conv1DParams,
    SequenceCNNConfig,
    batch_encode,
    conv1d_forward,
    global_avg_pool,
    global_max_pool,
    init_sequence_cnn_weights,
    one_hot_encode,
    predict_sequences,
    relu,
)


class TestOneHotEncode:
    """Tests for one_hot_encode."""

    def test_dna_shape_and_values(self):
        """DNA encoding produces (len, 5) rows with exactly one hot value."""
        encoded = one_hot_encode("ACGT")

        assert encoded.shape == (4, 5)
        assert encoded.dtype == np.float32
        for i, base in enumerate("ACGT"):
            assert encoded[i].sum() == 1.0
            assert encoded[i].argmax() == i  # A=0, C=1, G=2, T=3

    def test_unknown_base_falls_back_to_N(self):
        """Bases outside the alphabet map to the N channel."""
        encoded = one_hot_encode("AXC")

        assert encoded[0].argmax() == 0  # A
        assert encoded[1].argmax() == 4  # X -> N channel
        assert encoded[2].argmax() == 1  # C

    def test_lowercase_normalised(self):
        """Lowercase input is uppercased before encoding."""
        assert np.array_equal(one_hot_encode("acgt"), one_hot_encode("ACGT"))

    def test_rna_alphabet(self):
        """RNA encoding uses ACGU with 5 channels."""
        encoded = one_hot_encode("ACGU", alphabet="RNA")

        assert encoded.shape == (4, 5)
        for i, base in enumerate("ACGU"):
            assert encoded[i].argmax() == i  # A=0, C=1, G=2, U=3


class TestBatchEncode:
    """Tests for batch_encode."""

    def test_pads_to_longest_sequence(self):
        """Without max_length, batch pads to the longest sequence."""
        batch = batch_encode(["ACG", "A"])

        assert batch.shape == (2, 3, 5)
        # Second sequence is padded with zeros beyond its length
        assert batch[1, 0].sum() == 1.0  # 'A' encoded
        assert batch[1, 1:].sum() == 0.0

    def test_truncates_to_max_length(self):
        """Sequences longer than max_length are truncated."""
        batch = batch_encode(["ACGTACGT"], max_length=4)

        assert batch.shape == (1, 4, 5)

    def test_explicit_max_length_pads(self):
        """Explicit max_length pads shorter sequences."""
        batch = batch_encode(["A"], max_length=3)

        assert batch.shape == (1, 3, 5)
        assert batch[0, 1:].sum() == 0.0


class TestConv1D:
    """Tests for conv1d_forward, activations, and pooling."""

    def test_output_shape(self):
        """Conv output shape follows (batch, out_len, out_ch)."""
        weights = np.ones((2, 3, 2), dtype=np.float32)
        bias = np.array([0.5, -0.5], dtype=np.float32)
        x = np.ones((3, 7, 3), dtype=np.float32)

        out = conv1d_forward(x, Conv1DParams(weights=weights, bias=bias), stride=1)

        assert out.shape == (3, 6, 2)
        # All-ones weights and input: each window sums to 6 (+ bias)
        assert np.allclose(out[:, :, 0], 6.5)
        assert np.allclose(out[:, :, 1], 5.5)

    def test_stride(self):
        """Stride reduces output length."""
        weights = np.ones((1, 1, 1), dtype=np.float32)
        bias = np.zeros(1, dtype=np.float32)
        x = np.ones((1, 6, 1), dtype=np.float32)

        out = conv1d_forward(x, Conv1DParams(weights=weights, bias=bias), stride=2)

        assert out.shape == (1, 3, 1)

    def test_relu(self):
        """ReLU zeroes negatives and keeps positives."""
        x = np.array([-1.0, 0.0, 2.5])
        assert np.array_equal(relu(x), np.array([0.0, 0.0, 2.5]))

    def test_global_max_pool(self):
        """Max pooling reduces (batch, seq, channels) to (batch, channels)."""
        x = np.array([[[1.0, 5.0], [3.0, 2.0]]])
        assert np.array_equal(global_max_pool(x), np.array([[3.0, 5.0]]))

    def test_global_avg_pool(self):
        """Average pooling reduces (batch, seq, channels) to (batch, channels)."""
        x = np.array([[[1.0, 5.0], [3.0, 2.0]]])
        assert np.allclose(global_avg_pool(x), np.array([[2.0, 3.5]]))


class TestSequenceCNN:
    """Tests for weight initialization and end-to-end inference."""

    def test_init_weights_shapes(self):
        """Initialized weights match the configured layer sizes."""
        config = SequenceCNNConfig(n_filters=[4, 8], kernel_sizes=[3, 2], n_classes=2)
        weights = init_sequence_cnn_weights(config)

        assert len(weights.conv_layers) == 2
        assert weights.conv_layers[0].weights.shape == (4, 5, 3)  # 5 = DNA channels
        assert weights.conv_layers[1].weights.shape == (8, 4, 2)
        assert weights.dense_weights.shape == (8, 2)
        assert weights.dense_bias.shape == (2,)

    def test_init_weights_reproducible(self):
        """Same random_state yields identical weights."""
        config = SequenceCNNConfig(n_filters=[2], kernel_sizes=[3])
        w1 = init_sequence_cnn_weights(config, random_state=7)
        w2 = init_sequence_cnn_weights(config, random_state=7)

        assert np.array_equal(w1.conv_layers[0].weights, w2.conv_layers[0].weights)
        assert np.array_equal(w1.dense_weights, w2.dense_weights)

    def test_predict_sequences_shapes(self):
        """Inference returns (batch, n_classes) logits."""
        config = SequenceCNNConfig(n_filters=[2], kernel_sizes=[3], n_classes=2)
        weights = init_sequence_cnn_weights(config, random_state=42)

        logits = predict_sequences(["ACGT", "TTTT"], weights, config, max_length=8)

        assert logits.shape == (2, 2)

    def test_predict_sequences_deterministic(self):
        """Same weights and input produce identical predictions."""
        config = SequenceCNNConfig(n_filters=[2], kernel_sizes=[3], n_classes=3)
        weights = init_sequence_cnn_weights(config, random_state=0)

        logits_a = predict_sequences(["ACGTAC"], weights, config)
        logits_b = predict_sequences(["ACGTAC"], weights, config)

        assert np.array_equal(logits_a, logits_b)

    def test_avg_pooling_variant(self):
        """The 'avg' pooling strategy runs end-to-end."""
        config = SequenceCNNConfig(n_filters=[2], kernel_sizes=[3], pool="avg")
        weights = init_sequence_cnn_weights(config, random_state=1)

        logits = predict_sequences(["ACGT"], weights, config)

        assert logits.shape == (1, 2)


def test_batch_encode_empty_raises():
    """Empty batch has no max length to infer and must fail loudly."""
    with pytest.raises(ValueError):
        batch_encode([])
