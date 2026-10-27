"""Regression tests for sequences shorter than the convolution kernel in the NumPy sequence CNN.

``conv1d_forward`` computed ``out_len = (seq_len - kernel_size) // stride + 1`` without checking the input. A sequence
a few positions shorter than the kernel gave a negative length and NumPy's ``ValueError: negative dimensions are not
allowed``; one position shorter gave an empty output that blew up later in the pooling step as ``zero-size array to
reduction operation maximum which has no identity``. Neither says that the sequence is too short, let alone for which
kernel. With the default configuration (kernels 8 and 4) a sequence of 8 to 10 positions passed the first layer and
failed in the second.

The error now names the sequence length and the kernel size (``conv1d_forward``) or the minimum length the whole stack
needs (``predict_sequences``). The convolution oracle is a direct sliding-window computation written in the test. No
test doubles.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.lib.stride_tricks import sliding_window_view

from metainformant.ml.deep_learning.sequences import (
    Conv1DParams,
    SequenceCNNConfig,
    batch_encode,
    conv1d_forward,
    init_sequence_cnn_weights,
    predict_sequences,
)


def _params(out_channels: int = 3, in_channels: int = 4, kernel_size: int = 5, seed: int = 0) -> Conv1DParams:
    rng = np.random.default_rng(seed)
    return Conv1DParams(
        weights=rng.normal(size=(out_channels, in_channels, kernel_size)).astype(np.float32),
        bias=rng.normal(size=out_channels).astype(np.float32),
    )


@pytest.mark.parametrize("seq_len", [1, 3, 4])  # kernel size 5: one short, and further short
def test_input_shorter_than_the_kernel_gets_a_clear_message(seq_len: int) -> None:
    x = np.zeros((2, seq_len, 4), dtype=np.float32)

    with pytest.raises(ValueError, match=rf"length {seq_len} is shorter than the convolution kernel size 5"):
        conv1d_forward(x, _params(kernel_size=5))


def test_input_exactly_as_long_as_the_kernel_gives_one_output_position() -> None:
    params = _params(kernel_size=5)
    x = np.random.default_rng(1).normal(size=(2, 5, 4)).astype(np.float32)

    out = conv1d_forward(x, params)

    assert out.shape == (2, 1, 3)
    # oracle: out[b, 0, c] = sum over (position, channel) of x[b] * w[c] + bias[c]
    expected = np.einsum("bki,cik->bc", x, params.weights) + params.bias
    assert out[:, 0, :] == pytest.approx(expected, rel=1e-5, abs=1e-5)


def test_longer_input_matches_a_sliding_window_oracle() -> None:
    params = _params(kernel_size=3)
    x = np.random.default_rng(2).normal(size=(2, 9, 4)).astype(np.float32)

    out = conv1d_forward(x, params)

    windows = sliding_window_view(x, 3, axis=1)  # (batch, positions, channels, kernel)
    expected = np.einsum("bpik,cik->bpc", windows, params.weights) + params.bias
    assert out.shape == (2, 7, 3)
    assert out == pytest.approx(expected, rel=1e-5, abs=1e-5)


DEFAULT_CONFIG = SequenceCNNConfig()  # kernel sizes [8, 4]: the stack needs 8 + 4 - 1 = 11 positions
DEFAULT_WEIGHTS = init_sequence_cnn_weights(DEFAULT_CONFIG, random_state=0)


@pytest.mark.parametrize("length", [3, 7, 8, 10])
def test_predict_sequences_names_the_length_the_whole_stack_needs(length: int) -> None:
    sequences = ["ACGT" * 5 + "A"] * 2  # 21 bases, cut or padded to `length` below

    with pytest.raises(ValueError, match=rf"{length} positions long.*kernel sizes \[8, 4\].*at least 11"):
        predict_sequences(sequences, DEFAULT_WEIGHTS, DEFAULT_CONFIG, max_length=length)


def test_short_sequences_in_a_batch_are_reported_when_nothing_pads_them_up() -> None:
    with pytest.raises(ValueError, match=r"9 positions long.*at least 11"):
        predict_sequences(["ACGTACGTA", "ACG"], DEFAULT_WEIGHTS, DEFAULT_CONFIG)  # padded to the longest, 9


def test_the_shortest_sequence_the_network_accepts_is_predicted() -> None:
    sequences = ["ACGTACGTACG", "TTTTGGGGCCC"]  # 11 bases

    logits = predict_sequences(sequences, DEFAULT_WEIGHTS, DEFAULT_CONFIG)

    assert logits.shape == (2, DEFAULT_CONFIG.n_classes)
    assert np.all(np.isfinite(logits))
    # the same sequences zero-padded to 14 positions are accepted as well
    assert predict_sequences(sequences, DEFAULT_WEIGHTS, DEFAULT_CONFIG, max_length=14).shape == (2, 2)


def test_batch_encode_pads_to_the_requested_length() -> None:
    assert batch_encode(["ACG"], max_length=11).shape == (1, 11, 5)
