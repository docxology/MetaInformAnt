"""Tests for DNA codon usage functions."""

from __future__ import annotations

import pytest

from metainformant.dna.expression import codon


def test_codon_counts_exact() -> None:
    """Each codon of the input is counted exactly once."""
    seq = "ATGAAATTTGGGCCC"
    counts = codon.codon_counts(seq)
    assert counts == {"ATG": 1, "AAA": 1, "TTT": 1, "GGG": 1, "CCC": 1}


def test_codon_counts_rejects_bad_length() -> None:
    with pytest.raises(ValueError):
        codon.codon_counts("ATGA")


def test_codon_counts_empty_sequence() -> None:
    assert codon.codon_counts("") == {}


def test_codon_frequencies_sum_to_one() -> None:
    seq = "ATGAAATTTGGGCCC"
    freqs = codon.codon_frequencies(seq)
    assert len(freqs) == 5
    assert sum(freqs.values()) == pytest.approx(1.0)
    assert freqs["ATG"] == pytest.approx(0.2)


def test_codon_usage_matches_frequencies() -> None:
    """codon_usage reports the same relative usage as codon_frequencies."""
    seq = "ATGAAATTTGGGCCC"
    assert codon.codon_usage(seq) == codon.codon_frequencies(seq)
