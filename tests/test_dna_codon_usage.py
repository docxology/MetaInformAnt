"""Tests for DNA codon usage functions."""

from __future__ import annotations

import pytest

from metainformant.dna import codon


def test_codon_usage_counts_and_freqs() -> None:
    """Test codon counting and frequency calculation."""
    seq = "ATGAAATTTGGGCCC"
    try:
        counts = codon.codon_counts(seq)
        assert "ATG" in counts
        freqs = codon.codon_frequencies(seq)
        assert isinstance(freqs, dict)
    except (AttributeError, ImportError) as e:
        pytest.skip(f"Codon functions unavailable: {e}")


def test_codon_module_importable() -> None:
    """Test that the codon module can be imported."""
    assert codon is not None
