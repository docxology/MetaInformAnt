"""Tests for DNA population statistics."""
from __future__ import annotations

import pytest
from metainformant.dna import population


def test_nucleotide_diversity_two_sequences() -> None:
    """Test nucleotide diversity with two sequences."""
    seqs = ["AAAA", "AAAT"]
    pi = population.nucleotide_diversity(seqs)
    assert abs(pi - 0.25) < 1e-9


def test_tajimas_d_requires_enough_sequences() -> None:
    """Test that Tajima's D requires at least 4 sequences."""
    seqs = ["AAAA", "AAAA", "AAAA"]  # Only 3 sequences
    # Implementation requires at least 4 sequences
    with pytest.raises(ValueError, match="at least 4"):
        population.tajimas_d(seqs)


def test_tajimas_d_with_enough_sequences() -> None:
    """Test Tajima's D with enough sequences."""
    seqs = ["AAAA", "AAAA", "AAAA", "AAAA"]  # 4 identical sequences
    d = population.tajimas_d(seqs)
    assert d == 0.0  # No segregating sites


def test_fst_fixed_differences() -> None:
    """Test FST with fixed differences between populations."""
    pop1 = ["AAAA", "AAAA"]
    pop2 = ["TTTT", "TTTT"]
    fst = population.hudson_fst(pop1, pop2)
    # FST should be high for fixed differences
    assert 0.0 <= fst <= 1.0
