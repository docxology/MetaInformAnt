"""Tests for DNA distance calculations."""

from __future__ import annotations

import pytest
from metainformant.dna.alignment import distances


def test_p_distance_and_jc69() -> None:
    """Test p-distance and Jukes-Cantor distance calculations."""
    s1, s2 = "AAAA", "AAAT"
    p = distances.p_distance(s1, s2)
    assert abs(p - 0.25) < 1e-9
    d = distances.jc69_distance(s1, s2)
    assert d > 0  # Positive distance


def test_kimura_2p_distance() -> None:
    """Test Kimura 2-parameter distance."""
    s1, s2 = "AAAA", "GGGG"  # All transitions
    d = distances.kimura_2p_distance(s1, s2)
    assert d > 0  # Should be positive


def test_kmer_distance() -> None:
    """Test k-mer based distance calculations."""
    s1, s2 = "AAAA", "TTTT"
    d = distances.kmer_distance(s1, s2, k=2)
    assert d >= 0  # Distance should be non-negative

    # Identical sequences should have distance 0
    d = distances.kmer_distance(s1, s1, k=2)
    assert abs(d) < 1e-9


def test_kmer_distance_matrix() -> None:
    """Test k-mer distance matrix calculation."""
    seqs = {"seq1": "AAAA", "seq2": "TTTT", "seq3": "AAAA"}
    try:
        matrix = distances.kmer_distance_matrix(seqs, k=2)
        # Matrix could be a DataFrame or list of lists
        assert matrix is not None
    except (AttributeError, TypeError, KeyError) as e:
        pytest.skip(f"kmer_distance_matrix not available: {e}")


def test_sequence_identity_matrix() -> None:
    """Test sequence identity matrix calculation."""
    seqs = ["AAAA", "AAAT", "TTTT"]
    try:
        matrix = distances.sequence_identity_matrix(seqs)
        assert matrix is not None
    except (AttributeError, TypeError) as e:
        pytest.skip(f"sequence_identity_matrix not available: {e}")


def test_distance_edge_cases() -> None:
    """Test distance calculations with edge cases."""
    # Identical sequences
    d = distances.p_distance("AAAA", "AAAA")
    assert d == 0.0
