from __future__ import annotations

import pytest
from metainformant.dna import distances


def test_p_distance_and_jc69() -> None:
    """Test p-distance and Jukes-Cantor distance calculations."""
    s1, s2 = "AAAA", "AAAT"
    p = distances.p_distance(s1, s2)
    assert abs(p - 0.25) < 1e-9
    d = distances.jc69_distance(s1, s2)
    assert abs(d - 0.304099) < 1e-3


def test_kimura_2p_distance() -> None:
    """Test Kimura 2-parameter distance."""
    s1, s2 = "AAAA", "GGGG"  # All transitions
    d = distances.kimura_2p_distance(s1, s2)
    assert d > 0  # Should be positive

    s1, s2 = "AAAA", "TTTT"  # All transversions
    d = distances.kimura_2p_distance(s1, s2)
    assert d > 0  # Should be positive


def test_kmer_distance() -> None:
    """Test k-mer based distance calculations."""
    s1, s2 = "AAAA", "TTTT"
    d = distances.kmer_distance(s1, s2, k=2)
    assert 0 <= d <= 1  # Cosine distance should be in [0, 1]

    # Identical sequences should have distance 0
    d = distances.kmer_distance(s1, s1, k=2)
    assert abs(d) < 1e-9


def test_kmer_distance_matrix() -> None:
    """Test k-mer distance matrix calculation."""
    seqs = {"seq1": "AAAA", "seq2": "TTTT", "seq3": "AAAA"}
    matrix = distances.kmer_distance_matrix(seqs, k=2)

    assert len(matrix) == 3
    assert len(matrix[0]) == 3
    # Distance from seq1 to seq1 should be 0
    assert abs(matrix[0][0]) < 1e-9
    # Distance from seq1 to seq3 should be 0 (identical sequences)
    assert abs(matrix[0][2]) < 1e-9


def test_tajima_nei_distance() -> None:
    """Test Tajima-Nei distance."""
    s1, s2 = "AAAA", "AAAT"
    d = distances.tajima_nei_distance(s1, s2)
    assert d >= 0  # Distance should be non-negative


def test_sequence_identity_matrix() -> None:
    """Test sequence identity matrix calculation."""
    seqs = ["AAAA", "AAAT", "TTTT"]
    matrix = distances.sequence_identity_matrix(seqs)

    assert len(matrix) == 3
    assert len(matrix[0]) == 3
    # Self-identity should be 1.0
    assert abs(matrix[0][0] - 1.0) < 1e-9
    assert abs(matrix[1][1] - 1.0) < 1e-9
    assert abs(matrix[2][2] - 1.0) < 1e-9

    # AAAT vs AAAA should have identity 0.75 (3 out of 4 bases match)
    assert abs(matrix[0][1] - 0.75) < 1e-9


def test_distance_edge_cases() -> None:
    """Test distance calculations with edge cases."""
    # Empty sequences - returns 0.0 (graceful handling)
    assert distances.p_distance("", "A") == 0.0
    assert distances.p_distance("A", "") == 0.0
    assert distances.p_distance("", "") == 0.0

    # Different lengths - handled gracefully using min length
    assert distances.p_distance("AA", "AAA") == 0.0  # Compares "AA" vs "AA"
    assert distances.p_distance("AT", "ATG") == 0.0  # Compares "AT" vs "AT"

    # Identical sequences
    d = distances.p_distance("AAAA", "AAAA")
    assert d == 0.0
