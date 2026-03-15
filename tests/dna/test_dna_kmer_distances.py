"""Tests for DNA k-mer distances."""

from __future__ import annotations

from metainformant.dna.alignment import distances


def test_kmer_distance_identical_zero() -> None:
    """Test k-mer distance is zero for identical sequences."""
    a = "ACGTACGT"
    b = "ACGTACGT"
    # kmer_distance takes k param only (no metric)
    d_ab = distances.kmer_distance(a, b, k=2)
    assert abs(d_ab) < 1e-12


def test_kmer_distance_different_positive() -> None:
    """Test k-mer distance is positive for different sequences."""
    a = "ACGTACGT"
    c = "TTTTTTTT"
    d_ac = distances.kmer_distance(a, c, k=2)
    assert d_ac > 0


def test_kmer_distance_symmetric() -> None:
    """Test k-mer distance is symmetric."""
    a = "ACGTACGT"
    b = "AAAACCCC"
    d_ab = distances.kmer_distance(a, b, k=2)
    d_ba = distances.kmer_distance(b, a, k=2)
    assert abs(d_ab - d_ba) < 1e-12
