"""Tests for DNA sequence utility functions."""
from __future__ import annotations

from metainformant.dna import sequences


def test_reverse_complement_and_gc_content() -> None:
    """Test reverse complement and GC content calculation."""
    s = "ACGTACGTTT"
    rc = sequences.reverse_complement(s)
    assert rc == "AAACGTACGT"  # reverse complement of s
    gc = sequences.gc_content(s)
    assert abs(gc - 0.4) < 1e-9


def test_sequence_complexity() -> None:
    """Test sequence complexity calculation."""
    complex_seq = "ATCGATCGATCG"
    complexity = sequences.calculate_sequence_complexity(complex_seq)
    assert 0 <= complexity <= 1


def test_gc_content_edge_cases() -> None:
    """Test GC content with edge cases."""
    # All GC
    gc = sequences.gc_content("GCGCGC")
    assert gc == 1.0
    
    # No GC
    gc = sequences.gc_content("ATATAT")
    assert gc == 0.0
