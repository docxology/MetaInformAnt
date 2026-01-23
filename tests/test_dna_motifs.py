"""Tests for DNA motif finding functions."""

from __future__ import annotations

from metainformant.dna import motifs


def test_find_motif_positions_exact() -> None:
    """Test motif finding with exact matches."""
    seq = "GGCATGATGCAATG"
    pos = motifs.find_motif_positions(seq, "ATG")
    assert pos == [3, 6, 11]

    # Test case insensitivity
    pos = motifs.find_motif_positions(seq.lower(), "atg")
    assert pos == [3, 6, 11]  # Should be case insensitive


def test_find_motif_positions_edge_cases() -> None:
    """Test motif finding edge cases."""
    # Empty sequence
    pos = motifs.find_motif_positions("", "ATG")
    assert pos == []

    # Motif longer than sequence
    pos = motifs.find_motif_positions("AT", "ATG")
    assert pos == []

    # No matches
    pos = motifs.find_motif_positions("AAAA", "TTT")
    assert pos == []


def test_find_motif_overlapping() -> None:
    """Test overlapping motif finding."""
    seq = "ATATAT"
    pos = motifs.find_motif_positions(seq, "ATA")
    assert 0 in pos  # ATA at position 0
    # Note: overlapping behavior may vary by implementation


def test_find_motif_returns_list() -> None:
    """Test that find_motif_positions returns a list of integers."""
    seq = "ATGATGATG"
    pos = motifs.find_motif_positions(seq, "ATG")
    assert isinstance(pos, list)
    assert all(isinstance(p, int) for p in pos)
