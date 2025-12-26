from __future__ import annotations

from metainformant.dna import motifs


def test_find_motif_positions_iupac() -> None:
    """Test motif finding with IUPAC codes."""
    seq = "GGCATGATGCAATG"
    # ATN should match ATG and ATA etc
    pos = motifs.find_motif_positions(seq, "ATN")
    # positions of 'ATG' occur at index 3 and 6 and 10; ATN also includes index 11? (ATG at 3, 6, 10)
    assert 3 in pos and 6 in pos and 10 in pos


def test_find_motif_positions_exact() -> None:
    """Test motif finding with exact matches."""
    seq = "GGCATGATGCAATG"
    pos = motifs.find_motif_positions(seq, "ATG")
    assert pos == [3, 6, 11]

    # Test case insensitivity
    pos = motifs.find_motif_positions(seq.lower(), "atg")
    assert pos == [3, 6, 11]  # Should be case insensitive


def test_find_motif_positions_iupac_variations() -> None:
    """Test motif finding with various IUPAC codes."""
    seq = "GGCATGATGCAATG"

    # Test R (A or G)
    pos = motifs.find_motif_positions(seq, "RG")
    assert 1 in pos  # GG at position 1

    # Test Y (C or T)
    pos = motifs.find_motif_positions(seq, "TY")
    assert 8 in pos  # CA at position 8

    # Test N (any nucleotide)
    pos = motifs.find_motif_positions(seq, "NNN")
    assert len(pos) > 0  # Should find many matches


def test_find_motif_positions_edge_cases() -> None:
    """Test motif finding edge cases."""
    # Empty sequence
    pos = motifs.find_motif_positions("", "ATG")
    assert pos == []

    # Motif longer than sequence
    pos = motifs.find_motif_positions("AT", "ATG")
    assert pos == []

    # No matches
    pos = motifs.find_motif_positions("AAAA", "TTTT")
    assert pos == []

    # Overlapping matches
    seq = "ATATAT"
    pos = motifs.find_motif_positions(seq, "ATA")
    assert 0 in pos  # ATA at position 0
    assert 2 in pos  # ATA at position 2 (overlapping)
