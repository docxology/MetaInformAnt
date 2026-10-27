"""Tests for DNA translation functions."""

from __future__ import annotations

from metainformant.dna.expression import translation


def test_translate_basic() -> None:
    """Test basic DNA translation."""
    # ATG = M (Met), GCC = A (Ala)
    assert translation.translate_dna("ATGGCC") == "MA"


def test_translate_empty() -> None:
    """Test translation of empty sequence."""
    assert translation.translate_dna("") == ""


def test_translate_rna() -> None:
    """Test RNA translation."""
    # Direct RNA translation
    assert translation.translate("AUGGCC") == "MA"


def test_find_orfs_simple() -> None:
    """Each AUG opens an ORF that ends at the first in-frame stop codon."""
    rna_seq = "AUGAAAUUUUGAAUGCCCCUAG"  # AUG AAA UUU UGA | AUG CCC CUA G
    orfs = translation.find_orfs(rna_seq, min_length=2)
    assert orfs == [(0, 12, "+1"), (12, 22, "+1")]


def test_find_orfs_reports_each_orf_once_in_its_own_frame() -> None:
    # The only start is at index 1 (frame +2); unaligned matching used to
    # report it again as frame +1.
    assert translation.find_orfs("CAUGAAAUAA", min_length=1) == [(1, 10, "+2")]
    rna = "AUGAAAUUUUAAAUGGGCCCCUAGAUGAAAUUUAAA"
    assert translation.find_orfs(rna, min_length=3) == [(0, 12, "+1"), (12, 24, "+1"), (24, 36, "+1")]


def test_find_start_codons() -> None:
    """Test finding start codons."""
    rna_seq = "CCAUGAAAUGCCC"
    positions = translation.find_start_codons(rna_seq)
    assert 2 in positions  # First AUG
    assert 7 in positions  # Second AUG


def test_find_stop_codons() -> None:
    """Test finding stop codons."""
    rna_seq = "CCUAAUAGUGACC"
    positions = translation.find_stop_codons(rna_seq)
    assert len(positions) >= 1  # Should find UAA, UAG, or UGA


def test_get_genetic_code() -> None:
    """Test getting genetic code table."""
    code = translation.get_genetic_code()
    assert isinstance(code, dict)
    assert code.get("AUG") == "M"  # Standard start codon
    assert code.get("UAA") == "*"  # Stop codon


def test_back_translate() -> None:
    """Test back-translation of protein to DNA."""
    dna = translation.back_translate("M")
    assert "ATG" in dna  # Met codon (with U->T conversion)
