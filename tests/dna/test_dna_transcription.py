"""Tests for DNA transcription functions."""

from __future__ import annotations

import pytest

from metainformant.dna.expression import transcription


def test_transcribe_basic() -> None:
    """Test basic DNA to RNA transcription."""
    assert transcription.transcribe("ATGC") == "AUGC"


def test_transcribe_handles_lowercase_and_empty() -> None:
    """Test transcription handles edge cases."""
    assert transcription.transcribe("") == ""
    # Function converts to uppercase, so result is uppercase
    assert transcription.transcribe("atgc") == "AUGC"


def test_transcribe_rejects_invalid_characters() -> None:
    """Non-DNA characters must be rejected, not silently transcribed."""
    with pytest.raises(ValueError):
        transcription.transcribe("ATGX")


def test_transcribe_reverse_complement() -> None:
    """Reverse complement of ATGC (GCAT) transcribes to GCAU."""
    assert transcription.transcribe_reverse_complement("ATGC") == "GCAU"
    assert transcription.transcribe_reverse_complement("") == ""


def test_transcribe_with_introns_removes_intervals() -> None:
    """Intron intervals are excised before transcription."""
    dna = "AAATTTCCCGGG"
    # Remove TTT (3..6): exons AAA + CCCGGG
    assert transcription.transcribe_with_introns(dna, [(3, 6)]) == transcription.transcribe("AAACCCGGG")
    # Two introns
    assert transcription.transcribe_with_introns(dna, [(3, 6), (9, 12)]) == transcription.transcribe("AAACCC")
    # No introns -> plain transcription
    assert transcription.transcribe_with_introns(dna, []) == transcription.transcribe(dna)


def test_transcribe_with_introns_contained_interval() -> None:
    """A contained/unsorted intron must not re-emit already-removed DNA."""
    dna = "AAAGGGTTTAAACCC"
    # (0, 10) removes the first 10 bases; (2, 5) is contained in it.
    assert transcription.transcribe_with_introns(dna, [(0, 10), (2, 5)]) == transcription.transcribe("AACCC")


def test_transcribe_with_introns_rejects_malformed_intervals() -> None:
    """start > end (or out-of-bounds coordinates) must raise ValueError."""
    with pytest.raises(ValueError):
        transcription.transcribe_with_introns("AAACCC", [(5, 2)])
    with pytest.raises(ValueError):
        transcription.transcribe_with_introns("AAACCC", [(0, 10)])


def test_find_transcription_start_sites_offset() -> None:
    """TSS is reported len(pattern)+30 bp downstream of each promoter hit."""
    seq = "TATA" + "N" * 30 + "ATG"
    assert transcription.find_transcription_start_sites(seq) == [34]
    # Sequence too short for the offset: no hit reported.
    assert transcription.find_transcription_start_sites("TATAAA") == []
    # Case-insensitive.
    assert transcription.find_transcription_start_sites("tata" + "N" * 30 + "ATG") == [34]


def test_find_transcription_start_sites_custom_pattern() -> None:
    seq = "CAAT" + "N" * 30 + "G"
    assert transcription.find_transcription_start_sites(seq, promoter_pattern="CAAT") == [34]


def test_calculate_transcription_efficiency_bands() -> None:
    """Efficiency scoring bands are pinned: TATA 0.4, GC 0.3, CAAT 0.3."""
    # Short sequences score 0.
    assert transcription.calculate_transcription_efficiency("ATG") == 0.0
    # TATA only (uniform A promoter, GC content 0 -> no GC bonus).
    tata_only = "TATA" + "A" * 196
    assert transcription.calculate_transcription_efficiency(tata_only) == pytest.approx(0.4)
    # TATA + CAAT (inside the first 200 bp), GC content outside 0.4-0.6.
    assert transcription.calculate_transcription_efficiency("TATA" + "A" * 190 + "CAAT") == pytest.approx(0.7)
    # GC-rich promoter region (48/100 GC) with TATA: 0.4 + 0.3.
    gc_promoter = "TATA" + "G" * 48 + "A" * 48 + "C" * 10
    assert transcription.calculate_transcription_efficiency(gc_promoter) == pytest.approx(0.7)
