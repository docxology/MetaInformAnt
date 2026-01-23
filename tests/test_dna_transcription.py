"""Tests for DNA transcription functions."""

from metainformant.dna.expression import transcription


def test_transcribe_basic() -> None:
    """Test basic DNA to RNA transcription."""
    assert transcription.transcribe("ATGC") == "AUGC"


def test_transcribe_handles_lowercase_and_empty() -> None:
    """Test transcription handles edge cases."""
    assert transcription.transcribe("") == ""
    # Function converts to uppercase, so result is uppercase
    assert transcription.transcribe("atgc") == "AUGC"


def test_transcribe_reverse_complement() -> None:
    """Test transcription of reverse complement."""
    # This function exists in the module
    result = transcription.transcribe_reverse_complement("ATGC")
    assert isinstance(result, str)
    assert len(result) > 0


def test_find_transcription_start_sites() -> None:
    """Test finding transcription start sites."""
    # Sequence with TATA box
    seq = "AAAAAATATAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGT"
    sites = transcription.find_transcription_start_sites(seq)
    assert isinstance(sites, list)


def test_calculate_transcription_efficiency() -> None:
    """Test transcription efficiency calculation."""
    # Short sequence should have low efficiency
    assert transcription.calculate_transcription_efficiency("ATG") == 0.0

    # Sequence with TATA box should have higher efficiency
    seq_with_tata = "TATA" + "A" * 200
    efficiency = transcription.calculate_transcription_efficiency(seq_with_tata)
    assert efficiency > 0.0
