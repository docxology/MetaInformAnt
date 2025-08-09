from metainformant.dna import sequences, transcription


def test_transcribe_dna_to_rna_basic() -> None:
    assert transcription.transcribe_dna_to_rna("ATGC") == "AUGC"


def test_transcribe_handles_lowercase_and_empty() -> None:
    assert transcription.transcribe_dna_to_rna("") == ""
    assert transcription.transcribe_dna_to_rna("atgc") == "augc"


def test_reverse_transcribe_rna_to_dna() -> None:
    assert transcription.reverse_transcribe_rna_to_dna("AUGC") == "ATGC"


