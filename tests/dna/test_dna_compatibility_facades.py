"""Tests for supported DNA sequence modules."""

from __future__ import annotations

from metainformant.dna.sequence.core import gc_content
from metainformant.dna.transcription import transcribe_dna_to_rna
from metainformant.dna.translation import translate_dna


def test_sequence_core_computes_real_outputs() -> None:
    """The supported sequence core computes real outputs."""
    assert gc_content("GCGC") == 1.0


def test_transcription_and_translation_facades_work() -> None:
    """Supported transcription and translation modules compute real outputs."""
    assert transcribe_dna_to_rna("ATGC") == "AUGC"
    assert translate_dna("ATGTAA", to_stop=True) == "M"
