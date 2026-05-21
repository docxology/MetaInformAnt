"""Compatibility tests for historic DNA facade modules."""

from __future__ import annotations

from metainformant.dna.sequence.core import gc_content as canonical_gc_content
from metainformant.dna.sequence.core import read_fasta as canonical_read_fasta
from metainformant.dna.sequences import gc_content, read_fasta
from metainformant.dna.transcription import transcribe_dna_to_rna
from metainformant.dna.translation import translate_dna


def test_sequences_facade_reexports_sequence_core() -> None:
    """The old dna.sequences path remains a real compatibility facade."""
    assert gc_content is canonical_gc_content
    assert read_fasta is canonical_read_fasta
    assert gc_content("GCGC") == 1.0


def test_transcription_and_translation_facades_work() -> None:
    """Historic transcription and translation facades still compute real outputs."""
    assert transcribe_dna_to_rna("ATGC") == "AUGC"
    assert translate_dna("ATGTAA", to_stop=True) == "M"
