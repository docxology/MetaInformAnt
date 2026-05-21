"""Compatibility re-export for DNA transcription utilities."""

from __future__ import annotations

from metainformant.dna.expression.transcription import *  # noqa: F403
from metainformant.dna.expression.transcription import transcribe


def transcribe_dna_to_rna(dna_seq: str) -> str:
    """Transcribe DNA sequence to RNA."""
    return transcribe(dna_seq)
