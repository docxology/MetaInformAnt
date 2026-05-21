"""Compatibility re-export for DNA translation utilities."""

from __future__ import annotations

from metainformant.dna.expression.translation import *  # noqa: F403
from metainformant.dna.expression.translation import translate_dna as _translate_dna


def translate_dna(dna_seq: str, genetic_code: int = 1, *, to_stop: bool = False) -> str:
    """Translate DNA to protein, optionally truncating at the first stop codon."""
    protein = _translate_dna(dna_seq, genetic_code=genetic_code)
    if to_stop:
        return protein.split("*", 1)[0]
    return protein
