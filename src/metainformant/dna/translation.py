"""Compatibility re-export for DNA translation utilities."""

from __future__ import annotations

from metainformant.dna.expression.translation import (  # noqa: F401
    back_translate,
    calculate_cai,
    find_orfs,
    find_start_codons,
    find_stop_codons,
    get_genetic_code,
    optimize_codons,
    six_frame_translation,
    translate,
    translate_dna as _translate_dna,
)


def translate_dna(dna_seq: str, genetic_code: int = 1, *, to_stop: bool = False) -> str:
    """Translate DNA to protein, optionally truncating at the first stop codon."""
    protein = _translate_dna(dna_seq, genetic_code=genetic_code)
    if to_stop:
        return protein.split("*", 1)[0]
    return protein
