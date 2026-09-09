"""Shared sequence primitives used across domains.

Hosts the canonical DNA complement tables, sequence validation, and
reverse-complement so that any domain (``dna``, ``metagenomics``,
``longread``, ``networks``, ...) can use them without a cross-domain
import (see ``scripts/quality/check_module_boundaries.py``).
"""

from __future__ import annotations

COMPLEMENT_TABLE = str.maketrans("ATCGatcg", "TAGCtagc")
COMPLEMENT_TABLE_UPPER = str.maketrans("ATCG", "TAGC")

_VALID_DNA_CHARS = frozenset("ATCGNUWSMKRYBDHVatcgnuwsmkrybdhv-")


def validate_dna_sequence(seq: str) -> bool:
    """Validate that a sequence contains only valid DNA characters.

    Args:
        seq: Sequence to validate

    Returns:
        True if sequence is valid DNA, False otherwise
    """
    if not seq:
        return False

    # Allow IUPAC ambiguity codes
    return all(c in _VALID_DNA_CHARS for c in seq)


def reverse_complement(seq: str) -> str:
    """Generate the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        Reverse complement sequence

    Raises:
        ValueError: If sequence contains invalid characters
    """
    if not seq:
        return ""

    # Validate sequence
    if not validate_dna_sequence(seq):
        raise ValueError(f"Invalid DNA sequence: {seq}")

    # Reverse and complement
    return seq.translate(COMPLEMENT_TABLE)[::-1]
