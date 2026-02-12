"""DNA transcription utilities.

This module provides functions for transcribing DNA to RNA.
"""

from __future__ import annotations

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def transcribe(dna_seq: str) -> str:
    """Transcribe DNA sequence to RNA.

    Replaces T with U in DNA sequence.

    Args:
        dna_seq: DNA sequence string

    Returns:
        RNA sequence string

    Raises:
        ValueError: If sequence contains invalid characters
    """
    if not dna_seq:
        return ""

    dna_upper = dna_seq.upper()

    # Validate DNA characters
    valid_chars = set("ATCGN")
    if not all(c in valid_chars for c in dna_upper):
        invalid = set(dna_upper) - valid_chars
        raise ValueError(f"Invalid DNA characters: {invalid}")

    # Transcribe T -> U
    rna_seq = dna_upper.replace("T", "U")

    return rna_seq


def transcribe_reverse_complement(dna_seq: str) -> str:
    """Transcribe the reverse complement of DNA to RNA.

    Args:
        dna_seq: DNA sequence string

    Returns:
        RNA sequence from reverse complement
    """
    from metainformant.dna.sequence.core import reverse_complement

    rev_comp = reverse_complement(dna_seq)
    return transcribe(rev_comp)


def transcribe_with_introns(dna_seq: str, introns: list[tuple[int, int]]) -> str:
    """Transcribe DNA after removing introns.

    Args:
        dna_seq: DNA sequence string
        introns: List of (start, end) intron positions

    Returns:
        RNA sequence after intron removal
    """
    if not introns:
        return transcribe(dna_seq)

    # Remove introns from DNA
    exons = []
    prev_end = 0

    for start, end in sorted(introns):
        if start > prev_end:
            exons.append(dna_seq[prev_end:start])
        prev_end = end

    # Add final exon
    if prev_end < len(dna_seq):
        exons.append(dna_seq[prev_end:])

    exon_sequence = "".join(exons)
    return transcribe(exon_sequence)


def find_transcription_start_sites(dna_seq: str, promoter_pattern: str = "TATA") -> list[int]:
    """Find potential transcription start sites.

    Args:
        dna_seq: DNA sequence string
        promoter_pattern: Promoter motif to search for

    Returns:
        List of transcription start positions
    """
    positions = []

    # Look for promoter motifs
    pattern_upper = promoter_pattern.upper()
    seq_upper = dna_seq.upper()

    start = 0
    while True:
        pos = seq_upper.find(pattern_upper, start)
        if pos == -1:
            break

        # Assume TSS is ~25-35 bp downstream of TATA box
        tss_pos = pos + len(pattern_upper) + 30
        if tss_pos < len(dna_seq):
            positions.append(tss_pos)

        start = pos + 1

    return positions


def calculate_transcription_efficiency(dna_seq: str) -> float:
    """Calculate estimated transcription efficiency.

    This is a simplified model based on promoter strength and sequence features.

    Args:
        dna_seq: DNA sequence string

    Returns:
        Transcription efficiency score (0.0 to 1.0)
    """
    if not dna_seq or len(dna_seq) < 50:
        return 0.0

    seq_upper = dna_seq.upper()
    score = 0.0

    # Check for TATA box
    if "TATA" in seq_upper[:100]:
        score += 0.4

    # Check for GC content in promoter region
    promoter = seq_upper[:100]
    gc_count = promoter.count("G") + promoter.count("C")
    gc_content = gc_count / len(promoter) if promoter else 0

    if 0.4 <= gc_content <= 0.6:
        score += 0.3

    # Check for CAAT box
    if "CAAT" in seq_upper[:200]:
        score += 0.3

    return min(score, 1.0)  # Cap at 1.0
