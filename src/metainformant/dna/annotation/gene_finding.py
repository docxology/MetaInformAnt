"""Gene finding and ORF prediction utilities.

This module provides tools for predicting open reading frames (ORFs),
translating sequences, and related sequence utility functions used
across the gene annotation pipeline.
"""

from __future__ import annotations

import re
from typing import Any, Dict

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Standard genetic code (NCBI translation table 1) - DNA codons
GENETIC_CODES: Dict[int, Dict[str, str]] = {
    1: {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGT": "C",
        "TGC": "C",
        "TGA": "*",
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    },
}

# Start codons per genetic code
START_CODONS: Dict[int, set] = {
    1: {"ATG"},
}

# Stop codons per genetic code
STOP_CODONS: Dict[int, set] = {
    1: {"TAA", "TAG", "TGA"},
}

# IUPAC ambiguity codes for nucleotides
IUPAC_MAP: Dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}

# Consensus regulatory motifs with IUPAC codes
REGULATORY_MOTIFS: Dict[str, Dict[str, Any]] = {
    "TATA_box": {
        "pattern": "TATAAW",
        "description": "TATA box promoter element",
        "type": "promoter",
        "typical_position": (-30, -25),
    },
    "CAAT_box": {
        "pattern": "CCAAT",
        "description": "CAAT box promoter element",
        "type": "promoter",
        "typical_position": (-80, -70),
    },
    "GC_box": {
        "pattern": "GGGCGG",
        "description": "GC box (Sp1 binding site)",
        "type": "promoter",
        "typical_position": (-110, -40),
    },
    "Kozak_sequence": {
        "pattern": "GCCRMVATGG",
        "description": "Kozak consensus sequence for translation initiation",
        "type": "translation_initiation",
        "typical_position": (-6, 4),
    },
    "polyadenylation_signal": {
        "pattern": "AATAAA",
        "description": "Polyadenylation signal",
        "type": "terminator",
        "typical_position": (10, 30),
    },
    "CCAAT_enhancer": {
        "pattern": "CCAAT",
        "description": "CCAAT/enhancer binding element",
        "type": "enhancer",
        "typical_position": (-80, -60),
    },
    "E_box": {
        "pattern": "CANNTG",
        "description": "E-box enhancer element (bHLH binding)",
        "type": "enhancer",
        "typical_position": None,
    },
    "AP1_site": {
        "pattern": "TGASTCA",
        "description": "AP-1 transcription factor binding site",
        "type": "enhancer",
        "typical_position": None,
    },
    "CRE_element": {
        "pattern": "TGACGTCA",
        "description": "cAMP response element",
        "type": "enhancer",
        "typical_position": None,
    },
}


def _iupac_to_regex(pattern: str) -> str:
    """Convert an IUPAC nucleotide pattern to a regular expression.

    Args:
        pattern: IUPAC nucleotide pattern string.

    Returns:
        Regular expression string.
    """
    regex_parts = []
    for char in pattern.upper():
        if char in IUPAC_MAP:
            regex_parts.append(IUPAC_MAP[char])
        else:
            regex_parts.append(re.escape(char))
    return "".join(regex_parts)


def _reverse_complement(sequence: str) -> str:
    """Compute the reverse complement of a DNA sequence.

    Args:
        sequence: DNA sequence string.

    Returns:
        Reverse complement string.
    """
    complement_table = str.maketrans("ACGTacgt", "TGCAtgca")
    return sequence.translate(complement_table)[::-1]


def _translate_sequence(sequence: str, genetic_code: int = 1) -> str:
    """Translate a DNA sequence to a protein sequence.

    Args:
        sequence: DNA sequence (length must be divisible by 3).
        genetic_code: NCBI genetic code table number.

    Returns:
        Protein sequence string.
    """
    code = GENETIC_CODES.get(genetic_code, GENETIC_CODES[1])
    protein = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i : i + 3].upper()
        aa = code.get(codon, "X")
        protein.append(aa)
    return "".join(protein)


def predict_orfs(
    sequence: str,
    min_length: int = 100,
    genetic_code: int = 1,
) -> list[dict]:
    """Find all open reading frames across all 6 reading frames.

    Scans the forward and reverse complement strands across all three reading
    frames to identify ORFs bounded by start (ATG) and stop codons. Returns
    ORFs meeting the minimum nucleotide length threshold.

    Args:
        sequence: DNA sequence to scan for ORFs.
        min_length: Minimum ORF length in nucleotides (default 100).
        genetic_code: NCBI genetic code table number (default 1, standard).

    Returns:
        List of dicts, each containing:
            - frame: Reading frame (1, 2, 3 for forward; -1, -2, -3 for reverse)
            - start: Start position (0-based) in the original sequence
            - end: End position (0-based, exclusive) in the original sequence
            - length: Length in nucleotides
            - sequence: Nucleotide sequence of the ORF
            - protein: Translated protein sequence

    Raises:
        ValueError: If sequence is empty or contains invalid characters.
    """
    if not sequence:
        raise ValueError("Sequence must not be empty")

    seq_upper = sequence.upper()
    valid_bases = set("ACGTNRYSWKMBDHV")
    invalid = set(seq_upper) - valid_bases
    if invalid:
        raise ValueError(f"Invalid characters in sequence: {invalid}")

    start_codons = START_CODONS.get(genetic_code, START_CODONS[1])
    stop_codons = STOP_CODONS.get(genetic_code, STOP_CODONS[1])
    orfs: list[dict] = []

    strands = [
        (seq_upper, [1, 2, 3]),
        (_reverse_complement(seq_upper), [-1, -2, -3]),
    ]

    for strand_seq, frames in strands:
        for frame_offset, frame_label in zip(range(3), frames):
            i = frame_offset
            while i <= len(strand_seq) - 3:
                codon = strand_seq[i : i + 3]
                if codon in start_codons:
                    # Scan forward for stop codon
                    orf_start = i
                    j = i + 3
                    found_stop = False
                    while j <= len(strand_seq) - 3:
                        next_codon = strand_seq[j : j + 3]
                        if next_codon in stop_codons:
                            orf_end = j + 3  # Include stop codon
                            orf_length = orf_end - orf_start
                            if orf_length >= min_length:
                                orf_seq = strand_seq[orf_start:orf_end]
                                protein = _translate_sequence(orf_seq, genetic_code)
                                # Map coordinates back to original sequence for reverse strand
                                if frame_label > 0:
                                    orig_start = orf_start
                                    orig_end = orf_end
                                else:
                                    orig_start = len(seq_upper) - orf_end
                                    orig_end = len(seq_upper) - orf_start
                                orfs.append(
                                    {
                                        "frame": frame_label,
                                        "start": orig_start,
                                        "end": orig_end,
                                        "length": orf_length,
                                        "sequence": orf_seq,
                                        "protein": protein,
                                    }
                                )
                            found_stop = True
                            i = j + 3  # Continue scanning after stop
                            break
                        j += 3
                    if not found_stop:
                        i += 3
                else:
                    i += 3

    # Sort by length descending
    orfs.sort(key=lambda x: x["length"], reverse=True)
    logger.debug(f"Found {len(orfs)} ORFs (min_length={min_length}) in sequence of {len(sequence)} bp")
    return orfs
