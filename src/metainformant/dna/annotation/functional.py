"""Functional annotation utilities for DNA variants and protein sequences.

This module provides tools for classifying variant effects (synonymous,
nonsynonymous, nonsense, frameshift), predicting variant impact using
substitution matrices and physicochemical distances, computing positional
conservation from alignments, and identifying functional protein domains
via motif-based scanning.
"""

from __future__ import annotations

import math
import re
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Standard genetic code (DNA codons, NCBI table 1)
GENETIC_CODE: Dict[str, str] = {
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
}

# BLOSUM62 substitution matrix (subset of common amino acid pairs)
# Scores represent log-odds of substitution in related proteins
BLOSUM62: Dict[Tuple[str, str], int] = {}

_BLOSUM62_RAW = {
    "A": {
        "A": 4,
        "R": -1,
        "N": -2,
        "D": -2,
        "C": 0,
        "Q": -1,
        "E": -1,
        "G": 0,
        "H": -2,
        "I": -1,
        "L": -1,
        "K": -1,
        "M": -1,
        "F": -2,
        "P": -1,
        "S": 1,
        "T": 0,
        "W": -3,
        "Y": -2,
        "V": 0,
    },
    "R": {
        "R": 5,
        "N": 0,
        "D": -2,
        "C": -3,
        "Q": 1,
        "E": 0,
        "G": -2,
        "H": 0,
        "I": -3,
        "L": -2,
        "K": 2,
        "M": -1,
        "F": -3,
        "P": -2,
        "S": -1,
        "T": -1,
        "W": -3,
        "Y": -2,
        "V": -3,
    },
    "N": {
        "N": 6,
        "D": 1,
        "C": -3,
        "Q": 0,
        "E": 0,
        "G": 0,
        "H": 1,
        "I": -3,
        "L": -3,
        "K": 0,
        "M": -2,
        "F": -3,
        "P": -2,
        "S": 1,
        "T": 0,
        "W": -4,
        "Y": -2,
        "V": -3,
    },
    "D": {
        "D": 6,
        "C": -3,
        "Q": 0,
        "E": 2,
        "G": -1,
        "H": -1,
        "I": -3,
        "L": -4,
        "K": -1,
        "M": -3,
        "F": -3,
        "P": -1,
        "S": 0,
        "T": -1,
        "W": -4,
        "Y": -3,
        "V": -3,
    },
    "C": {
        "C": 9,
        "Q": -3,
        "E": -4,
        "G": -3,
        "H": -3,
        "I": -1,
        "L": -1,
        "K": -3,
        "M": -1,
        "F": -2,
        "P": -3,
        "S": -1,
        "T": -1,
        "W": -2,
        "Y": -2,
        "V": -1,
    },
    "Q": {
        "Q": 5,
        "E": 2,
        "G": -2,
        "H": 0,
        "I": -3,
        "L": -2,
        "K": 1,
        "M": 0,
        "F": -3,
        "P": -1,
        "S": 0,
        "T": -1,
        "W": -2,
        "Y": -1,
        "V": -2,
    },
    "E": {
        "E": 5,
        "G": -2,
        "H": 0,
        "I": -3,
        "L": -3,
        "K": 1,
        "M": -2,
        "F": -3,
        "P": -1,
        "S": 0,
        "T": -1,
        "W": -3,
        "Y": -2,
        "V": -2,
    },
    "G": {
        "G": 6,
        "H": -2,
        "I": -4,
        "L": -4,
        "K": -2,
        "M": -3,
        "F": -3,
        "P": -2,
        "S": 0,
        "T": -2,
        "W": -2,
        "Y": -3,
        "V": -3,
    },
    "H": {"H": 8, "I": -3, "L": -3, "K": -1, "M": -2, "F": -1, "P": -2, "S": -1, "T": -2, "W": -2, "Y": 2, "V": -3},
    "I": {"I": 4, "L": 2, "K": -3, "M": 1, "F": 0, "P": -3, "S": -2, "T": -1, "W": -3, "Y": -1, "V": 3},
    "L": {"L": 4, "K": -2, "M": 2, "F": 0, "P": -3, "S": -2, "T": -1, "W": -2, "Y": -1, "V": 1},
    "K": {"K": 5, "M": -1, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -2, "V": -2},
    "M": {"M": 5, "F": 0, "P": -2, "S": -1, "T": -1, "W": -1, "Y": -1, "V": 1},
    "F": {"F": 6, "P": -4, "S": -2, "T": -2, "W": 1, "Y": 3, "V": -1},
    "P": {"P": 7, "S": -1, "T": -1, "W": -4, "Y": -3, "V": -2},
    "S": {"S": 4, "T": 1, "W": -3, "Y": -2, "V": -2},
    "T": {"T": 5, "W": -2, "Y": -2, "V": 0},
    "W": {"W": 11, "Y": 2, "V": -3},
    "Y": {"Y": 7, "V": -1},
    "V": {"V": 4},
}

# Populate symmetric BLOSUM62 matrix
for aa1, scores in _BLOSUM62_RAW.items():
    for aa2, score in scores.items():
        BLOSUM62[(aa1, aa2)] = score
        BLOSUM62[(aa2, aa1)] = score

# Grantham distances between amino acid pairs (physicochemical property differences)
# Higher values = more different amino acids
GRANTHAM_DISTANCES: Dict[Tuple[str, str], int] = {
    ("S", "R"): 110,
    ("S", "L"): 145,
    ("S", "P"): 74,
    ("S", "T"): 58,
    ("S", "A"): 99,
    ("S", "V"): 124,
    ("S", "G"): 56,
    ("S", "I"): 142,
    ("S", "F"): 155,
    ("S", "Y"): 144,
    ("S", "C"): 112,
    ("S", "H"): 89,
    ("S", "Q"): 68,
    ("S", "N"): 46,
    ("S", "K"): 121,
    ("S", "D"): 65,
    ("S", "E"): 80,
    ("S", "M"): 135,
    ("S", "W"): 177,
    ("R", "L"): 102,
    ("R", "P"): 103,
    ("R", "T"): 71,
    ("R", "A"): 112,
    ("R", "V"): 96,
    ("R", "G"): 125,
    ("R", "I"): 97,
    ("R", "F"): 97,
    ("R", "Y"): 77,
    ("R", "C"): 180,
    ("R", "H"): 29,
    ("R", "Q"): 43,
    ("R", "N"): 86,
    ("R", "K"): 26,
    ("R", "D"): 96,
    ("R", "E"): 54,
    ("R", "M"): 91,
    ("R", "W"): 101,
    ("L", "P"): 98,
    ("L", "T"): 92,
    ("L", "A"): 96,
    ("L", "V"): 32,
    ("L", "G"): 138,
    ("L", "I"): 5,
    ("L", "F"): 22,
    ("L", "Y"): 36,
    ("L", "C"): 198,
    ("L", "H"): 99,
    ("L", "Q"): 113,
    ("L", "N"): 153,
    ("L", "K"): 107,
    ("L", "D"): 172,
    ("L", "E"): 138,
    ("L", "M"): 15,
    ("L", "W"): 61,
    ("P", "T"): 38,
    ("P", "A"): 27,
    ("P", "V"): 68,
    ("P", "G"): 42,
    ("P", "I"): 95,
    ("P", "F"): 114,
    ("P", "Y"): 110,
    ("P", "C"): 169,
    ("P", "H"): 77,
    ("P", "Q"): 76,
    ("P", "N"): 91,
    ("P", "K"): 103,
    ("P", "D"): 108,
    ("P", "E"): 93,
    ("P", "M"): 87,
    ("P", "W"): 147,
    ("T", "A"): 58,
    ("T", "V"): 69,
    ("T", "G"): 59,
    ("T", "I"): 89,
    ("T", "F"): 103,
    ("T", "Y"): 92,
    ("T", "C"): 149,
    ("T", "H"): 47,
    ("T", "Q"): 42,
    ("T", "N"): 65,
    ("T", "K"): 78,
    ("T", "D"): 85,
    ("T", "E"): 65,
    ("T", "M"): 81,
    ("T", "W"): 128,
    ("A", "V"): 64,
    ("A", "G"): 60,
    ("A", "I"): 94,
    ("A", "F"): 113,
    ("A", "Y"): 112,
    ("A", "C"): 195,
    ("A", "H"): 86,
    ("A", "Q"): 91,
    ("A", "N"): 111,
    ("A", "K"): 106,
    ("A", "D"): 126,
    ("A", "E"): 107,
    ("A", "M"): 84,
    ("A", "W"): 148,
    ("V", "G"): 109,
    ("V", "I"): 29,
    ("V", "F"): 50,
    ("V", "Y"): 55,
    ("V", "C"): 192,
    ("V", "H"): 84,
    ("V", "Q"): 96,
    ("V", "N"): 133,
    ("V", "K"): 97,
    ("V", "D"): 152,
    ("V", "E"): 121,
    ("V", "M"): 21,
    ("V", "W"): 88,
    ("G", "I"): 135,
    ("G", "F"): 153,
    ("G", "Y"): 147,
    ("G", "C"): 159,
    ("G", "H"): 98,
    ("G", "Q"): 87,
    ("G", "N"): 80,
    ("G", "K"): 127,
    ("G", "D"): 94,
    ("G", "E"): 98,
    ("G", "M"): 127,
    ("G", "W"): 184,
    ("I", "F"): 21,
    ("I", "Y"): 33,
    ("I", "C"): 198,
    ("I", "H"): 94,
    ("I", "Q"): 109,
    ("I", "N"): 149,
    ("I", "K"): 102,
    ("I", "D"): 168,
    ("I", "E"): 134,
    ("I", "M"): 10,
    ("I", "W"): 61,
    ("F", "Y"): 22,
    ("F", "C"): 205,
    ("F", "H"): 100,
    ("F", "Q"): 116,
    ("F", "N"): 158,
    ("F", "K"): 102,
    ("F", "D"): 177,
    ("F", "E"): 140,
    ("F", "M"): 28,
    ("F", "W"): 40,
    ("Y", "C"): 194,
    ("Y", "H"): 83,
    ("Y", "Q"): 99,
    ("Y", "N"): 143,
    ("Y", "K"): 85,
    ("Y", "D"): 160,
    ("Y", "E"): 122,
    ("Y", "M"): 36,
    ("Y", "W"): 37,
    ("C", "H"): 174,
    ("C", "Q"): 154,
    ("C", "N"): 139,
    ("C", "K"): 202,
    ("C", "D"): 154,
    ("C", "E"): 170,
    ("C", "M"): 196,
    ("C", "W"): 215,
    ("H", "Q"): 24,
    ("H", "N"): 68,
    ("H", "K"): 32,
    ("H", "D"): 81,
    ("H", "E"): 40,
    ("H", "M"): 87,
    ("H", "W"): 115,
    ("Q", "N"): 46,
    ("Q", "K"): 53,
    ("Q", "D"): 61,
    ("Q", "E"): 29,
    ("Q", "M"): 101,
    ("Q", "W"): 130,
    ("N", "K"): 94,
    ("N", "D"): 23,
    ("N", "E"): 42,
    ("N", "M"): 142,
    ("N", "W"): 174,
    ("K", "D"): 101,
    ("K", "E"): 56,
    ("K", "M"): 95,
    ("K", "W"): 110,
    ("D", "E"): 45,
    ("D", "M"): 160,
    ("D", "W"): 181,
    ("E", "M"): 126,
    ("E", "W"): 152,
    ("M", "W"): 67,
}

# Make symmetric
_symmetric_grantham: Dict[Tuple[str, str], int] = {}
for (aa1, aa2), dist in GRANTHAM_DISTANCES.items():
    _symmetric_grantham[(aa1, aa2)] = dist
    _symmetric_grantham[(aa2, aa1)] = dist
    _symmetric_grantham[(aa1, aa1)] = 0
    _symmetric_grantham[(aa2, aa2)] = 0
GRANTHAM_DISTANCES = _symmetric_grantham

# Kyte-Doolittle hydrophobicity scale
HYDROPHOBICITY: Dict[str, float] = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
}


def annotate_variants(
    variants: list[dict],
    reference: str,
) -> list[dict]:
    """Classify variants as synonymous, nonsynonymous, nonsense, or frameshift.

    Each variant is classified based on its effect on the coding sequence
    when compared to the reference. The reference is treated as a coding
    sequence (reading frame starts at position 0).

    Args:
        variants: List of variant dicts, each containing:
            - pos: 0-based position in the reference
            - ref: Reference allele string
            - alt: Alternate allele string
        reference: Reference coding DNA sequence.

    Returns:
        List of annotated variant dicts, each containing all input fields plus:
            - effect: One of "synonymous", "nonsynonymous", "nonsense",
              "frameshift", "start_loss", "stop_loss", "no_change"
            - ref_codon: Reference codon at the variant position
            - alt_codon: Alternate codon after applying the variant
            - ref_aa: Reference amino acid
            - alt_aa: Alternate amino acid
            - codon_position: Position within the codon (0, 1, or 2)

    Raises:
        ValueError: If reference is empty or variants list is None.
    """
    if not reference:
        raise ValueError("Reference sequence must not be empty")
    if variants is None:
        raise ValueError("Variants list must not be None")

    ref_upper = reference.upper()
    annotated: list[dict] = []

    for variant in variants:
        pos = variant.get("pos", 0)
        ref_allele = variant.get("ref", "").upper()
        alt_allele = variant.get("alt", "").upper()

        result = dict(variant)

        # Check for indel (frameshift candidate)
        ref_len = len(ref_allele)
        alt_len = len(alt_allele)

        if ref_len != alt_len and (ref_len - alt_len) % 3 != 0:
            result["effect"] = "frameshift"
            result["ref_codon"] = ""
            result["alt_codon"] = ""
            result["ref_aa"] = ""
            result["alt_aa"] = ""
            result["codon_position"] = pos % 3
            annotated.append(result)
            continue

        if ref_len != alt_len and (ref_len - alt_len) % 3 == 0:
            # In-frame indel (not frameshift but still a complex variant)
            result["effect"] = "inframe_indel"
            result["ref_codon"] = ""
            result["alt_codon"] = ""
            result["ref_aa"] = ""
            result["alt_aa"] = ""
            result["codon_position"] = pos % 3
            annotated.append(result)
            continue

        # SNP or MNP (same length substitution)
        # Apply variant to reference to get altered sequence
        alt_seq = ref_upper[:pos] + alt_allele + ref_upper[pos + ref_len :]

        # Determine affected codons
        codon_start = (pos // 3) * 3
        codon_end = ((pos + ref_len - 1) // 3) * 3 + 3

        ref_codons_region = ref_upper[codon_start:codon_end]
        alt_codons_region = alt_seq[codon_start:codon_end]

        # Translate both
        ref_aas = _translate_region(ref_codons_region)
        alt_aas = _translate_region(alt_codons_region)

        ref_codon = ref_upper[codon_start : codon_start + 3] if codon_start + 3 <= len(ref_upper) else ""
        alt_codon = alt_seq[codon_start : codon_start + 3] if codon_start + 3 <= len(alt_seq) else ""

        result["ref_codon"] = ref_codon
        result["alt_codon"] = alt_codon
        result["ref_aa"] = ref_aas
        result["alt_aa"] = alt_aas
        result["codon_position"] = pos % 3

        # Classify effect
        if ref_aas == alt_aas:
            result["effect"] = "synonymous" if ref_codon != alt_codon else "no_change"
        elif "*" in alt_aas and "*" not in ref_aas:
            result["effect"] = "nonsense"
        elif "*" in ref_aas and "*" not in alt_aas:
            result["effect"] = "stop_loss"
        elif ref_codon == "ATG" and alt_codon != "ATG" and codon_start == 0:
            result["effect"] = "start_loss"
        else:
            result["effect"] = "nonsynonymous"

        annotated.append(result)

    logger.debug(f"Annotated {len(annotated)} variants")
    return annotated


def _translate_region(dna_region: str) -> str:
    """Translate a DNA region to amino acids.

    Args:
        dna_region: DNA sequence (length should be divisible by 3).

    Returns:
        Amino acid sequence string.
    """
    protein = []
    for i in range(0, len(dna_region) - 2, 3):
        codon = dna_region[i : i + 3]
        aa = GENETIC_CODE.get(codon, "X")
        protein.append(aa)
    return "".join(protein)


def predict_variant_impact(
    variant: dict,
    context: dict | None = None,
) -> dict:
    """Score variant impact using BLOSUM62, Grantham distance, and conservation.

    Combines multiple scoring methods to produce an overall impact assessment.
    Higher scores indicate more damaging variants.

    Args:
        variant: Variant dict containing at minimum:
            - ref_aa: Reference amino acid (single letter)
            - alt_aa: Alternate amino acid (single letter)
            - effect: Variant effect classification
        context: Optional context dict containing:
            - conservation: Conservation score at this position (0.0-1.0)
            - domain: Whether position is in a functional domain (bool)

    Returns:
        Dict containing:
            - blosum62_score: BLOSUM62 substitution score
            - grantham_distance: Grantham physicochemical distance
            - impact_score: Combined impact score (0.0-1.0, higher = more damaging)
            - impact_category: "benign", "possibly_damaging", or "probably_damaging"
            - conservation_weight: Conservation contribution (if provided)
            - status: "success"

    Raises:
        ValueError: If required variant fields are missing.
    """
    ref_aa = variant.get("ref_aa", "")
    alt_aa = variant.get("alt_aa", "")
    effect = variant.get("effect", "unknown")

    if not ref_aa or not alt_aa:
        # Cannot score without amino acid information
        return {
            "blosum62_score": 0,
            "grantham_distance": 0,
            "impact_score": 1.0 if effect in ("nonsense", "frameshift") else 0.0,
            "impact_category": "probably_damaging" if effect in ("nonsense", "frameshift") else "benign",
            "conservation_weight": 0.0,
            "status": "success",
        }

    # Use only first amino acid if multiple
    ref_aa_single = ref_aa[0] if ref_aa else "X"
    alt_aa_single = alt_aa[0] if alt_aa else "X"

    # BLOSUM62 score
    blosum_score = BLOSUM62.get((ref_aa_single, alt_aa_single), -4)
    # Normalize: typical range is -4 to +11, map to 0-1 (inverted: negative = damaging)
    blosum_normalized = max(0.0, min(1.0, (4 - blosum_score) / 15.0))

    # Grantham distance
    grantham = GRANTHAM_DISTANCES.get((ref_aa_single, alt_aa_single), 100)
    # Normalize: range 0-215, higher = more different
    grantham_normalized = min(1.0, grantham / 215.0)

    # Conservation weight
    conservation_score = 0.0
    if context and "conservation" in context:
        conservation_score = float(context["conservation"])

    domain_bonus = 0.0
    if context and context.get("domain", False):
        domain_bonus = 0.15

    # Combined impact score
    # Weight: 35% BLOSUM, 35% Grantham, 20% conservation, 10% domain
    if conservation_score > 0:
        impact_score = (
            0.35 * blosum_normalized
            + 0.35 * grantham_normalized
            + 0.20 * conservation_score
            + 0.10 * (1.0 if domain_bonus > 0 else 0.0)
        )
    else:
        # Without conservation data, weight BLOSUM and Grantham more
        impact_score = 0.45 * blosum_normalized + 0.45 * grantham_normalized + 0.10 * (1.0 if domain_bonus > 0 else 0.0)

    impact_score = min(1.0, impact_score + domain_bonus)

    # Nonsense and frameshift are always probably damaging
    if effect in ("nonsense", "frameshift", "start_loss"):
        impact_score = max(impact_score, 0.9)

    # Categorize
    if impact_score >= 0.7:
        category = "probably_damaging"
    elif impact_score >= 0.4:
        category = "possibly_damaging"
    else:
        category = "benign"

    logger.debug(
        f"Variant impact: {ref_aa_single}->{alt_aa_single}, "
        f"BLOSUM62={blosum_score}, Grantham={grantham}, "
        f"score={impact_score:.3f} ({category})"
    )

    return {
        "blosum62_score": blosum_score,
        "grantham_distance": grantham,
        "impact_score": round(impact_score, 4),
        "impact_category": category,
        "conservation_weight": round(conservation_score, 4),
        "status": "success",
    }


def compute_conservation_score(
    alignment: list[str],
    position: int,
) -> float:
    """Compute Shannon entropy-based conservation at an alignment position.

    Lower entropy indicates higher conservation. The returned score is
    normalized to 0.0-1.0 where 1.0 means perfectly conserved (entropy = 0)
    and 0.0 means maximally variable.

    Args:
        alignment: List of aligned sequences (same length, may contain gaps "-").
        position: Column position (0-based) in the alignment.

    Returns:
        Conservation score (0.0-1.0, where 1.0 = perfectly conserved).

    Raises:
        ValueError: If alignment is empty, position is out of range, or
            sequences have different lengths.
    """
    if not alignment:
        raise ValueError("Alignment must not be empty")

    seq_length = len(alignment[0])
    for seq in alignment:
        if len(seq) != seq_length:
            raise ValueError("All sequences in alignment must have the same length")

    if position < 0 or position >= seq_length:
        raise ValueError(f"Position {position} out of range for alignment of length {seq_length}")

    # Extract column
    column = [seq[position].upper() for seq in alignment]

    # Count residues (excluding gaps)
    residue_counts: Dict[str, int] = defaultdict(int)
    total = 0
    for residue in column:
        if residue != "-" and residue != ".":
            residue_counts[residue] += 1
            total += 1

    if total == 0:
        return 0.0  # All gaps

    # Shannon entropy
    entropy = 0.0
    for count in residue_counts.values():
        freq = count / total
        if freq > 0:
            entropy -= freq * math.log2(freq)

    # Maximum possible entropy for the observed alphabet size
    # For DNA: max = log2(4) = 2.0; for protein: max = log2(20) = 4.32
    unique_residues = len(residue_counts)
    if unique_residues <= 1:
        return 1.0  # Perfectly conserved

    max_entropy = math.log2(min(unique_residues, 20))  # Cap at protein alphabet

    # Normalize: 0 entropy = 1.0 conservation, max entropy = 0.0 conservation
    if max_entropy > 0:
        conservation = 1.0 - (entropy / max_entropy)
    else:
        conservation = 1.0

    return round(max(0.0, min(1.0, conservation)), 4)


def identify_functional_domains(
    protein_sequence: str,
) -> list[dict]:
    """Identify functional domains in a protein sequence using motif scanning.

    Scans for common structural and functional motifs including:
    - Zinc finger domains (C2H2 pattern: C-X{2,4}-C-X{3}-[LIVMFYWC]-X{8}-H-X{3,5}-H)
    - Leucine zipper domains (L-X{6}-L repeats)
    - Signal peptides (hydrophobic N-terminal region)
    - Transmembrane helices (hydrophobic stretches via Kyte-Doolittle)
    - Nuclear localization signals (KKKRK-like patterns)
    - N-glycosylation sites (N-X-[ST] where X != P)

    Args:
        protein_sequence: Amino acid sequence (single-letter codes).

    Returns:
        List of dicts, each containing:
            - domain: Domain name
            - start: Start position (0-based)
            - end: End position (0-based, exclusive)
            - sequence: Matched sequence
            - score: Confidence score (0.0-1.0)
            - description: Human-readable description

    Raises:
        ValueError: If protein_sequence is empty.
    """
    if not protein_sequence:
        raise ValueError("Protein sequence must not be empty")

    seq = protein_sequence.upper()
    domains: list[dict] = []

    # 1. Zinc finger C2H2 motif
    zf_pattern = re.compile(r"C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H")
    for match in zf_pattern.finditer(seq):
        domains.append(
            {
                "domain": "zinc_finger_C2H2",
                "start": match.start(),
                "end": match.end(),
                "sequence": match.group(),
                "score": 0.9,
                "description": "C2H2-type zinc finger domain",
            }
        )

    # 2. Leucine zipper (L-X6 repeated 4+ times)
    for i in range(len(seq) - 27):
        leucine_positions = []
        for j in range(4):
            pos = i + j * 7
            if pos < len(seq) and seq[pos] == "L":
                leucine_positions.append(pos)

        if len(leucine_positions) >= 4:
            end_pos = min(leucine_positions[-1] + 7, len(seq))
            region = seq[i:end_pos]
            # Check it's not already captured
            overlap = any(d["domain"] == "leucine_zipper" and d["start"] <= i < d["end"] for d in domains)
            if not overlap:
                domains.append(
                    {
                        "domain": "leucine_zipper",
                        "start": i,
                        "end": end_pos,
                        "sequence": region,
                        "score": 0.7,
                        "description": "Leucine zipper dimerization motif",
                    }
                )

    # 3. Signal peptide (hydrophobic N-terminal 15-30 residues)
    if len(seq) >= 15:
        signal_region = seq[:30] if len(seq) >= 30 else seq
        hydro_scores = [HYDROPHOBICITY.get(aa, 0.0) for aa in signal_region]

        # Check for hydrophobic core (positions ~2-15)
        if len(hydro_scores) >= 15:
            core_start = 1
            core_end = min(15, len(hydro_scores))
            core_avg = sum(hydro_scores[core_start:core_end]) / (core_end - core_start)

            if core_avg > 1.5 and seq[0] == "M":
                # Look for signal peptidase cleavage site (small-X-small pattern)
                cleavage_pos = core_end
                for k in range(core_end, min(core_end + 10, len(seq) - 2)):
                    if seq[k] in "AGS" and seq[k + 2] in "AGS":
                        cleavage_pos = k + 3
                        break

                domains.append(
                    {
                        "domain": "signal_peptide",
                        "start": 0,
                        "end": cleavage_pos,
                        "sequence": seq[:cleavage_pos],
                        "score": 0.6 + min(0.3, (core_avg - 1.5) / 3.0),
                        "description": "N-terminal signal peptide",
                    }
                )

    # 4. Transmembrane helices (hydrophobic stretches of 18-25 residues)
    window_size = 21
    threshold = 1.6
    if len(seq) >= window_size:
        in_tm = False
        tm_start = 0

        for i in range(len(seq) - window_size + 1):
            window = seq[i : i + window_size]
            avg_hydro = sum(HYDROPHOBICITY.get(aa, 0.0) for aa in window) / window_size

            if avg_hydro > threshold and not in_tm:
                in_tm = True
                tm_start = i
            elif avg_hydro <= threshold and in_tm:
                tm_end = i + window_size - 1
                tm_length = tm_end - tm_start
                if 18 <= tm_length <= 30:
                    domains.append(
                        {
                            "domain": "transmembrane_helix",
                            "start": tm_start,
                            "end": tm_end,
                            "sequence": seq[tm_start:tm_end],
                            "score": 0.75,
                            "description": "Predicted transmembrane helix",
                        }
                    )
                in_tm = False

    # 5. Nuclear localization signal (NLS)
    # Monopartite: K-K/R-X-K/R or similar clusters
    nls_pattern = re.compile(r"[KR]{4,}")
    for match in nls_pattern.finditer(seq):
        if match.end() - match.start() >= 4:
            domains.append(
                {
                    "domain": "nuclear_localization_signal",
                    "start": match.start(),
                    "end": match.end(),
                    "sequence": match.group(),
                    "score": 0.65,
                    "description": "Putative nuclear localization signal",
                }
            )

    # Bipartite NLS: KR-X{10,12}-K[KR]{2}
    bipartite_nls = re.compile(r"[KR]{2}.{10,12}[KR][KR]{2}")
    for match in bipartite_nls.finditer(seq):
        domains.append(
            {
                "domain": "bipartite_nls",
                "start": match.start(),
                "end": match.end(),
                "sequence": match.group(),
                "score": 0.6,
                "description": "Putative bipartite nuclear localization signal",
            }
        )

    # 6. N-glycosylation sites (N-X-[ST] where X != P)
    for i in range(len(seq) - 2):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in "ST":
            domains.append(
                {
                    "domain": "n_glycosylation_site",
                    "start": i,
                    "end": i + 3,
                    "sequence": seq[i : i + 3],
                    "score": 0.8,
                    "description": "Potential N-linked glycosylation site (N-X-S/T)",
                }
            )

    # Sort by start position
    domains.sort(key=lambda x: (x["start"], x["domain"]))
    logger.debug(f"Identified {len(domains)} functional domains/motifs in {len(seq)} aa sequence")
    return domains
