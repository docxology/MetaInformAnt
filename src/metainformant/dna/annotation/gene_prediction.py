"""Gene annotation and prediction utilities.

This module provides tools for predicting open reading frames (ORFs), annotating
coding regions, finding regulatory elements, masking repeats, annotating CpG islands,
computing codon usage statistics, and identifying splice sites in DNA sequences.
"""

from __future__ import annotations

import math
import re
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core import logging

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


def annotate_coding_regions(
    sequence: str,
    orfs: list[dict] | None = None,
) -> dict:
    """Classify regions of a sequence as coding or non-coding.

    Uses ORF predictions and codon usage statistics (hexamer scoring) to
    distinguish coding from non-coding regions. If ORFs are not provided,
    they are predicted with default parameters.

    Args:
        sequence: DNA sequence to annotate.
        orfs: Pre-computed ORFs (optional; predicted if not provided).

    Returns:
        Dict containing:
            - regions: List of {start, end, type, score} dicts
            - coding_fraction: Fraction of sequence in coding regions
            - total_orfs: Number of ORFs used
            - hexamer_scores: Per-region hexamer log-likelihood scores
            - status: "success"

    Raises:
        ValueError: If sequence is empty.
    """
    if not sequence:
        raise ValueError("Sequence must not be empty")

    seq_upper = sequence.upper()

    if orfs is None:
        orfs = predict_orfs(seq_upper, min_length=90)

    # Build coding mask
    coding_mask = [False] * len(seq_upper)
    for orf in orfs:
        start = orf["start"]
        end = orf["end"]
        for idx in range(max(0, start), min(end, len(seq_upper))):
            coding_mask[idx] = True

    # Compute hexamer scores for coding vs non-coding discrimination
    hexamer_scores = _compute_hexamer_scores(seq_upper, coding_mask)

    # Build contiguous regions
    regions: list[dict] = []
    if not coding_mask:
        return {
            "regions": [],
            "coding_fraction": 0.0,
            "total_orfs": 0,
            "hexamer_scores": {},
            "status": "success",
        }

    current_type = "coding" if coding_mask[0] else "non-coding"
    region_start = 0

    for i in range(1, len(coding_mask)):
        new_type = "coding" if coding_mask[i] else "non-coding"
        if new_type != current_type:
            region_seq = seq_upper[region_start:i]
            score = hexamer_scores.get(f"{region_start}-{i}", 0.0)
            regions.append(
                {
                    "start": region_start,
                    "end": i,
                    "type": current_type,
                    "score": round(score, 4),
                }
            )
            current_type = new_type
            region_start = i

    # Final region
    region_seq = seq_upper[region_start : len(seq_upper)]
    score = hexamer_scores.get(f"{region_start}-{len(seq_upper)}", 0.0)
    regions.append(
        {
            "start": region_start,
            "end": len(seq_upper),
            "type": current_type,
            "score": round(score, 4),
        }
    )

    coding_bp = sum(1 for x in coding_mask if x)
    coding_fraction = coding_bp / len(seq_upper) if seq_upper else 0.0

    logger.debug(
        f"Annotated {len(regions)} regions: {coding_fraction:.1%} coding " f"({coding_bp}/{len(seq_upper)} bp)"
    )

    return {
        "regions": regions,
        "coding_fraction": round(coding_fraction, 4),
        "total_orfs": len(orfs),
        "hexamer_scores": hexamer_scores,
        "status": "success",
    }


def _compute_hexamer_scores(
    sequence: str,
    coding_mask: list[bool],
) -> Dict[str, float]:
    """Compute hexamer log-likelihood ratio scores for regions.

    Uses coding vs non-coding hexamer frequency distributions to score
    regions. Positive scores suggest coding; negative suggest non-coding.

    Args:
        sequence: DNA sequence.
        coding_mask: Boolean mask indicating coding positions.

    Returns:
        Dict mapping "start-end" region keys to hexamer scores.
    """
    if len(sequence) < 6:
        return {}

    # Build hexamer frequency tables from coding and non-coding regions
    coding_hexamers: Dict[str, int] = defaultdict(int)
    noncoding_hexamers: Dict[str, int] = defaultdict(int)
    coding_total = 0
    noncoding_total = 0

    for i in range(len(sequence) - 5):
        hexamer = sequence[i : i + 6]
        if any(c not in "ACGT" for c in hexamer):
            continue
        if coding_mask[i]:
            coding_hexamers[hexamer] += 1
            coding_total += 1
        else:
            noncoding_hexamers[hexamer] += 1
            noncoding_total += 1

    # Score each contiguous region
    scores: Dict[str, float] = {}
    if coding_total == 0 or noncoding_total == 0:
        return scores

    # Identify regions
    current_type = coding_mask[0] if coding_mask else False
    region_start = 0

    for i in range(1, len(coding_mask)):
        if coding_mask[i] != current_type:
            key = f"{region_start}-{i}"
            scores[key] = _score_region_hexamers(
                sequence[region_start:i],
                coding_hexamers,
                noncoding_hexamers,
                coding_total,
                noncoding_total,
            )
            current_type = coding_mask[i]
            region_start = i

    key = f"{region_start}-{len(coding_mask)}"
    scores[key] = _score_region_hexamers(
        sequence[region_start : len(coding_mask)],
        coding_hexamers,
        noncoding_hexamers,
        coding_total,
        noncoding_total,
    )

    return scores


def _score_region_hexamers(
    region_seq: str,
    coding_hexamers: Dict[str, int],
    noncoding_hexamers: Dict[str, int],
    coding_total: int,
    noncoding_total: int,
) -> float:
    """Score a single region using hexamer log-likelihood ratio.

    Args:
        region_seq: Sequence of the region to score.
        coding_hexamers: Hexamer counts from coding regions.
        noncoding_hexamers: Hexamer counts from non-coding regions.
        coding_total: Total coding hexamers.
        noncoding_total: Total non-coding hexamers.

    Returns:
        Log-likelihood ratio score (positive = coding-like).
    """
    if len(region_seq) < 6:
        return 0.0

    score = 0.0
    count = 0
    pseudocount = 1.0

    for i in range(len(region_seq) - 5):
        hexamer = region_seq[i : i + 6]
        if any(c not in "ACGT" for c in hexamer):
            continue

        coding_freq = (coding_hexamers.get(hexamer, 0) + pseudocount) / (coding_total + pseudocount * 4096)
        noncoding_freq = (noncoding_hexamers.get(hexamer, 0) + pseudocount) / (noncoding_total + pseudocount * 4096)

        if noncoding_freq > 0:
            score += math.log(coding_freq / noncoding_freq)
            count += 1

    return score / count if count > 0 else 0.0


def find_regulatory_elements(
    sequence: str,
    elements: list[str] | None = None,
) -> list[dict]:
    """Search for regulatory elements using consensus motif patterns.

    Scans the sequence for known promoter elements (TATA box, CAAT box, GC box),
    enhancer elements (E-box, AP-1, CRE), and other regulatory motifs using
    IUPAC-aware pattern matching.

    Args:
        sequence: DNA sequence to scan.
        elements: List of element names to search for (default: all known motifs).
            Valid names: TATA_box, CAAT_box, GC_box, Kozak_sequence,
            polyadenylation_signal, CCAAT_enhancer, E_box, AP1_site, CRE_element.

    Returns:
        List of dicts, each containing:
            - name: Element name
            - type: Element type (promoter, enhancer, terminator, etc.)
            - start: Start position (0-based)
            - end: End position (0-based, exclusive)
            - sequence: Matched sequence
            - strand: "+" or "-"
            - description: Human-readable description
            - score: Match confidence (1.0 for exact match)

    Raises:
        ValueError: If sequence is empty or element names are invalid.
    """
    if not sequence:
        raise ValueError("Sequence must not be empty")

    seq_upper = sequence.upper()
    results: list[dict] = []

    if elements is None:
        motifs_to_search = REGULATORY_MOTIFS
    else:
        motifs_to_search = {}
        for name in elements:
            if name not in REGULATORY_MOTIFS:
                raise ValueError(
                    f"Unknown regulatory element: {name}. " f"Valid elements: {list(REGULATORY_MOTIFS.keys())}"
                )
            motifs_to_search[name] = REGULATORY_MOTIFS[name]

    for name, motif_info in motifs_to_search.items():
        pattern = _iupac_to_regex(motif_info["pattern"])
        regex = re.compile(pattern, re.IGNORECASE)

        # Search forward strand
        for match in regex.finditer(seq_upper):
            results.append(
                {
                    "name": name,
                    "type": motif_info["type"],
                    "start": match.start(),
                    "end": match.end(),
                    "sequence": match.group(),
                    "strand": "+",
                    "description": motif_info["description"],
                    "score": 1.0,
                }
            )

        # Search reverse strand
        rev_comp = _reverse_complement(seq_upper)
        for match in regex.finditer(rev_comp):
            # Convert reverse-strand coordinates to forward-strand
            fwd_start = len(seq_upper) - match.end()
            fwd_end = len(seq_upper) - match.start()
            results.append(
                {
                    "name": name,
                    "type": motif_info["type"],
                    "start": fwd_start,
                    "end": fwd_end,
                    "sequence": match.group(),
                    "strand": "-",
                    "description": motif_info["description"],
                    "score": 1.0,
                }
            )

    # Sort by position
    results.sort(key=lambda x: (x["start"], x["name"]))
    logger.debug(f"Found {len(results)} regulatory elements in {len(sequence)} bp sequence")
    return results


def mask_repeats(
    sequence: str,
    min_length: int = 10,
    min_copies: int = 2,
) -> dict:
    """Identify and soft-mask tandem repeats and simple sequence repeats.

    Scans for tandem repeats (exact adjacent copies of a unit) and simple
    sequence repeats (SSRs/microsatellites) of unit lengths 1-6 bp. Masked
    positions are converted to lowercase in the output.

    Args:
        sequence: DNA sequence to scan for repeats.
        min_length: Minimum total repeat region length in bp (default 10).
        min_copies: Minimum number of repeat copies (default 2).

    Returns:
        Dict containing:
            - masked_sequence: Sequence with repeats soft-masked (lowercase)
            - repeats: List of {start, end, unit, copies, length} dicts
            - total_masked_bp: Total bases masked
            - masked_fraction: Fraction of sequence that is masked
            - status: "success"

    Raises:
        ValueError: If sequence is empty or parameters are invalid.
    """
    if not sequence:
        raise ValueError("Sequence must not be empty")
    if min_length < 1:
        raise ValueError("min_length must be >= 1")
    if min_copies < 2:
        raise ValueError("min_copies must be >= 2")

    seq_upper = sequence.upper()
    repeats: list[dict] = []
    mask = [False] * len(seq_upper)

    # Scan for tandem repeats with unit sizes 1-6 (microsatellites)
    for unit_size in range(1, 7):
        i = 0
        while i <= len(seq_upper) - unit_size * min_copies:
            unit = seq_upper[i : i + unit_size]
            if any(c not in "ACGT" for c in unit):
                i += 1
                continue

            # Count consecutive copies
            copies = 1
            j = i + unit_size
            while j + unit_size <= len(seq_upper):
                if seq_upper[j : j + unit_size] == unit:
                    copies += 1
                    j += unit_size
                else:
                    break

            total_length = copies * unit_size
            if copies >= min_copies and total_length >= min_length:
                # Check for overlap with existing longer-unit repeats
                already_masked = all(mask[k] for k in range(i, i + total_length))
                if not already_masked:
                    repeats.append(
                        {
                            "start": i,
                            "end": i + total_length,
                            "unit": unit,
                            "copies": copies,
                            "length": total_length,
                        }
                    )
                    for k in range(i, i + total_length):
                        mask[k] = True
                i = j  # Skip past this repeat
            else:
                i += 1

    # Build masked sequence (lowercase for masked positions)
    masked_chars = []
    for idx, char in enumerate(sequence):
        if mask[idx]:
            masked_chars.append(char.lower())
        else:
            masked_chars.append(char.upper())
    masked_sequence = "".join(masked_chars)

    total_masked = sum(mask)
    masked_fraction = total_masked / len(sequence) if sequence else 0.0

    # Sort repeats by position
    repeats.sort(key=lambda x: x["start"])

    logger.debug(
        f"Masked {total_masked}/{len(sequence)} bp ({masked_fraction:.1%}) " f"across {len(repeats)} repeat regions"
    )

    return {
        "masked_sequence": masked_sequence,
        "repeats": repeats,
        "total_masked_bp": total_masked,
        "masked_fraction": round(masked_fraction, 4),
        "status": "success",
    }


def annotate_cpg_islands(
    sequence: str,
    min_length: int = 200,
    min_gc: float = 0.5,
    min_obs_exp: float = 0.6,
) -> list[dict]:
    """Find CpG islands using Gardiner-Garden and Frommer criteria.

    A CpG island is defined as a region of at least `min_length` bp with
    GC content >= `min_gc` and observed/expected CpG ratio >= `min_obs_exp`.

    The sliding window approach scans the sequence and merges overlapping
    regions that meet all three criteria.

    Args:
        sequence: DNA sequence to scan.
        min_length: Minimum island length in bp (default 200).
        min_gc: Minimum GC fraction (default 0.5).
        min_obs_exp: Minimum observed/expected CpG ratio (default 0.6).

    Returns:
        List of dicts, each containing:
            - start: Start position (0-based)
            - end: End position (0-based, exclusive)
            - length: Island length in bp
            - gc_content: GC fraction in the island
            - obs_exp_cpg: Observed/expected CpG ratio
            - cpg_count: Number of CpG dinucleotides
            - sequence: Island sequence

    Raises:
        ValueError: If sequence is empty or parameters are out of range.
    """
    if not sequence:
        raise ValueError("Sequence must not be empty")
    if min_length < 1:
        raise ValueError("min_length must be >= 1")
    if not 0.0 <= min_gc <= 1.0:
        raise ValueError("min_gc must be between 0.0 and 1.0")
    if not 0.0 <= min_obs_exp <= 2.0:
        raise ValueError("min_obs_exp must be between 0.0 and 2.0")

    seq_upper = sequence.upper()

    if len(seq_upper) < min_length:
        return []

    # Sliding window scan
    window_size = min_length
    step = max(1, window_size // 10)
    candidate_regions: list[Tuple[int, int]] = []

    for start in range(0, len(seq_upper) - window_size + 1, step):
        window = seq_upper[start : start + window_size]
        gc_content = (window.count("G") + window.count("C")) / len(window)

        if gc_content < min_gc:
            continue

        cpg_count = window.count("CG")
        c_count = window.count("C")
        g_count = window.count("G")

        # Observed/expected CpG = (CpG * N) / (C * G)
        if c_count > 0 and g_count > 0:
            obs_exp = (cpg_count * len(window)) / (c_count * g_count)
        else:
            obs_exp = 0.0

        if obs_exp >= min_obs_exp:
            candidate_regions.append((start, start + window_size))

    if not candidate_regions:
        return []

    # Merge overlapping candidate regions
    merged = _merge_intervals(candidate_regions)

    # Extend merged regions to maximize length while maintaining criteria
    islands: list[dict] = []
    for region_start, region_end in merged:
        # Try to extend in both directions
        best_start, best_end = _extend_cpg_region(seq_upper, region_start, region_end, min_gc, min_obs_exp)
        region_length = best_end - best_start

        if region_length < min_length:
            continue

        island_seq = seq_upper[best_start:best_end]
        gc_content = (island_seq.count("G") + island_seq.count("C")) / len(island_seq)
        cpg_count = island_seq.count("CG")
        c_count = island_seq.count("C")
        g_count = island_seq.count("G")

        if c_count > 0 and g_count > 0:
            obs_exp = (cpg_count * len(island_seq)) / (c_count * g_count)
        else:
            obs_exp = 0.0

        # Final check against criteria
        if gc_content >= min_gc and obs_exp >= min_obs_exp:
            islands.append(
                {
                    "start": best_start,
                    "end": best_end,
                    "length": region_length,
                    "gc_content": round(gc_content, 4),
                    "obs_exp_cpg": round(obs_exp, 4),
                    "cpg_count": cpg_count,
                    "sequence": island_seq,
                }
            )

    # Deduplicate overlapping islands (keep longer)
    islands = _deduplicate_islands(islands)
    islands.sort(key=lambda x: x["start"])

    logger.debug(f"Found {len(islands)} CpG islands in {len(sequence)} bp sequence")
    return islands


def _merge_intervals(intervals: list[Tuple[int, int]]) -> list[Tuple[int, int]]:
    """Merge overlapping or adjacent intervals.

    Args:
        intervals: List of (start, end) tuples.

    Returns:
        Merged list of (start, end) tuples.
    """
    if not intervals:
        return []

    sorted_intervals = sorted(intervals)
    merged = [sorted_intervals[0]]

    for start, end in sorted_intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    return merged


def _extend_cpg_region(
    sequence: str,
    start: int,
    end: int,
    min_gc: float,
    min_obs_exp: float,
) -> Tuple[int, int]:
    """Extend a CpG region while maintaining criteria.

    Args:
        sequence: Full DNA sequence.
        start: Current region start.
        end: Current region end.
        min_gc: Minimum GC content.
        min_obs_exp: Minimum observed/expected CpG ratio.

    Returns:
        Extended (start, end) tuple.
    """
    step = 10
    best_start = start
    best_end = end

    # Extend left
    new_start = start
    while new_start - step >= 0:
        new_start -= step
        region = sequence[new_start:best_end]
        gc = (region.count("G") + region.count("C")) / len(region)
        cpg = region.count("CG")
        c = region.count("C")
        g = region.count("G")
        obs_exp = (cpg * len(region)) / (c * g) if c > 0 and g > 0 else 0.0

        if gc >= min_gc and obs_exp >= min_obs_exp:
            best_start = new_start
        else:
            break

    # Extend right
    new_end = end
    while new_end + step <= len(sequence):
        new_end += step
        region = sequence[best_start:new_end]
        gc = (region.count("G") + region.count("C")) / len(region)
        cpg = region.count("CG")
        c = region.count("C")
        g = region.count("G")
        obs_exp = (cpg * len(region)) / (c * g) if c > 0 and g > 0 else 0.0

        if gc >= min_gc and obs_exp >= min_obs_exp:
            best_end = new_end
        else:
            break

    return best_start, best_end


def _deduplicate_islands(islands: list[dict]) -> list[dict]:
    """Remove overlapping islands, keeping the longer one.

    Args:
        islands: List of island dicts with start/end keys.

    Returns:
        Deduplicated list of islands.
    """
    if len(islands) <= 1:
        return islands

    # Sort by length descending to prioritize longer islands
    sorted_islands = sorted(islands, key=lambda x: x["length"], reverse=True)
    kept: list[dict] = []
    used_positions: set = set()

    for island in sorted_islands:
        positions = set(range(island["start"], island["end"]))
        overlap = positions & used_positions
        if len(overlap) < len(positions) * 0.5:  # Less than 50% overlap
            kept.append(island)
            used_positions.update(positions)

    return kept


def compute_codon_usage(
    sequence: str,
    genetic_code: int = 1,
) -> dict:
    """Compute codon usage table with RSCU values and Codon Adaptation Index.

    Calculates the frequency of each codon, its Relative Synonymous Codon
    Usage (RSCU), and the overall Codon Adaptation Index (CAI) for the sequence.

    RSCU = (observed frequency of codon) / (expected frequency if all synonymous
    codons were used equally).

    Args:
        sequence: DNA coding sequence (should be divisible by 3).
        genetic_code: NCBI genetic code table number (default 1).

    Returns:
        Dict containing:
            - codon_counts: Dict mapping codons to raw counts
            - codon_frequencies: Dict mapping codons to frequencies
            - rscu: Dict mapping codons to RSCU values
            - amino_acid_counts: Dict mapping amino acids to total counts
            - cai: Codon Adaptation Index value (0.0-1.0)
            - total_codons: Total number of codons analyzed
            - gc3: GC content at third codon position
            - status: "success"

    Raises:
        ValueError: If sequence is empty.
    """
    if not sequence:
        raise ValueError("Sequence must not be empty")

    seq_upper = sequence.upper()
    code = GENETIC_CODES.get(genetic_code, GENETIC_CODES[1])

    # Count codons
    codon_counts: Dict[str, int] = defaultdict(int)
    total_codons = 0
    gc3_count = 0
    total_third_pos = 0

    for i in range(0, len(seq_upper) - 2, 3):
        codon = seq_upper[i : i + 3]
        if len(codon) == 3 and codon in code:
            codon_counts[codon] += 1
            total_codons += 1
            # Track GC at third position
            if codon[2] in "GC":
                gc3_count += 1
            total_third_pos += 1

    if total_codons == 0:
        return {
            "codon_counts": {},
            "codon_frequencies": {},
            "rscu": {},
            "amino_acid_counts": {},
            "cai": 0.0,
            "total_codons": 0,
            "gc3": 0.0,
            "status": "success",
        }

    # Compute frequencies
    codon_frequencies = {codon: count / total_codons for codon, count in codon_counts.items()}

    # Group codons by amino acid
    aa_to_codons: Dict[str, list[str]] = defaultdict(list)
    for codon, aa in code.items():
        aa_to_codons[aa].append(codon)

    # Compute amino acid counts
    amino_acid_counts: Dict[str, int] = defaultdict(int)
    for codon, count in codon_counts.items():
        aa = code.get(codon, "X")
        amino_acid_counts[aa] += count

    # Compute RSCU
    rscu: Dict[str, float] = {}
    for aa, synonymous in aa_to_codons.items():
        total_aa_count = sum(codon_counts.get(c, 0) for c in synonymous)
        n_synonymous = len(synonymous)
        for codon in synonymous:
            observed = codon_counts.get(codon, 0)
            expected = total_aa_count / n_synonymous if n_synonymous > 0 else 0
            rscu[codon] = observed / expected if expected > 0 else 0.0

    # Compute CAI using RSCU values
    cai_values: list[float] = []
    for codon, count in codon_counts.items():
        aa = code.get(codon, "X")
        if aa == "*":
            continue
        synonymous = aa_to_codons.get(aa, [])
        max_rscu = max(rscu.get(c, 0) for c in synonymous) if synonymous else 0
        if max_rscu > 0:
            w = rscu.get(codon, 0) / max_rscu
            if w > 0:
                cai_values.extend([math.log(w)] * count)

    cai_value = math.exp(sum(cai_values) / len(cai_values)) if cai_values else 0.0
    gc3 = gc3_count / total_third_pos if total_third_pos > 0 else 0.0

    logger.debug(f"Computed codon usage for {total_codons} codons, CAI={cai_value:.3f}")

    return {
        "codon_counts": dict(codon_counts),
        "codon_frequencies": {k: round(v, 6) for k, v in codon_frequencies.items()},
        "rscu": {k: round(v, 4) for k, v in rscu.items()},
        "amino_acid_counts": dict(amino_acid_counts),
        "cai": round(cai_value, 4),
        "total_codons": total_codons,
        "gc3": round(gc3, 4),
        "status": "success",
    }


def find_splice_sites(
    sequence: str,
    model: str = "consensus",
) -> list[dict]:
    """Identify potential GT-AG splice donor and acceptor sites.

    Searches for canonical splice site dinucleotides (GT for donor, AG for
    acceptor) and scores them based on the surrounding consensus sequence
    context. Supports the consensus scoring model.

    Donor consensus: [C/A]AG|GT[A/G]AGT (| marks the splice site)
    Acceptor consensus: [C/T]{10,}N[C/T]AG| (polypyrimidine tract + AG)

    Args:
        sequence: DNA sequence to scan for splice sites.
        model: Scoring model to use. Currently supports "consensus" (default).

    Returns:
        List of dicts, each containing:
            - position: Splice site position (0-based, at the intron boundary)
            - type: "donor" or "acceptor"
            - score: Confidence score (0.0-1.0)
            - sequence_context: Surrounding sequence context
            - dinucleotide: The splice site dinucleotide (GT or AG)
            - strand: "+" (forward strand only)

    Raises:
        ValueError: If sequence is empty or model is unsupported.
    """
    if not sequence:
        raise ValueError("Sequence must not be empty")
    if model not in ("consensus",):
        raise ValueError(f"Unsupported model: {model}. Supported: 'consensus'")

    seq_upper = sequence.upper()
    results: list[dict] = []

    # Donor site consensus scoring weights (positions -3 to +6 relative to GT)
    # Based on Shapiro and Senapathy (1987) consensus
    donor_consensus = {
        -3: {"C": 0.35, "A": 0.35, "G": 0.15, "T": 0.15},
        -2: {"A": 0.6, "G": 0.1, "C": 0.1, "T": 0.2},
        -1: {"G": 0.8, "A": 0.1, "C": 0.05, "T": 0.05},
        # Position 0,1 = GT (required)
        2: {"A": 0.6, "G": 0.2, "T": 0.1, "C": 0.1},
        3: {"A": 0.7, "G": 0.1, "T": 0.1, "C": 0.1},
        4: {"G": 0.5, "T": 0.2, "A": 0.2, "C": 0.1},
        5: {"T": 0.5, "A": 0.2, "G": 0.2, "C": 0.1},
    }

    # Acceptor site consensus scoring (positions relative to AG)
    # Polypyrimidine tract typically 15-40 nt upstream
    acceptor_pyrimidine_weight = 0.6

    # Scan for donor sites (GT dinucleotide)
    for i in range(len(seq_upper) - 1):
        if seq_upper[i : i + 2] == "GT":
            score = _score_donor_site(seq_upper, i, donor_consensus)
            if score > 0.3:  # Minimum threshold
                context_start = max(0, i - 5)
                context_end = min(len(seq_upper), i + 8)
                results.append(
                    {
                        "position": i,
                        "type": "donor",
                        "score": round(score, 4),
                        "sequence_context": seq_upper[context_start:context_end],
                        "dinucleotide": "GT",
                        "strand": "+",
                    }
                )

    # Scan for acceptor sites (AG dinucleotide)
    for i in range(1, len(seq_upper)):
        if seq_upper[i - 1 : i + 1] == "AG":
            score = _score_acceptor_site(seq_upper, i - 1, acceptor_pyrimidine_weight)
            if score > 0.3:  # Minimum threshold
                context_start = max(0, i - 20)
                context_end = min(len(seq_upper), i + 5)
                results.append(
                    {
                        "position": i - 1,
                        "type": "acceptor",
                        "score": round(score, 4),
                        "sequence_context": seq_upper[context_start:context_end],
                        "dinucleotide": "AG",
                        "strand": "+",
                    }
                )

    # Sort by position then score
    results.sort(key=lambda x: (x["position"], -x["score"]))
    logger.debug(
        f"Found {len(results)} potential splice sites "
        f"({sum(1 for r in results if r['type'] == 'donor')} donors, "
        f"{sum(1 for r in results if r['type'] == 'acceptor')} acceptors)"
    )
    return results


def _score_donor_site(
    sequence: str,
    gt_position: int,
    consensus: Dict[int, Dict[str, float]],
) -> float:
    """Score a GT donor splice site against the consensus.

    Args:
        sequence: Full DNA sequence.
        gt_position: Position of the G in the GT dinucleotide.
        consensus: Position weight matrix for donor site.

    Returns:
        Score between 0.0 and 1.0.
    """
    score = 0.0
    max_score = 0.0

    for offset, weights in consensus.items():
        pos = gt_position + offset
        if 0 <= pos < len(sequence):
            base = sequence[pos]
            score += weights.get(base, 0.0)
            max_score += max(weights.values())
        else:
            max_score += max(weights.values())

    return score / max_score if max_score > 0 else 0.0


def _score_acceptor_site(
    sequence: str,
    ag_position: int,
    pyrimidine_weight: float,
) -> float:
    """Score an AG acceptor splice site based on polypyrimidine tract.

    The score is based on the fraction of pyrimidines (C, T) in the
    upstream region (typically 5-30 nt before the AG).

    Args:
        sequence: Full DNA sequence.
        ag_position: Position of the A in the AG dinucleotide.
        pyrimidine_weight: Weight for pyrimidine content scoring.

    Returns:
        Score between 0.0 and 1.0.
    """
    # Check upstream polypyrimidine tract (20 nt window before AG)
    tract_start = max(0, ag_position - 25)
    tract_end = max(0, ag_position - 3)  # Leave small gap before AG

    if tract_end <= tract_start:
        return 0.0

    tract = sequence[tract_start:tract_end]
    if not tract:
        return 0.0

    pyrimidine_count = sum(1 for b in tract if b in "CT")
    pyrimidine_fraction = pyrimidine_count / len(tract)

    # Score based on pyrimidine richness
    base_score = pyrimidine_fraction * pyrimidine_weight

    # Bonus for branch point consensus (YNYURAY ~20-50 nt upstream)
    branch_start = max(0, ag_position - 50)
    branch_end = max(0, ag_position - 18)
    if branch_end > branch_start:
        branch_region = sequence[branch_start:branch_end]
        # Look for branch point A (the adenosine that forms the lariat)
        if "A" in branch_region:
            base_score += 0.1

    # Bonus for the base immediately before AG being C or T
    if ag_position > 0 and sequence[ag_position - 1] in "CT":
        base_score += 0.1

    return min(1.0, base_score)
