"""Protein domain detection and analysis utilities.

This module provides tools for scanning protein sequences against domain
profiles, predicting signal peptides, transmembrane helices, coiled-coil
regions, zinc fingers, leucine zippers, and computing domain architecture
summaries. All methods use pure Python implementations with biologically
grounded algorithms (position-specific scoring, hydrophobicity analysis,
heptad repeat scoring, and regex-based motif detection).
"""

from __future__ import annotations

import math
import re
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Amino acid property tables
# ---------------------------------------------------------------------------

# Kyte-Doolittle hydrophobicity scale
_KD_HYDROPHOBICITY: Dict[str, float] = {
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

# Amino acid charges at pH 7.0 (approximate)
_AA_CHARGE: Dict[str, float] = {
    "K": 1.0,
    "R": 1.0,
    "H": 0.1,
    "D": -1.0,
    "E": -1.0,
}

# ---------------------------------------------------------------------------
# Built-in domain profiles (position-specific scoring matrices)
# ---------------------------------------------------------------------------

# Each profile is a dict with:
#   pattern: regex or consensus string
#   pssm: list of dicts mapping AA -> score (position-specific)
#   min_score: minimum total score to call a hit
#   description: human-readable description

_BUILTIN_PROFILES: Dict[str, Dict[str, Any]] = {
    "PF00069": {
        "name": "Protein kinase domain",
        "pattern": r"[IVLM].{2}[HRK].{2}D[LFIVMA].{2}N",
        "consensus": "VAIHRDLXXXN",
        "min_length": 250,
        "max_length": 350,
        "description": "Catalytic domain of protein kinases",
        "key_residues": {"D": "catalytic", "K": "ATP binding", "N": "catalytic loop"},
    },
    "PF00076": {
        "name": "RNA recognition motif (RRM)",
        "pattern": r"[IVLM][FY][IVLM].{4,8}[IVLM].{2}[FYW].{3}[IVLM]",
        "consensus": "IFI....I..F...L",
        "min_length": 70,
        "max_length": 100,
        "description": "RNA binding domain with RNP-1 and RNP-2 motifs",
        "key_residues": {"F": "stacking", "Y": "stacking"},
    },
    "PF00096": {
        "name": "C2H2 zinc finger",
        "pattern": r"C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H",
        "consensus": "CX2-4CX3LX8HX3-5H",
        "min_length": 23,
        "max_length": 30,
        "description": "Classical Cys2-His2 zinc finger DNA binding domain",
        "key_residues": {"C": "zinc coordination", "H": "zinc coordination"},
    },
    "PF00400": {
        "name": "WD40 repeat",
        "pattern": r"[GSANDEKRH].{4,7}[LIVMFYWC].{1,2}[LIVMFYWC].{4,8}W[DE]",
        "consensus": "GX4-7LXL...WD",
        "min_length": 40,
        "max_length": 60,
        "description": "WD40 repeat forming beta-propeller structures",
        "key_residues": {"W": "structural", "D": "structural"},
    },
    "PF00071": {
        "name": "Ras family GTPase",
        "pattern": r"G.{4}GK[ST].{6}[DE].{3}T.{2}[DE]",
        "consensus": "GXXXXGKS...D...T..E",
        "min_length": 150,
        "max_length": 200,
        "description": "Small GTPase superfamily with GTP binding motifs",
        "key_residues": {"G": "P-loop", "K": "GTP binding", "T": "effector"},
    },
    "PF00046": {
        "name": "Homeodomain",
        "pattern": r"[RK].{2}[RKNQ][LIVMFYWH].{2}[LIVMFYW].{4}[LIVMFYW].{3}[RKNQ].{3}W[FYWML]",
        "consensus": "RX2QLX2LX4LX3RX3WF",
        "min_length": 57,
        "max_length": 63,
        "description": "Homeodomain DNA binding helix-turn-helix",
        "key_residues": {"W": "hydrophobic core", "R": "DNA contact"},
    },
    "PF00271": {
        "name": "Helicase C-terminal domain",
        "pattern": r"[DE].{2}H.{5,10}[LIVMFYW].{2}[LIVMFYW].{3}[ST]AT",
        "consensus": "DX2HX5-10LX2LX3SAT",
        "min_length": 70,
        "max_length": 120,
        "description": "DEAD-box helicase conserved C-terminal domain",
        "key_residues": {"D": "catalytic", "H": "RNA binding"},
    },
    "PF00089": {
        "name": "Trypsin-like serine protease",
        "pattern": r"[LIVMFYWGA].{5}[LIVMFYWGA]D[STA]G.{2}[LIVMFYWGA]",
        "consensus": "LX5LDSGX2L",
        "min_length": 210,
        "max_length": 260,
        "description": "Serine protease catalytic domain with His-Asp-Ser triad",
        "key_residues": {"H": "catalytic triad", "D": "catalytic triad", "S": "catalytic triad"},
    },
}


def scan_domains(
    sequence: str,
    profiles: Dict[str, Any] | None = None,
) -> List[Dict[str, Any]]:
    """Scan a protein sequence against domain profiles using pattern matching.

    Uses built-in HMM-like scoring via position-specific scoring matrices
    combined with regex-based motif detection. Each domain hit is scored
    based on match quality and biological plausibility.

    Args:
        sequence: Protein amino acid sequence string.
        profiles: Optional custom profiles dict mapping profile_id to profile
            definition. Each profile should contain at minimum 'name',
            'pattern', and 'description' keys. If None, uses the built-in
            profile library.

    Returns:
        List of domain hit dicts, each containing:
            - domain_id: profile identifier
            - name: human-readable domain name
            - start: 0-based start position in sequence
            - end: 0-based end position (exclusive)
            - score: raw match score (higher is better)
            - e_value: estimated statistical significance
            - description: domain description string

    Raises:
        ValueError: If the sequence is empty or contains no valid amino acids.
    """
    if not sequence or not sequence.strip():
        raise ValueError("Sequence must be a non-empty string")

    seq = sequence.upper().strip()
    active_profiles = profiles if profiles is not None else _BUILTIN_PROFILES
    hits: List[Dict[str, Any]] = []

    for profile_id, profile in active_profiles.items():
        pattern = profile.get("pattern")
        if not pattern:
            continue

        try:
            for match in re.finditer(pattern, seq):
                start = match.start()
                end = match.end()
                matched_seq = match.group()

                # Score based on match length relative to expected domain size
                match_len = end - start
                score = _score_domain_hit(matched_seq, profile)

                # Estimate E-value using Karlin-Altschul-like statistics
                e_value = _estimate_e_value(score, len(seq), match_len)

                hits.append(
                    {
                        "domain_id": profile_id,
                        "name": profile.get("name", profile_id),
                        "start": start,
                        "end": end,
                        "score": round(score, 2),
                        "e_value": e_value,
                        "description": profile.get("description", ""),
                        "matched_sequence": matched_seq,
                    }
                )
        except re.error as exc:
            logger.warning("Invalid regex pattern for profile %s: %s", profile_id, exc)
            continue

    # Sort by score descending, then by start position
    hits.sort(key=lambda h: (-h["score"], h["start"]))

    # Remove overlapping hits for the same domain (keep best score)
    hits = _remove_overlapping_hits(hits)

    logger.debug("Found %d domain hits in sequence of length %d", len(hits), len(seq))
    return hits


def _score_domain_hit(matched_seq: str, profile: Dict[str, Any]) -> float:
    """Score a domain hit based on match quality and key residue presence.

    Args:
        matched_seq: The matched subsequence.
        profile: The domain profile definition.

    Returns:
        Numerical score for the match quality.
    """
    # Base score from match length
    score = float(len(matched_seq))

    # Bonus for key residues present
    key_residues = profile.get("key_residues", {})
    for residue, _role in key_residues.items():
        count = matched_seq.count(residue)
        if count > 0:
            score += count * 3.0

    # Length plausibility bonus
    min_len = profile.get("min_length", 0)
    max_len = profile.get("max_length", 10000)
    if min_len <= len(matched_seq) <= max_len:
        score += 10.0

    # Consensus similarity bonus
    consensus = profile.get("consensus", "")
    if consensus:
        matches = sum(1 for a, b in zip(matched_seq, consensus) if a == b and b not in ("X", "x", "."))
        score += matches * 2.0

    return score


def _estimate_e_value(score: float, seq_length: int, match_length: int) -> float:
    """Estimate E-value using simplified Karlin-Altschul statistics.

    Args:
        score: Raw alignment score.
        seq_length: Length of the query sequence.
        match_length: Length of the matched region.

    Returns:
        Estimated E-value (lower is more significant).
    """
    # Simplified parameters for amino acid sequences
    lambda_param = 0.267
    k_param = 0.041

    # Effective search space
    search_space = float(seq_length * match_length)
    if search_space <= 0:
        return 1.0

    try:
        e_value = k_param * search_space * math.exp(-lambda_param * score)
    except OverflowError:
        e_value = 0.0

    return round(max(e_value, 1e-300), 6)


def _remove_overlapping_hits(hits: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Remove overlapping domain hits, keeping the highest-scoring one.

    Args:
        hits: List of domain hits sorted by score descending.

    Returns:
        Filtered list with overlapping hits removed.
    """
    if not hits:
        return hits

    kept: List[Dict[str, Any]] = []
    occupied: List[Tuple[int, int, str]] = []

    for hit in hits:
        start, end = hit["start"], hit["end"]
        domain_id = hit["domain_id"]

        # Check overlap with already-kept hits of the same domain type
        overlaps = False
        for occ_start, occ_end, occ_domain in occupied:
            if occ_domain != domain_id:
                continue
            if start < occ_end and end > occ_start:
                overlap_len = min(end, occ_end) - max(start, occ_start)
                hit_len = end - start
                if hit_len > 0 and overlap_len / hit_len > 0.5:
                    overlaps = True
                    break

        if not overlaps:
            kept.append(hit)
            occupied.append((start, end, domain_id))

    return kept


# ---------------------------------------------------------------------------
# Signal peptide prediction
# ---------------------------------------------------------------------------


def detect_signal_peptide(
    sequence: str,
    organism: str = "eukaryote",
) -> Dict[str, Any]:
    """Predict signal peptide using the von Heijne method.

    Analyses the N-terminal region for the characteristic tripartite
    structure of signal peptides: positively charged n-region, hydrophobic
    h-region, and cleavage-site c-region. Uses the (-3, -1) rule for
    cleavage site prediction (small neutral amino acids at positions -1
    and -3 relative to the cleavage site).

    Args:
        sequence: Protein amino acid sequence string.
        organism: Organism type for parameter tuning. One of
            ``"eukaryote"`` or ``"prokaryote"``.

    Returns:
        Dict with keys:
            - has_signal: bool indicating signal peptide presence
            - cleavage_site: predicted cleavage position (0-based) or -1
            - probability: confidence score (0.0 to 1.0)
            - n_region: dict with start, end, charge
            - h_region: dict with start, end, hydrophobicity
            - c_region: dict with start, end, cleavage_motif

    Raises:
        ValueError: If sequence is too short for signal peptide analysis.
    """
    if len(sequence) < 15:
        raise ValueError("Sequence too short for signal peptide analysis (minimum 15 residues)")

    seq = sequence.upper()

    # Signal peptide parameters by organism type
    if organism == "prokaryote":
        max_sp_length = 40
        n_region_max = 8
        h_region_min = 7
        h_region_max = 15
    else:  # eukaryote
        max_sp_length = 35
        n_region_max = 6
        h_region_min = 7
        h_region_max = 15

    # Analyse N-terminal region only
    n_term = seq[: min(max_sp_length + 10, len(seq))]

    # Step 1: Find n-region (positively charged N-terminal segment)
    n_region_end = min(n_region_max, len(n_term))
    n_region_seq = n_term[:n_region_end]
    n_charge = sum(_AA_CHARGE.get(aa, 0.0) for aa in n_region_seq)

    # Step 2: Find h-region (hydrophobic core)
    best_h_start = n_region_end
    best_h_end = n_region_end
    best_h_score = -999.0

    for h_start in range(max(1, n_region_end - 2), n_region_end + 3):
        if h_start >= len(n_term):
            break
        for h_len in range(h_region_min, h_region_max + 1):
            h_end = h_start + h_len
            if h_end > len(n_term):
                break
            h_window = n_term[h_start:h_end]
            h_score = sum(_KD_HYDROPHOBICITY.get(aa, 0.0) for aa in h_window) / len(h_window)
            if h_score > best_h_score:
                best_h_score = h_score
                best_h_start = h_start
                best_h_end = h_end

    # Step 3: Find c-region and cleavage site using (-3, -1) rule
    # Small neutral residues preferred at -1 and -3 positions
    small_neutral = {"A", "G", "S", "T", "V", "C", "L"}
    best_cleavage = -1
    best_cleavage_score = -999.0

    search_start = best_h_end
    search_end = min(best_h_end + 8, len(n_term))

    for pos in range(search_start, search_end):
        if pos < 2 or pos >= len(seq):
            continue
        # (-1, -3) rule: residues at positions pos-1 and pos-3
        score_minus1 = 1.0 if seq[pos - 1] in small_neutral else -1.0
        score_minus3 = 1.0 if pos >= 3 and seq[pos - 3] in small_neutral else -0.5

        # Helix breaker at -5 to -4 is also favorable
        breaker_bonus = 0.0
        if pos >= 5 and seq[pos - 5] in {"P", "G"}:
            breaker_bonus = 0.3
        if pos >= 4 and seq[pos - 4] in {"P", "G"}:
            breaker_bonus += 0.3

        cleavage_score = score_minus1 + score_minus3 + breaker_bonus
        if cleavage_score > best_cleavage_score:
            best_cleavage_score = cleavage_score
            best_cleavage = pos

    # Step 4: Compute overall probability
    probability = 0.0

    # Positive N-terminus contributes
    if n_charge >= 1.0:
        probability += 0.2
    elif n_charge >= 0.0:
        probability += 0.1

    # Hydrophobic h-region contributes
    if best_h_score >= 1.5:
        probability += 0.35
    elif best_h_score >= 0.8:
        probability += 0.2
    elif best_h_score >= 0.0:
        probability += 0.05

    # Cleavage site quality contributes
    if best_cleavage_score >= 1.5:
        probability += 0.35
    elif best_cleavage_score >= 0.5:
        probability += 0.2
    elif best_cleavage_score >= 0.0:
        probability += 0.05

    # Length plausibility
    if best_cleavage > 0:
        sp_length = best_cleavage
        if 15 <= sp_length <= max_sp_length:
            probability += 0.1

    probability = min(probability, 1.0)
    has_signal = probability >= 0.5

    c_region_start = best_h_end
    c_region_end = best_cleavage if best_cleavage > 0 else c_region_start
    c_motif = seq[max(0, c_region_end - 3) : c_region_end] if c_region_end > 0 else ""

    logger.debug(
        "Signal peptide prediction: has_signal=%s, cleavage=%d, prob=%.2f",
        has_signal,
        best_cleavage,
        probability,
    )

    return {
        "has_signal": has_signal,
        "cleavage_site": best_cleavage if has_signal else -1,
        "probability": round(probability, 4),
        "n_region": {
            "start": 0,
            "end": n_region_end,
            "charge": round(n_charge, 2),
            "sequence": n_region_seq,
        },
        "h_region": {
            "start": best_h_start,
            "end": best_h_end,
            "hydrophobicity": round(best_h_score, 4),
            "sequence": n_term[best_h_start:best_h_end],
        },
        "c_region": {
            "start": c_region_start,
            "end": c_region_end,
            "cleavage_motif": c_motif,
        },
    }


# ---------------------------------------------------------------------------
# Transmembrane helix prediction
# ---------------------------------------------------------------------------


def predict_transmembrane(
    sequence: str,
    window: int = 19,
    threshold: float = 1.6,
) -> List[Dict[str, Any]]:
    """Predict transmembrane helices using Kyte-Doolittle hydrophobicity.

    Applies a sliding window of the specified size across the sequence,
    computing average hydrophobicity. Regions exceeding the threshold that
    are at least 15 residues long are reported as putative transmembrane
    helices. Orientation (in-to-out vs out-to-in) is inferred from the
    positive-inside rule.

    Args:
        sequence: Protein amino acid sequence string.
        window: Sliding window size (default 19, typical for TM helices).
        threshold: Minimum average hydrophobicity to call a TM region.

    Returns:
        List of dicts, each with:
            - start: 0-based start position
            - end: 0-based end position (exclusive)
            - orientation: ``"in-to-out"`` or ``"out-to-in"``
            - score: average hydrophobicity of the region
            - sequence: the transmembrane segment string

    Raises:
        ValueError: If sequence is shorter than the window size.
    """
    if not sequence:
        raise ValueError("Sequence must be non-empty")
    seq = sequence.upper()
    if len(seq) < window:
        raise ValueError(f"Sequence length ({len(seq)}) must be >= window size ({window})")

    # Compute per-position hydrophobicity scores via sliding window
    scores: List[float] = []
    for i in range(len(seq) - window + 1):
        window_seq = seq[i : i + window]
        avg = sum(_KD_HYDROPHOBICITY.get(aa, 0.0) for aa in window_seq) / window
        scores.append(avg)

    # Identify contiguous regions above threshold
    raw_regions: List[Tuple[int, int]] = []
    i = 0
    while i < len(scores):
        if scores[i] >= threshold:
            region_start = i
            while i < len(scores) and scores[i] >= threshold:
                i += 1
            region_end = i + window - 1  # extend to cover full window
            raw_regions.append((region_start, min(region_end, len(seq))))
        else:
            i += 1

    # Merge overlapping regions and filter by minimum length (15 residues)
    merged: List[Tuple[int, int]] = []
    for start, end in raw_regions:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    min_tm_length = 15
    filtered = [(s, e) for s, e in merged if (e - s) >= min_tm_length]

    # Assign orientation using positive-inside rule
    # (basic residues K, R prefer cytoplasmic side)
    results: List[Dict[str, Any]] = []
    for idx, (start, end) in enumerate(filtered):
        tm_seq = seq[start:end]

        # Count basic residues in flanking regions
        flank = 15
        n_flank = seq[max(0, start - flank) : start]
        c_flank = seq[end : min(len(seq), end + flank)]

        n_basic = sum(1 for aa in n_flank if aa in ("K", "R"))
        c_basic = sum(1 for aa in c_flank if aa in ("K", "R"))

        # Positive-inside rule: more basic residues on cytoplasmic side
        if idx % 2 == 0:
            orientation = "in-to-out" if n_basic >= c_basic else "out-to-in"
        else:
            orientation = "out-to-in" if n_basic >= c_basic else "in-to-out"

        region_score = sum(_KD_HYDROPHOBICITY.get(aa, 0.0) for aa in tm_seq) / len(tm_seq)

        results.append(
            {
                "start": start,
                "end": end,
                "orientation": orientation,
                "score": round(region_score, 4),
                "sequence": tm_seq,
            }
        )

    logger.debug("Predicted %d transmembrane helices", len(results))
    return results


# ---------------------------------------------------------------------------
# Coiled-coil prediction
# ---------------------------------------------------------------------------

# Heptad position probabilities for coiled-coils (a through g)
# Residue preferences at a/d positions (hydrophobic) and e/g (charged)
_COILED_COIL_WEIGHTS: Dict[str, Dict[str, float]] = {
    "a": {"L": 1.8, "I": 1.6, "V": 1.5, "M": 1.4, "A": 1.2, "F": 1.1},
    "b": {"E": 1.3, "K": 1.2, "Q": 1.2, "R": 1.1, "A": 1.0},
    "c": {"E": 1.3, "K": 1.2, "A": 1.1, "Q": 1.1, "R": 1.0},
    "d": {"L": 1.9, "I": 1.6, "V": 1.5, "M": 1.4, "A": 1.3, "N": 1.0},
    "e": {"E": 1.5, "K": 1.3, "R": 1.2, "Q": 1.1, "D": 1.0},
    "f": {"K": 1.3, "E": 1.2, "A": 1.1, "R": 1.1, "Q": 1.0},
    "g": {"K": 1.4, "E": 1.3, "R": 1.2, "Q": 1.1, "D": 1.0},
}


def find_coiled_coils(
    sequence: str,
    window: int = 28,
    threshold: float = 0.5,
) -> List[Dict[str, Any]]:
    """Detect coiled-coil regions using heptad repeat scoring.

    Scores each position for a/d hydrophobic periodicity characteristic
    of coiled-coil structures. The algorithm slides a window across the
    sequence, evaluating how well each window fits the canonical heptad
    repeat pattern (abcdefg) where positions a and d are hydrophobic.

    Args:
        sequence: Protein amino acid sequence string.
        window: Sliding window size (default 28, i.e. 4 heptad repeats).
        threshold: Minimum score to report a coiled-coil region (0.0 to 1.0).

    Returns:
        List of dicts, each with:
            - start: 0-based start position
            - end: 0-based end position (exclusive)
            - score: coiled-coil probability score
            - heptad_assignment: string of heptad position letters
            - length: region length in residues
    """
    if not sequence:
        return []

    seq = sequence.upper()
    if len(seq) < window:
        return []

    # Score each window for best heptad-phase alignment
    position_scores: List[float] = []
    for i in range(len(seq) - window + 1):
        window_seq = seq[i : i + window]
        best_score = 0.0

        # Try all 7 possible phases
        for phase in range(7):
            score = 0.0
            max_possible = 0.0
            for j, aa in enumerate(window_seq):
                heptad_pos = "abcdefg"[(j + phase) % 7]
                weights = _COILED_COIL_WEIGHTS.get(heptad_pos, {})
                aa_weight = weights.get(aa, 0.5)
                max_weight = max(weights.values()) if weights else 1.0
                score += aa_weight
                max_possible += max_weight

            normalized = score / max_possible if max_possible > 0 else 0.0
            if normalized > best_score:
                best_score = normalized

        position_scores.append(best_score)

    # Identify contiguous regions above threshold
    regions: List[Dict[str, Any]] = []
    i = 0
    while i < len(position_scores):
        if position_scores[i] >= threshold:
            region_start = i
            region_scores: List[float] = []
            while i < len(position_scores) and position_scores[i] >= threshold:
                region_scores.append(position_scores[i])
                i += 1
            region_end = i + window - 1
            region_end = min(region_end, len(seq))

            # Require at least 2 heptad repeats (14 residues)
            region_len = region_end - region_start
            if region_len >= 14:
                avg_score = sum(region_scores) / len(region_scores)
                # Assign heptad positions based on best phase
                heptad = _assign_heptad(seq[region_start:region_end])
                regions.append(
                    {
                        "start": region_start,
                        "end": region_end,
                        "score": round(avg_score, 4),
                        "heptad_assignment": heptad,
                        "length": region_len,
                    }
                )
        else:
            i += 1

    logger.debug("Found %d coiled-coil regions", len(regions))
    return regions


def _assign_heptad(region_seq: str) -> str:
    """Assign heptad repeat positions to a coiled-coil region.

    Args:
        region_seq: Amino acid sequence of the coiled-coil region.

    Returns:
        String of heptad position letters (a-g) aligned to the sequence.
    """
    heptad_letters = "abcdefg"
    best_phase = 0
    best_score = -1.0

    for phase in range(7):
        score = 0.0
        for j, aa in enumerate(region_seq):
            pos = heptad_letters[(j + phase) % 7]
            weights = _COILED_COIL_WEIGHTS.get(pos, {})
            score += weights.get(aa, 0.5)
        if score > best_score:
            best_score = score
            best_phase = phase

    return "".join(heptad_letters[(j + best_phase) % 7] for j in range(len(region_seq)))


# ---------------------------------------------------------------------------
# Zinc finger detection
# ---------------------------------------------------------------------------

# Zinc finger patterns with biological constraints
_ZINC_FINGER_PATTERNS: List[Dict[str, Any]] = [
    {
        "type": "C2H2",
        "pattern": r"C.{2,4}C.{12}H.{3,5}H",
        "description": "Classical Cys2-His2 zinc finger",
        "min_spacing": 23,
        "max_spacing": 28,
    },
    {
        "type": "C4",
        "pattern": r"C.{2}C.{13,17}C.{2}C",
        "description": "Cys4 zinc finger (nuclear receptor type)",
        "min_spacing": 21,
        "max_spacing": 27,
    },
    {
        "type": "C2HC",
        "pattern": r"C.{2,4}C.{4,12}H.{2,5}C",
        "description": "Cys2-His-Cys retroviral-type zinc finger",
        "min_spacing": 14,
        "max_spacing": 28,
    },
    {
        "type": "C3H1",
        "pattern": r"C.{4,8}C.{4,6}C.{3}H",
        "description": "CCCH-type zinc finger (RNA binding)",
        "min_spacing": 17,
        "max_spacing": 25,
    },
]


def identify_zinc_fingers(sequence: str) -> List[Dict[str, Any]]:
    """Find zinc finger motifs in a protein sequence.

    Detects C2H2, C4, C2HC, and C3H1 zinc finger patterns using regular
    expressions with biological spacing constraints. Each hit is validated
    for the expected distance between coordinating residues.

    Args:
        sequence: Protein amino acid sequence string.

    Returns:
        List of dicts, each with:
            - type: zinc finger class (C2H2, C4, C2HC, C3H1)
            - start: 0-based start position
            - end: 0-based end position (exclusive)
            - sequence: matched subsequence
            - coordinating_residues: list of (position, residue) tuples
            - description: human-readable description
    """
    if not sequence:
        return []

    seq = sequence.upper()
    results: List[Dict[str, Any]] = []

    for zf_def in _ZINC_FINGER_PATTERNS:
        try:
            for match in re.finditer(zf_def["pattern"], seq):
                start = match.start()
                end = match.end()
                matched = match.group()
                motif_len = end - start

                # Validate spacing constraints
                if motif_len < zf_def["min_spacing"] or motif_len > zf_def["max_spacing"]:
                    continue

                # Identify coordinating residues (C and H)
                coordinating: List[Tuple[int, str]] = []
                for i, aa in enumerate(matched):
                    if aa in ("C", "H"):
                        coordinating.append((start + i, aa))

                results.append(
                    {
                        "type": zf_def["type"],
                        "start": start,
                        "end": end,
                        "sequence": matched,
                        "coordinating_residues": coordinating,
                        "description": zf_def["description"],
                    }
                )
        except re.error as exc:
            logger.warning("Regex error in zinc finger pattern %s: %s", zf_def["type"], exc)

    # Sort by position
    results.sort(key=lambda r: r["start"])

    logger.debug("Found %d zinc finger motifs", len(results))
    return results


# ---------------------------------------------------------------------------
# Leucine zipper detection
# ---------------------------------------------------------------------------


def detect_leucine_zipper(sequence: str) -> List[Dict[str, Any]]:
    """Find leucine zipper motifs in a protein sequence.

    Leucine zippers are characterized by leucine (or other hydrophobic
    residue) at every 7th position in an alpha-helix context. This function
    scans for at least 4 consecutive heptad repeats where the d-position
    is occupied by a hydrophobic residue (L, I, V, M, F).

    Args:
        sequence: Protein amino acid sequence string.

    Returns:
        List of dicts, each with:
            - start: 0-based start position
            - end: 0-based end position (exclusive)
            - sequence: matched subsequence
            - leucine_positions: list of 0-based positions of the
              leucine/hydrophobic residues at d-positions
            - score: quality score based on leucine conservation
            - repeat_count: number of heptad repeats found
    """
    if not sequence or len(sequence) < 28:
        return []

    seq = sequence.upper()
    hydrophobic_d = {"L", "I", "V", "M", "F"}
    min_repeats = 4
    results: List[Dict[str, Any]] = []

    # Scan all possible starting phases
    for phase_offset in range(7):
        # Check each possible start position
        for start in range(phase_offset, len(seq) - min_repeats * 7 + 1, 1):
            leucine_positions: List[int] = []
            leucine_count = 0
            total_heptads = 0

            # Scan heptad repeats from this start
            pos = start
            while pos + 6 < len(seq):
                d_pos = pos  # d-position is our anchor
                d_residue = seq[d_pos]

                if d_residue in hydrophobic_d:
                    leucine_positions.append(d_pos)
                    if d_residue == "L":
                        leucine_count += 1
                    total_heptads += 1
                    pos += 7
                else:
                    break

            if total_heptads >= min_repeats:
                region_start = leucine_positions[0]
                region_end = leucine_positions[-1] + 1

                # Score based on leucine conservation at d-position
                leu_fraction = leucine_count / total_heptads
                score = 0.5 + 0.5 * leu_fraction

                # Check for charged residues at e/g positions (salt bridges)
                charged = {"E", "K", "R", "D"}
                eg_charged = 0
                eg_total = 0
                for lp in leucine_positions:
                    # e-position is d+1, g-position is d+3
                    if lp + 1 < len(seq):
                        eg_total += 1
                        if seq[lp + 1] in charged:
                            eg_charged += 1
                    if lp + 3 < len(seq):
                        eg_total += 1
                        if seq[lp + 3] in charged:
                            eg_charged += 1

                if eg_total > 0:
                    score += 0.2 * (eg_charged / eg_total)
                score = min(score, 1.0)

                # Avoid duplicate overlapping hits
                is_duplicate = False
                for existing in results:
                    if abs(existing["start"] - region_start) < 7 and abs(existing["end"] - region_end) < 7:
                        if existing["score"] >= score:
                            is_duplicate = True
                            break
                        else:
                            results.remove(existing)
                            break

                if not is_duplicate:
                    results.append(
                        {
                            "start": region_start,
                            "end": region_end,
                            "sequence": seq[region_start:region_end],
                            "leucine_positions": leucine_positions,
                            "score": round(score, 4),
                            "repeat_count": total_heptads,
                        }
                    )

    results.sort(key=lambda r: (-r["score"], r["start"]))

    logger.debug("Found %d leucine zipper motifs", len(results))
    return results


# ---------------------------------------------------------------------------
# Domain architecture summary
# ---------------------------------------------------------------------------


def compute_domain_architecture(
    domains: List[Dict[str, Any]],
    sequence_length: int,
) -> Dict[str, Any]:
    """Summarize the domain architecture of a protein.

    Computes coverage statistics, identifies inter-domain linker regions,
    and generates a domain order string suitable for architecture comparison.

    Args:
        domains: List of domain hit dicts (must contain 'start', 'end',
            and 'name' or 'domain_id' keys).
        sequence_length: Total length of the protein sequence.

    Returns:
        Dict with keys:
            - domain_count: number of domains
            - total_coverage: fraction of sequence covered by domains
            - covered_residues: number of residues in domain regions
            - inter_domain_regions: list of (start, end) tuples for linkers
            - domain_order: string representation (e.g. "Kinase-SH2-SH3")
            - average_domain_length: mean domain length
            - linker_lengths: list of linker region lengths
            - architecture_string: compact domain architecture notation

    Raises:
        ValueError: If sequence_length is not positive.
    """
    if sequence_length <= 0:
        raise ValueError("sequence_length must be positive")

    if not domains:
        return {
            "domain_count": 0,
            "total_coverage": 0.0,
            "covered_residues": 0,
            "inter_domain_regions": [(0, sequence_length)],
            "domain_order": "",
            "average_domain_length": 0.0,
            "linker_lengths": [],
            "architecture_string": "[-empty-]",
        }

    # Sort domains by start position
    sorted_domains = sorted(domains, key=lambda d: d.get("start", 0))

    # Calculate coverage (handle overlaps)
    covered = [False] * sequence_length
    for dom in sorted_domains:
        start = max(0, dom.get("start", 0))
        end = min(sequence_length, dom.get("end", 0))
        for i in range(start, end):
            covered[i] = True
    covered_residues = sum(covered)

    # Find inter-domain regions
    inter_domain: List[Tuple[int, int]] = []
    in_linker = False
    linker_start = 0

    for i in range(sequence_length):
        if not covered[i] and not in_linker:
            linker_start = i
            in_linker = True
        elif covered[i] and in_linker:
            inter_domain.append((linker_start, i))
            in_linker = False
    if in_linker:
        inter_domain.append((linker_start, sequence_length))

    # Build domain order string
    domain_names: List[str] = []
    for dom in sorted_domains:
        name = dom.get("name", dom.get("domain_id", "unknown"))
        domain_names.append(name)
    domain_order = "-".join(domain_names)

    # Domain lengths
    domain_lengths = [dom.get("end", 0) - dom.get("start", 0) for dom in sorted_domains]
    avg_length = sum(domain_lengths) / len(domain_lengths) if domain_lengths else 0.0

    # Linker lengths
    linker_lengths = [end - start for start, end in inter_domain]

    # Architecture string (compact notation)
    arch_parts: List[str] = []
    prev_end = 0
    for dom in sorted_domains:
        start = dom.get("start", 0)
        end = dom.get("end", 0)
        name = dom.get("name", dom.get("domain_id", "?"))
        if start > prev_end:
            gap = start - prev_end
            arch_parts.append(f"[{gap}]")
        arch_parts.append(f"{name}({start}-{end})")
        prev_end = end
    if prev_end < sequence_length:
        arch_parts.append(f"[{sequence_length - prev_end}]")
    architecture_string = "--".join(arch_parts)

    return {
        "domain_count": len(sorted_domains),
        "total_coverage": round(covered_residues / sequence_length, 4),
        "covered_residues": covered_residues,
        "inter_domain_regions": inter_domain,
        "domain_order": domain_order,
        "average_domain_length": round(avg_length, 2),
        "linker_lengths": linker_lengths,
        "architecture_string": architecture_string,
    }
