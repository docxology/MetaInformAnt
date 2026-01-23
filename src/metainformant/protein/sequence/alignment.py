"""Protein sequence alignment utilities.

This module provides functions for aligning protein sequences using
various algorithms and scoring matrices.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)


def global_align(seq1: str, seq2: str, match: int = 1, mismatch: int = -1, gap: int = -2) -> Dict[str, Any]:
    """Perform global alignment of two protein sequences using Needleman-Wunsch algorithm.

    Args:
        seq1: First protein sequence
        seq2: Second protein sequence
        match: Match score
        mismatch: Mismatch penalty
        gap: Gap penalty

    Returns:
        Dictionary with alignment results
    """
    logger.debug(f"Performing global alignment of sequences (len: {len(seq1)}, {len(seq2)})")

    # Initialize scoring matrix
    m, n = len(seq1), len(seq2)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize first row and column
    for i in range(1, m + 1):
        score_matrix[i][0] = i * gap
    for j in range(1, n + 1):
        score_matrix[0][j] = j * gap

    # Fill scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diagonal = score_matrix[i - 1][j - 1] + match_score
            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(diagonal, up, left)

    # Traceback to get aligned sequences
    aligned_seq1, aligned_seq2 = _traceback_global(seq1, seq2, score_matrix, match, mismatch, gap)

    # Calculate identity
    identity = calculate_alignment_identity({"aligned_seq1": aligned_seq1, "aligned_seq2": aligned_seq2})

    result = {
        "aligned_seq1": aligned_seq1,
        "aligned_seq2": aligned_seq2,
        "score": score_matrix[m][n],
        "identity": identity,
        "method": "global",
        "match_score": match,
        "mismatch_penalty": mismatch,
        "gap_penalty": gap,
    }

    return result


def _traceback_global(
    seq1: str, seq2: str, score_matrix: List[List[int]], match: int, mismatch: int, gap: int
) -> Tuple[str, str]:
    """Traceback through scoring matrix for global alignment."""
    aligned_seq1 = []
    aligned_seq2 = []

    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            match_score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            if score_matrix[i][j] == score_matrix[i - 1][j - 1] + match_score:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
                continue

        if i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + gap:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        elif j > 0 and score_matrix[i][j] == score_matrix[i][j - 1] + gap:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Reverse the alignments
    aligned_seq1.reverse()
    aligned_seq2.reverse()

    return "".join(aligned_seq1), "".join(aligned_seq2)


def local_align(seq1: str, seq2: str, match: int = 1, mismatch: int = -1, gap: int = -2) -> Dict[str, Any]:
    """Perform local alignment of two protein sequences using Smith-Waterman algorithm.

    Args:
        seq1: First protein sequence
        seq2: Second protein sequence
        match: Match score
        mismatch: Mismatch penalty
        gap: Gap penalty

    Returns:
        Dictionary with alignment results
    """
    logger.debug(f"Performing local alignment of sequences (len: {len(seq1)}, {len(seq2)})")

    # Initialize scoring matrix
    m, n = len(seq1), len(seq2)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]

    max_score = 0
    max_i, max_j = 0, 0

    # Fill scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diagonal = score_matrix[i - 1][j - 1] + match_score
            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(0, diagonal, up, left)  # Local alignment allows 0

            # Track maximum score position
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_i, max_j = i, j

    # Traceback from maximum score position
    aligned_seq1, aligned_seq2 = _traceback_local(seq1, seq2, score_matrix, max_i, max_j, match, mismatch, gap)

    # Calculate identity
    identity = calculate_alignment_identity({"aligned_seq1": aligned_seq1, "aligned_seq2": aligned_seq2})

    result = {
        "aligned_seq1": aligned_seq1,
        "aligned_seq2": aligned_seq2,
        "score": max_score,
        "identity": identity,
        "method": "local",
        "match_score": match,
        "mismatch_penalty": mismatch,
        "gap_penalty": gap,
    }

    return result


def _traceback_local(
    seq1: str, seq2: str, score_matrix: List[List[int]], start_i: int, start_j: int, match: int, mismatch: int, gap: int
) -> Tuple[str, str]:
    """Traceback through scoring matrix for local alignment."""
    aligned_seq1 = []
    aligned_seq2 = []

    i, j = start_i, start_j

    while score_matrix[i][j] > 0:
        if i > 0 and j > 0:
            match_score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            if score_matrix[i][j] == score_matrix[i - 1][j - 1] + match_score:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
                continue

        if i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + gap:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        elif j > 0 and score_matrix[i][j] == score_matrix[i][j - 1] + gap:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1
        else:
            break  # Should not happen in local alignment

    # Reverse the alignments
    aligned_seq1.reverse()
    aligned_seq2.reverse()

    return "".join(aligned_seq1), "".join(aligned_seq2)


def calculate_alignment_identity(alignment: Dict[str, Any]) -> float:
    """Calculate identity percentage from alignment result.

    Args:
        alignment: Alignment result dictionary

    Returns:
        Identity percentage (0.0 to 1.0)
    """
    aligned1 = alignment.get("aligned_seq1", "")
    aligned2 = alignment.get("aligned_seq2", "")

    if not aligned1 or not aligned2 or len(aligned1) != len(aligned2):
        return 0.0

    matches = 0
    total = 0

    for a, b in zip(aligned1, aligned2):
        if a != "-" and b != "-":
            total += 1
            if a == b:
                matches += 1

    return matches / total if total > 0 else 0.0


def alignment_statistics(alignment: Dict[str, Any]) -> Dict[str, float]:
    """Calculate comprehensive statistics for an alignment.

    Args:
        alignment: Alignment result dictionary

    Returns:
        Dictionary with alignment statistics
    """
    aligned1 = alignment.get("aligned_seq1", "")
    aligned2 = alignment.get("aligned_seq2", "")

    if not aligned1 or not aligned2:
        return {}

    identity = calculate_alignment_identity(alignment)

    # Count gaps
    gaps1 = aligned1.count("-")
    gaps2 = aligned2.count("-")

    # Calculate similarity (considering conservative substitutions)
    similar = 0
    total_positions = 0

    for a, b in zip(aligned1, aligned2):
        if a != "-" and b != "-":
            total_positions += 1
            if a == b:
                similar += 1
            # Could add similarity matrix logic here

    similarity = similar / total_positions if total_positions > 0 else 0.0

    return {
        "identity": identity,
        "similarity": similarity,
        "gaps_seq1": gaps1,
        "gaps_seq2": gaps2,
        "aligned_length": len(aligned1),
        "raw_length_seq1": len(alignment.get("seq1", "")),
        "raw_length_seq2": len(alignment.get("seq2", "")),
    }
