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


# BLOSUM62 substitution matrix
BLOSUM62: Dict[str, Dict[str, int]] = {}
_blosum62_data = {
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
        "A": -1,
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
        "A": -2,
        "R": 0,
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
        "A": -2,
        "R": -2,
        "N": 1,
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
        "A": 0,
        "R": -3,
        "N": -3,
        "D": -3,
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
        "A": -1,
        "R": 1,
        "N": 0,
        "D": 0,
        "C": -3,
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
        "A": -1,
        "R": 0,
        "N": 0,
        "D": 2,
        "C": -4,
        "Q": 2,
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
        "A": 0,
        "R": -2,
        "N": 0,
        "D": -1,
        "C": -3,
        "Q": -2,
        "E": -2,
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
    "H": {
        "A": -2,
        "R": 0,
        "N": 1,
        "D": -1,
        "C": -3,
        "Q": 0,
        "E": 0,
        "G": -2,
        "H": 8,
        "I": -3,
        "L": -3,
        "K": -1,
        "M": -2,
        "F": -1,
        "P": -2,
        "S": -1,
        "T": -2,
        "W": -2,
        "Y": 2,
        "V": -3,
    },
    "I": {
        "A": -1,
        "R": -3,
        "N": -3,
        "D": -3,
        "C": -1,
        "Q": -3,
        "E": -3,
        "G": -4,
        "H": -3,
        "I": 4,
        "L": 2,
        "K": -3,
        "M": 1,
        "F": 0,
        "P": -3,
        "S": -2,
        "T": -1,
        "W": -3,
        "Y": -1,
        "V": 3,
    },
    "L": {
        "A": -1,
        "R": -2,
        "N": -3,
        "D": -4,
        "C": -1,
        "Q": -2,
        "E": -3,
        "G": -4,
        "H": -3,
        "I": 2,
        "L": 4,
        "K": -2,
        "M": 2,
        "F": 0,
        "P": -3,
        "S": -2,
        "T": -1,
        "W": -2,
        "Y": -1,
        "V": 1,
    },
    "K": {
        "A": -1,
        "R": 2,
        "N": 0,
        "D": -1,
        "C": -3,
        "Q": 1,
        "E": 1,
        "G": -2,
        "H": -1,
        "I": -3,
        "L": -2,
        "K": 5,
        "M": -1,
        "F": -3,
        "P": -1,
        "S": 0,
        "T": -1,
        "W": -3,
        "Y": -2,
        "V": -2,
    },
    "M": {
        "A": -1,
        "R": -1,
        "N": -2,
        "D": -3,
        "C": -1,
        "Q": 0,
        "E": -2,
        "G": -3,
        "H": -2,
        "I": 1,
        "L": 2,
        "K": -1,
        "M": 5,
        "F": 0,
        "P": -2,
        "S": -1,
        "T": -1,
        "W": -1,
        "Y": -1,
        "V": 1,
    },
    "F": {
        "A": -2,
        "R": -3,
        "N": -3,
        "D": -3,
        "C": -2,
        "Q": -3,
        "E": -3,
        "G": -3,
        "H": -1,
        "I": 0,
        "L": 0,
        "K": -3,
        "M": 0,
        "F": 6,
        "P": -4,
        "S": -2,
        "T": -2,
        "W": 1,
        "Y": 3,
        "V": -1,
    },
    "P": {
        "A": -1,
        "R": -2,
        "N": -2,
        "D": -1,
        "C": -3,
        "Q": -1,
        "E": -1,
        "G": -2,
        "H": -2,
        "I": -3,
        "L": -3,
        "K": -1,
        "M": -2,
        "F": -4,
        "P": 7,
        "S": -1,
        "T": -1,
        "W": -4,
        "Y": -3,
        "V": -2,
    },
    "S": {
        "A": 1,
        "R": -1,
        "N": 1,
        "D": 0,
        "C": -1,
        "Q": 0,
        "E": 0,
        "G": 0,
        "H": -1,
        "I": -2,
        "L": -2,
        "K": 0,
        "M": -1,
        "F": -2,
        "P": -1,
        "S": 4,
        "T": 1,
        "W": -3,
        "Y": -2,
        "V": -2,
    },
    "T": {
        "A": 0,
        "R": -1,
        "N": 0,
        "D": -1,
        "C": -1,
        "Q": -1,
        "E": -1,
        "G": -2,
        "H": -2,
        "I": -1,
        "L": -1,
        "K": -1,
        "M": -1,
        "F": -2,
        "P": -1,
        "S": 1,
        "T": 5,
        "W": -2,
        "Y": -2,
        "V": 0,
    },
    "W": {
        "A": -3,
        "R": -3,
        "N": -4,
        "D": -4,
        "C": -2,
        "Q": -2,
        "E": -3,
        "G": -2,
        "H": -2,
        "I": -3,
        "L": -2,
        "K": -3,
        "M": -1,
        "F": 1,
        "P": -4,
        "S": -3,
        "T": -2,
        "W": 11,
        "Y": 2,
        "V": -3,
    },
    "Y": {
        "A": -2,
        "R": -2,
        "N": -2,
        "D": -3,
        "C": -2,
        "Q": -1,
        "E": -2,
        "G": -3,
        "H": 2,
        "I": -1,
        "L": -1,
        "K": -2,
        "M": -1,
        "F": 3,
        "P": -3,
        "S": -2,
        "T": -2,
        "W": 2,
        "Y": 7,
        "V": -1,
    },
    "V": {
        "A": 0,
        "R": -3,
        "N": -3,
        "D": -3,
        "C": -1,
        "Q": -2,
        "E": -2,
        "G": -3,
        "H": -3,
        "I": 3,
        "L": 1,
        "K": -2,
        "M": 1,
        "F": -1,
        "P": -2,
        "S": -2,
        "T": 0,
        "W": -3,
        "Y": -1,
        "V": 4,
    },
}
BLOSUM62 = _blosum62_data


def matrix_align(
    seq1: str,
    seq2: str,
    matrix: Optional[Dict[str, Dict[str, int]]] = None,
    gap_open: int = -11,
    gap_extend: int = -1,
    mode: str = "global",
) -> Dict[str, Any]:
    """Align two sequences using a substitution matrix with affine gap penalties.

    Args:
        seq1: First protein sequence
        seq2: Second protein sequence
        matrix: Substitution matrix (default: BLOSUM62)
        gap_open: Gap opening penalty
        gap_extend: Gap extension penalty
        mode: "global" or "local"

    Returns:
        Dictionary with alignment results
    """
    if matrix is None:
        matrix = BLOSUM62

    seq1 = seq1.upper()
    seq2 = seq2.upper()
    m, n = len(seq1), len(seq2)

    if m == 0 or n == 0:
        return {
            "aligned_seq1": seq1 or "-" * n,
            "aligned_seq2": seq2 or "-" * m,
            "score": 0,
            "identity": 0.0,
            "method": mode,
            "matrix": "BLOSUM62",
        }

    NEG_INF = float("-inf")

    # Three matrices for affine gap penalties: M (match), X (gap in seq2), Y (gap in seq1)
    M_mat = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    X_mat = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    Y_mat = [[NEG_INF] * (n + 1) for _ in range(m + 1)]

    # Initialize
    if mode == "global":
        M_mat[0][0] = 0
        for i in range(1, m + 1):
            X_mat[i][0] = gap_open + (i - 1) * gap_extend
        for j in range(1, n + 1):
            Y_mat[0][j] = gap_open + (j - 1) * gap_extend
    else:
        for i in range(m + 1):
            M_mat[i][0] = 0
        for j in range(n + 1):
            M_mat[0][j] = 0

    max_score = 0
    max_i, max_j = 0, 0

    # Fill
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            aa1, aa2 = seq1[i - 1], seq2[j - 1]
            sub_score = matrix.get(aa1, {}).get(aa2, -1)

            best_prev = max(M_mat[i - 1][j - 1], X_mat[i - 1][j - 1], Y_mat[i - 1][j - 1])
            m_val = best_prev + sub_score

            x_val = max(M_mat[i - 1][j] + gap_open, X_mat[i - 1][j] + gap_extend)
            y_val = max(M_mat[i][j - 1] + gap_open, Y_mat[i][j - 1] + gap_extend)

            if mode == "local":
                M_mat[i][j] = max(0, m_val)
                X_mat[i][j] = max(0, x_val)
                Y_mat[i][j] = max(0, y_val)
                cell_max = max(M_mat[i][j], X_mat[i][j], Y_mat[i][j])
                if cell_max > max_score:
                    max_score = cell_max
                    max_i, max_j = i, j
            else:
                M_mat[i][j] = m_val
                X_mat[i][j] = x_val
                Y_mat[i][j] = y_val

    # Traceback
    aligned1, aligned2 = [], []

    if mode == "global":
        i, j = m, n
        score = max(M_mat[m][n], X_mat[m][n], Y_mat[m][n])
    else:
        i, j = max_i, max_j
        score = max_score

    while i > 0 or j > 0:
        if mode == "local" and max(M_mat[i][j], X_mat[i][j], Y_mat[i][j]) <= 0:
            break

        best = max(M_mat[i][j], X_mat[i][j], Y_mat[i][j])

        if i > 0 and j > 0 and best == M_mat[i][j]:
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and best == X_mat[i][j]:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        elif j > 0:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1
        else:
            break

    aligned1.reverse()
    aligned2.reverse()

    a1 = "".join(aligned1)
    a2 = "".join(aligned2)
    identity = calculate_alignment_identity({"aligned_seq1": a1, "aligned_seq2": a2})

    return {
        "aligned_seq1": a1,
        "aligned_seq2": a2,
        "score": score if score != NEG_INF else 0,
        "identity": identity,
        "method": mode,
        "matrix": "BLOSUM62",
        "gap_open": gap_open,
        "gap_extend": gap_extend,
    }


def multi_sequence_alignment(sequences: List[str], method: str = "progressive") -> Dict[str, Any]:
    """Perform multiple sequence alignment using progressive alignment.

    Uses a guide tree based on pairwise distances to align sequences
    progressively.

    Args:
        sequences: List of protein sequences
        method: Alignment method ("progressive")

    Returns:
        Dictionary with MSA results
    """
    if not sequences:
        return {"aligned_sequences": [], "score": 0, "n_sequences": 0}

    if len(sequences) == 1:
        return {"aligned_sequences": sequences, "score": 0, "n_sequences": 1}

    n = len(sequences)

    # Step 1: Compute pairwise alignment scores
    distance_matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            result = global_align(sequences[i], sequences[j])
            identity = result.get("identity", 0.0)
            dist = 1.0 - identity
            distance_matrix[i][j] = dist
            distance_matrix[j][i] = dist

    # Step 2: Build guide tree using neighbor-joining (simplified UPGMA)
    # Start by aligning the two closest sequences
    min_dist = float("inf")
    best_i, best_j = 0, 1
    for i in range(n):
        for j in range(i + 1, n):
            if distance_matrix[i][j] < min_dist:
                min_dist = distance_matrix[i][j]
                best_i, best_j = i, j

    # Step 3: Progressive alignment
    # Align closest pair first, then add remaining sequences
    result = global_align(sequences[best_i], sequences[best_j])
    aligned = [result["aligned_seq1"], result["aligned_seq2"]]
    used = {best_i, best_j}
    order = [best_i, best_j]

    # Add remaining sequences one by one (closest to current profile)
    while len(used) < n:
        best_next = -1
        best_score = float("-inf")

        for k in range(n):
            if k in used:
                continue
            # Score against first aligned sequence as profile representative
            r = global_align(aligned[0].replace("-", ""), sequences[k])
            if r["score"] > best_score:
                best_score = r["score"]
                best_next = k

        if best_next == -1:
            break

        # Align new sequence against profile representative
        r = global_align(aligned[0].replace("-", ""), sequences[best_next])

        # Adjust existing alignments to accommodate new gaps
        new_aligned = []
        profile_aligned = r["aligned_seq1"]
        new_seq_aligned = r["aligned_seq2"]

        # Map gaps from new alignment back to existing MSA
        gap_positions = [i for i, c in enumerate(profile_aligned) if c == "-"]
        for existing in aligned:
            adjusted = list(existing)
            for pos in sorted(gap_positions):
                if pos <= len(adjusted):
                    adjusted.insert(pos, "-")
            new_aligned.append("".join(adjusted))

        new_aligned.append(new_seq_aligned)
        aligned = new_aligned
        used.add(best_next)
        order.append(best_next)

    # Pad to same length
    max_len = max(len(s) for s in aligned)
    aligned = [s.ljust(max_len, "-") for s in aligned]

    return {
        "aligned_sequences": aligned,
        "n_sequences": len(aligned),
        "alignment_length": max_len,
        "order": order,
    }


def needleman_wunsch(seq1: str, seq2: str, match: int = 1, mismatch: int = -1, gap: int = -2) -> Tuple[int, str, str]:
    """Needleman-Wunsch global alignment returning (score, aligned1, aligned2).

    Convenience wrapper around global_align for tuple-style return.

    Args:
        seq1: First sequence
        seq2: Second sequence
        match: Match score
        mismatch: Mismatch penalty
        gap: Gap penalty

    Returns:
        Tuple of (score, aligned_seq1, aligned_seq2)
    """
    result = global_align(seq1, seq2, match=match, mismatch=mismatch, gap=gap)
    return result["score"], result["aligned_seq1"], result["aligned_seq2"]


def pairwise_identity(seq1: str, seq2: str) -> float:
    """Calculate pairwise sequence identity between two sequences.

    Args:
        seq1: First sequence
        seq2: Second sequence

    Returns:
        Identity fraction (0.0 to 1.0)
    """
    result = global_align(seq1, seq2)
    return result["identity"]
