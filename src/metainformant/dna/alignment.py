"""Sequence alignment utilities for DNA sequences.

This module provides implementations of global (Needleman-Wunsch) and local
(Smith-Waterman) sequence alignment algorithms, along with alignment analysis
functions for identity calculation, conserved regions detection, and alignment
statistics.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)


@dataclass
class AlignmentResult:
    """Result of a sequence alignment operation.

    Attributes:
        seq1_aligned: First sequence with gaps inserted
        seq2_aligned: Second sequence with gaps inserted
        score: Alignment score
        score_matrix: Scoring matrix used for alignment
        traceback_matrix: Traceback matrix for path reconstruction
        start_positions: Tuple of (seq1_start, seq2_start) for local alignments
    """
    seq1_aligned: str
    seq2_aligned: str
    score: float
    score_matrix: np.ndarray
    traceback_matrix: np.ndarray
    start_positions: Optional[Tuple[int, int]] = None


def global_align(seq1: str, seq2: str, match: int = 1, mismatch: int = -1, gap: int = -2) -> AlignmentResult:
    """Perform global sequence alignment using Needleman-Wunsch algorithm.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence
        match: Score for matching nucleotides (default: 1)
        mismatch: Score for mismatching nucleotides (default: -1)
        gap: Score for gap insertion (default: -2)

    Returns:
        AlignmentResult containing aligned sequences and metadata

    Example:
        >>> result = global_align("ATCG", "ATCG")
        >>> result.score
        4
        >>> result.seq1_aligned
        'ATCG'
        >>> result.seq2_aligned
        'ATCG'
    """
    if not seq1 or not seq2:
        raise ValueError("Sequences cannot be empty")

    m, n = len(seq1), len(seq2)

    # Initialize scoring matrix
    score_matrix = np.zeros((m + 1, n + 1), dtype=float)
    traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)

    # Directions for traceback: 0=diagonal, 1=up, 2=left
    DIAGONAL, UP, LEFT = 0, 1, 2

    # Initialize first row and column with gap penalties
    for i in range(1, m + 1):
        score_matrix[i, 0] = score_matrix[i-1, 0] + gap
        traceback_matrix[i, 0] = UP

    for j in range(1, n + 1):
        score_matrix[0, j] = score_matrix[0, j-1] + gap
        traceback_matrix[0, j] = LEFT

    # Fill the matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculate scores for each direction
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            diagonal_score = score_matrix[i-1, j-1] + match_score
            up_score = score_matrix[i-1, j] + gap
            left_score = score_matrix[i, j-1] + gap

            # Choose the maximum score
            max_score = max(diagonal_score, up_score, left_score)
            score_matrix[i, j] = max_score

            # Set traceback direction
            if max_score == diagonal_score:
                traceback_matrix[i, j] = DIAGONAL
            elif max_score == up_score:
                traceback_matrix[i, j] = UP
            else:
                traceback_matrix[i, j] = LEFT

    # Traceback to get aligned sequences
    aligned_seq1, aligned_seq2 = _traceback_global(seq1, seq2, traceback_matrix)

    return AlignmentResult(
        seq1_aligned=aligned_seq1,
        seq2_aligned=aligned_seq2,
        score=score_matrix[m, n],
        score_matrix=score_matrix,
        traceback_matrix=traceback_matrix
    )


def local_align(seq1: str, seq2: str, match: int = 1, mismatch: int = -1, gap: int = -2) -> AlignmentResult:
    """Perform local sequence alignment using Smith-Waterman algorithm.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence
        match: Score for matching nucleotides (default: 1)
        mismatch: Score for mismatching nucleotides (default: -1)
        gap: Score for gap insertion (default: -2)

    Returns:
        AlignmentResult containing aligned sequences and metadata

    Example:
        >>> result = local_align("ATCGATCG", "GCTAGCTA")
        >>> result.score >= 0  # Local alignment always has non-negative score
        True
    """
    if not seq1 or not seq2:
        raise ValueError("Sequences cannot be empty")

    m, n = len(seq1), len(seq2)

    # Initialize scoring matrix
    score_matrix = np.zeros((m + 1, n + 1), dtype=float)
    traceback_matrix = np.zeros((m + 1, n + 1), dtype=int)

    # Directions for traceback: 0=diagonal, 1=up, 2=left, 3=stop
    DIAGONAL, UP, LEFT, STOP = 0, 1, 2, 3

    max_score = 0.0
    max_i, max_j = 0, 0

    # Fill the matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculate scores for each direction
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            diagonal_score = score_matrix[i-1, j-1] + match_score
            up_score = score_matrix[i-1, j] + gap
            left_score = score_matrix[i, j-1] + gap

            # Choose the maximum score (including 0 for local alignment)
            max_score_cell = max(diagonal_score, up_score, left_score, 0.0)
            score_matrix[i, j] = max_score_cell

            # Set traceback direction
            if max_score_cell == 0.0:
                traceback_matrix[i, j] = STOP
            elif max_score_cell == diagonal_score:
                traceback_matrix[i, j] = DIAGONAL
            elif max_score_cell == up_score:
                traceback_matrix[i, j] = UP
            else:
                traceback_matrix[i, j] = LEFT

            # Track maximum score position
            if max_score_cell > max_score:
                max_score = max_score_cell
                max_i, max_j = i, j

    # If no positive alignment found, return empty alignment
    if max_score == 0.0:
        return AlignmentResult(
            seq1_aligned="",
            seq2_aligned="",
            score=0.0,
            score_matrix=score_matrix,
            traceback_matrix=traceback_matrix,
            start_positions=(0, 0)
        )

    # Traceback from maximum score position
    aligned_seq1, aligned_seq2, start_pos = _traceback_local(seq1, seq2, traceback_matrix, max_i, max_j)

    return AlignmentResult(
        seq1_aligned=aligned_seq1,
        seq2_aligned=aligned_seq2,
        score=max_score,
        score_matrix=score_matrix,
        traceback_matrix=traceback_matrix,
        start_positions=start_pos
    )


def calculate_alignment_identity(alignment: AlignmentResult) -> float:
    """Calculate the identity percentage between aligned sequences.

    Args:
        alignment: AlignmentResult from global_align or local_align

    Returns:
        Identity as a percentage (0.0 to 100.0)

    Example:
        >>> result = global_align("ATCG", "ATCG")
        >>> calculate_alignment_identity(result)
        100.0
    """
    if not alignment.seq1_aligned or not alignment.seq2_aligned:
        return 0.0

    if len(alignment.seq1_aligned) != len(alignment.seq2_aligned):
        raise ValueError("Aligned sequences must have equal length")

    matches = 0
    total_positions = 0

    for a, b in zip(alignment.seq1_aligned, alignment.seq2_aligned):
        if a != '-' and b != '-':  # Only count non-gap positions
            total_positions += 1
            if a == b:
                matches += 1

    return (matches / total_positions * 100) if total_positions > 0 else 0.0


def find_conserved_regions(alignment: AlignmentResult, min_length: int = 5) -> List[Tuple[str, int, int]]:
    """Find conserved regions in an alignment.

    Args:
        alignment: AlignmentResult from global_align or local_align
        min_length: Minimum length of conserved region (default: 5)

    Returns:
        List of tuples (sequence, start_pos, end_pos) for conserved regions

    Example:
        >>> result = global_align("ATCGATCG", "ATCGATCG")
        >>> regions = find_conserved_regions(result, min_length=4)
        >>> len(regions) > 0
        True
    """
    if not alignment.seq1_aligned or not alignment.seq2_aligned:
        return []

    if len(alignment.seq1_aligned) != len(alignment.seq2_aligned):
        raise ValueError("Aligned sequences must have equal length")

    conserved_regions = []
    current_region_start = None
    current_region_seq = []

    for i, (a, b) in enumerate(zip(alignment.seq1_aligned, alignment.seq2_aligned)):
        if a == b and a != '-' and b != '-':
            # Matching position
            if current_region_start is None:
                current_region_start = i
            current_region_seq.append(a)
        else:
            # End of conserved region
            if current_region_start is not None and len(current_region_seq) >= min_length:
                conserved_regions.append((
                    ''.join(current_region_seq),
                    current_region_start,
                    i - 1
                ))
            current_region_start = None
            current_region_seq = []

    # Handle region that extends to end
    if current_region_start is not None and len(current_region_seq) >= min_length:
        conserved_regions.append((
            ''.join(current_region_seq),
            current_region_start,
            len(alignment.seq1_aligned) - 1
        ))

    return conserved_regions


def alignment_statistics(alignment: AlignmentResult) -> Dict[str, float]:
    """Calculate comprehensive statistics for an alignment.

    Args:
        alignment: AlignmentResult from global_align or local_align

    Returns:
        Dictionary containing alignment statistics

    Example:
        >>> result = global_align("ATCG", "ATCG")
        >>> stats = alignment_statistics(result)
        >>> stats['identity'] == 100.0
        True
    """
    if not alignment.seq1_aligned or not alignment.seq2_aligned:
        return {
            'identity': 0.0,
            'similarity': 0.0,
            'gaps': 0,
            'gap_percentage': 0.0,
            'alignment_length': 0,
            'matches': 0,
            'mismatches': 0
        }

    if len(alignment.seq1_aligned) != len(alignment.seq2_aligned):
        raise ValueError("Aligned sequences must have equal length")

    alignment_length = len(alignment.seq1_aligned)
    matches = 0
    mismatches = 0
    gaps = 0

    # Count matches, mismatches, and gaps
    for a, b in zip(alignment.seq1_aligned, alignment.seq2_aligned):
        if a == '-' or b == '-':
            gaps += 1
        elif a == b:
            matches += 1
        else:
            mismatches += 1

    # Calculate percentages
    identity = (matches / (alignment_length - gaps) * 100) if (alignment_length - gaps) > 0 else 0.0
    gap_percentage = (gaps / alignment_length * 100) if alignment_length > 0 else 0.0

    return {
        'identity': identity,
        'similarity': identity,  # For DNA, identity equals similarity
        'gaps': gaps,
        'gap_percentage': gap_percentage,
        'alignment_length': alignment_length,
        'matches': matches,
        'mismatches': mismatches
    }


def _traceback_global(seq1: str, seq2: str, traceback_matrix: np.ndarray) -> Tuple[str, str]:
    """Perform traceback for global alignment."""
    m, n = len(seq1), len(seq2)
    aligned_seq1 = []
    aligned_seq2 = []

    i, j = m, n
    DIAGONAL, UP, LEFT = 0, 1, 2

    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback_matrix[i, j] == DIAGONAL:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and traceback_matrix[i, j] == UP:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        elif j > 0 and traceback_matrix[i, j] == LEFT:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1
        else:
            # Should not reach here in global alignment
            break

    # Reverse the sequences since we built them backwards
    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


def _traceback_local(seq1: str, seq2: str, traceback_matrix: np.ndarray, start_i: int, start_j: int) -> Tuple[str, str, Tuple[int, int]]:
    """Perform traceback for local alignment."""
    aligned_seq1 = []
    aligned_seq2 = []

    i, j = start_i, start_j
    DIAGONAL, UP, LEFT, STOP = 0, 1, 2, 3

    while traceback_matrix[i, j] != STOP:
        if traceback_matrix[i, j] == DIAGONAL:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif traceback_matrix[i, j] == UP:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        elif traceback_matrix[i, j] == LEFT:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1
        else:
            break

        if i < 0 or j < 0:
            break

    # Reverse the sequences since we built them backwards
    aligned_seq1.reverse()
    aligned_seq2.reverse()

    # Calculate start positions in original sequences
    start_pos_seq1 = i if i >= 0 else 0
    start_pos_seq2 = j if j >= 0 else 0

    return ''.join(aligned_seq1), ''.join(aligned_seq2), (start_pos_seq1, start_pos_seq2)








