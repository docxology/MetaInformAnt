from __future__ import annotations

from dataclasses import dataclass

from Bio.Align import PairwiseAligner


@dataclass
class AlignmentResult:
    aligned_seq1: str
    aligned_seq2: str
    score: float


def global_align(
    seq1: str,
    seq2: str,
    *,
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    gap_score: float = -2.0,
    max_alignments: int = 1
) -> AlignmentResult:
    """Global alignment using Needleman-Wunsch algorithm with configurable scoring.

    Args:
        seq1: First sequence to align
        seq2: Second sequence to align
        match_score: Score for matching bases
        mismatch_score: Penalty for mismatching bases
        gap_score: Penalty for gaps
        max_alignments: Maximum number of alignments to return (best first)

    Returns:
        AlignmentResult with aligned sequences and score
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.gap_score = gap_score

    alignments = aligner.align(seq1, seq2)
    if len(alignments) == 0:
        # Return identity alignment if no alignment found
        return AlignmentResult(
            aligned_seq1=seq1,
            aligned_seq2=seq2,
            score=0.0
        )

    best_alignment = alignments[0]
    lines = best_alignment.format().splitlines()

    # Extract aligned sequences (handle case where format might differ)
    if len(lines) >= 3:
        gapped1 = lines[0].strip()
        gapped2 = lines[2].strip()
    else:
        # Fallback to original sequences
        gapped1 = seq1
        gapped2 = seq2

    return AlignmentResult(
        aligned_seq1=gapped1,
        aligned_seq2=gapped2,
        score=float(best_alignment.score)
    )


def local_align(seq1: str, seq2: str) -> AlignmentResult:
    aligner = PairwiseAligner()
    aligner.mode = "local"
    a = aligner.align(seq1, seq2)[0]
    lines = a.format().splitlines()
    gapped1 = lines[0].strip()
    gapped2 = lines[2].strip()
    return AlignmentResult(aligned_seq1=gapped1, aligned_seq2=gapped2, score=a.score)
