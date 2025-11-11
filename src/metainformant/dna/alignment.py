from __future__ import annotations

from dataclasses import dataclass

from Bio.Align import PairwiseAligner


@dataclass
class AlignmentResult:
    """Result of pairwise sequence alignment.
    
    Attributes:
        aligned_seq1: First sequence with gaps inserted
        aligned_seq2: Second sequence with gaps inserted
        score: Alignment score (higher is better)
    """
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

    # Extract aligned sequences from formatted output
    # Format: "target 0 ATCGATCG 8" -> extract "ATCGATCG"
    import re
    if len(lines) >= 3:
        # Extract sequence part (between numbers) using regex
        match1 = re.search(r'\d+\s+([A-Za-z-]+)\s+\d+', lines[0])
        match2 = re.search(r'\d+\s+([A-Za-z-]+)\s+\d+', lines[2])
        if match1 and match2:
            gapped1 = match1.group(1)
            gapped2 = match2.group(1)
        else:
            # Fallback: extract just sequence characters
            gapped1 = ''.join(c for c in lines[0] if c.isalpha() or c == '-')
            gapped2 = ''.join(c for c in lines[2] if c.isalpha() or c == '-')
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
    """Perform local (Smith-Waterman) sequence alignment.
    
    Args:
        seq1: First sequence to align
        seq2: Second sequence to align
        
    Returns:
        AlignmentResult with best local alignment
    """
    aligner = PairwiseAligner()
    aligner.mode = "local"
    a = aligner.align(seq1, seq2)[0]
    lines = a.format().splitlines()
    
    # Extract sequence part from formatted output
    import re
    if len(lines) >= 3:
        match1 = re.search(r'\d+\s+([A-Za-z-]+)\s+\d+', lines[0])
        match2 = re.search(r'\d+\s+([A-Za-z-]+)\s+\d+', lines[2])
        if match1 and match2:
            gapped1 = match1.group(1)
            gapped2 = match2.group(1)
        else:
            gapped1 = ''.join(c for c in lines[0] if c.isalpha() or c == '-')
            gapped2 = ''.join(c for c in lines[2] if c.isalpha() or c == '-')
    else:
        gapped1 = seq1
        gapped2 = seq2
    
    return AlignmentResult(aligned_seq1=gapped1, aligned_seq2=gapped2, score=a.score)


def calculate_alignment_identity(alignment: AlignmentResult) -> float:
    """Calculate percent identity of an alignment.

    Args:
        alignment: AlignmentResult from global_align or local_align

    Returns:
        Percent identity (0-100), calculated as matches / aligned positions (excluding gaps in both sequences)
    """
    seq1 = alignment.aligned_seq1
    seq2 = alignment.aligned_seq2

    if len(seq1) != len(seq2):
        return 0.0

    # Count matches and aligned positions (positions where neither sequence has a gap)
    matches = 0
    aligned_positions = 0

    for a, b in zip(seq1, seq2):
        if a != "-" and b != "-":
            aligned_positions += 1
            if a == b:
                matches += 1

    if aligned_positions == 0:
        return 0.0

    return (matches / aligned_positions) * 100


def find_conserved_regions(alignment: AlignmentResult, min_length: int = 5) -> list[tuple[str, int, int]]:
    """Find conserved regions in an alignment.

    Args:
        alignment: AlignmentResult from global_align or local_align
        min_length: Minimum length of conserved region

    Returns:
        List of (sequence, start, end) tuples for conserved regions
    """
    seq1 = alignment.aligned_seq1
    seq2 = alignment.aligned_seq2

    conserved_regions = []
    current_start = None

    for i in range(len(seq1)):
        if seq1[i] == seq2[i] and seq1[i] != "-":
            if current_start is None:
                current_start = i
        else:
            if current_start is not None:
                length = i - current_start
                if length >= min_length:
                    conserved_regions.append((seq1[current_start:i], current_start, i))
                current_start = None

    # Handle case where conserved region goes to end
    if current_start is not None:
        length = len(seq1) - current_start
        if length >= min_length:
            conserved_regions.append((seq1[current_start:], current_start, len(seq1)))

    return conserved_regions


def alignment_statistics(alignment: AlignmentResult) -> dict[str, float]:
    """Calculate comprehensive statistics for an alignment.

    Args:
        alignment: AlignmentResult from global_align or local_align

    Returns:
        Dictionary with alignment statistics
    """
    seq1 = alignment.aligned_seq1
    seq2 = alignment.aligned_seq2

    # Lengths
    len1 = len(seq1)
    len2 = len(seq2)

    # Gaps
    gaps1 = seq1.count("-")
    gaps2 = seq2.count("-")

    # Matches and mismatches (excluding gaps)
    matches = 0
    mismatches = 0

    for i in range(len1):
        char1 = seq1[i]
        char2 = seq2[i]
        if char1 != "-" and char2 != "-":
            if char1 == char2:
                matches += 1
            else:
                mismatches += 1

    # Identity (excluding gaps)
    total_positions = matches + mismatches
    identity = (matches / total_positions * 100) if total_positions > 0 else 0.0

    # Similarity (allowing for conservative substitutions)
    similarity = identity  # For DNA, identity and similarity are the same

    return {
        "length1": len1,
        "length2": len2,
        "gaps1": gaps1,
        "gaps2": gaps2,
        "matches": matches,
        "mismatches": mismatches,
        "identity": identity,
        "similarity": similarity,
        "score": alignment.score
    }
