from __future__ import annotations

from dataclasses import dataclass

from Bio.Align import PairwiseAligner


@dataclass
class AlignmentResult:
    aligned_seq1: str
    aligned_seq2: str
    score: float


def global_align(seq1: str, seq2: str) -> AlignmentResult:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    a = aligner.align(seq1, seq2)[0]
    lines = a.format().splitlines()
    gapped1 = lines[0].strip()
    gapped2 = lines[2].strip()
    return AlignmentResult(aligned_seq1=gapped1, aligned_seq2=gapped2, score=a.score)


def local_align(seq1: str, seq2: str) -> AlignmentResult:
    aligner = PairwiseAligner()
    aligner.mode = "local"
    a = aligner.align(seq1, seq2)[0]
    lines = a.format().splitlines()
    gapped1 = lines[0].strip()
    gapped2 = lines[2].strip()
    return AlignmentResult(aligned_seq1=gapped1, aligned_seq2=gapped2, score=a.score)
