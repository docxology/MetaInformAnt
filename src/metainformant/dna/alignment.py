from __future__ import annotations

from dataclasses import dataclass

from Bio import pairwise2


@dataclass
class AlignmentResult:
    aligned_seq1: str
    aligned_seq2: str
    score: float


def global_align(seq1: str, seq2: str) -> AlignmentResult:
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    a = alignments[0]
    return AlignmentResult(aligned_seq1=a.seqA, aligned_seq2=a.seqB, score=a.score)


def local_align(seq1: str, seq2: str) -> AlignmentResult:
    alignments = pairwise2.align.localxx(seq1, seq2, one_alignment_only=True)
    a = alignments[0]
    return AlignmentResult(aligned_seq1=a.seqA, aligned_seq2=a.seqB, score=a.score)


