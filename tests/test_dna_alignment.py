from __future__ import annotations

from pathlib import Path

from metainformant.dna import alignment, sequences


def test_global_alignment_score_and_lengths(tmp_path: Path) -> None:
    fasta_path = Path(__file__).parent / "data" / "dna" / "toy.fasta"
    seqs = sequences.read_fasta(str(fasta_path))
    a = seqs["A"]
    b = seqs["B"]

    result = alignment.global_align(a, b)

    assert result.score > 0
    assert len(result.aligned_seq1) == len(result.aligned_seq2)


def test_local_alignment_non_trivial() -> None:
    s1 = "ACGTACGT"
    s2 = "TTACGTAAT"
    result = alignment.local_align(s1, s2)
    assert result.score > 0
    assert len(result.aligned_seq1) == len(result.aligned_seq2)


