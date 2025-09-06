from __future__ import annotations

from pathlib import Path

from metainformant.dna import msa, sequences


def test_align_msa_returns_equal_length_alignment() -> None:
    fasta_path = Path(__file__).parent / "data" / "dna" / "toy.fasta"
    seqs = sequences.read_fasta(str(fasta_path))
    aln = msa.align_msa(seqs, method="auto")
    assert set(aln.keys()) == set(seqs.keys())
    lengths = {len(s) for s in aln.values()}
    assert len(lengths) == 1  # all aligned sequences have same length
