from __future__ import annotations

from pathlib import Path

from metainformant.dna import sequences


def test_read_fasta_parses_ids() -> None:
    fasta_path = Path(__file__).parent / "data" / "dna" / "toy.fasta"
    seqs = sequences.read_fasta(str(fasta_path))
    assert set(seqs.keys()) == {"A", "B", "C"}
    assert all(isinstance(v, str) and len(v) > 0 for v in seqs.values())
