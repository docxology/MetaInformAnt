from __future__ import annotations

from pathlib import Path

from metainformant.dna import sequences, phylogeny


def test_bootstrap_support_deterministic_with_seed() -> None:
    fasta_path = Path(__file__).parent / "data" / "dna" / "toy.fasta"
    seqs = sequences.read_fasta(str(fasta_path))
    s1 = phylogeny.bootstrap_support(seqs, n_replicates=25, method="nj", random_state=42)
    s2 = phylogeny.bootstrap_support(seqs, n_replicates=25, method="nj", random_state=42)
    assert s1 == s2
    assert all(0.0 <= v <= 1.0 for v in s1.values())


