from __future__ import annotations

from pathlib import Path

from metainformant.dna import sequences, phylogeny


def test_neighbor_joining_tree_from_kmer_distance() -> None:
    fasta_path = Path(__file__).parent / "data" / "dna" / "toy.fasta"
    seqs = sequences.read_fasta(str(fasta_path))
    tree = phylogeny.nj_tree_from_kmer(seqs, k=2, metric="cosine")
    assert tree is not None
    assert len(tree.get_terminals()) == len(seqs)


