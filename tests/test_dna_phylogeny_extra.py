from __future__ import annotations

from pathlib import Path

from metainformant.dna import phylogeny, sequences


def test_bootstrap_support_and_upgma() -> None:
    fasta_path = Path(__file__).parent / "data" / "dna" / "toy.fasta"
    seqs = sequences.read_fasta(str(fasta_path))
    tree_nj = phylogeny.neighbor_joining_tree(seqs)
    tree_upgma = phylogeny.upgma_tree(seqs)
    assert len(tree_nj.get_terminals()) == len(seqs)
    assert len(tree_upgma.get_terminals()) == len(seqs)

    # small bootstrap just to exercise path
    supports = phylogeny.bootstrap_support(seqs, n_replicates=10, method="nj")
    assert isinstance(supports, dict)
