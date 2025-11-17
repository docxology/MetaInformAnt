from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.dna import phylogeny, sequences


@pytest.mark.slow
def test_neighbor_joining_tree_from_toy_fasta() -> None:
    fasta_path = Path(__file__).parent / "data" / "dna" / "toy.fasta"
    seqs = sequences.read_fasta(str(fasta_path))

    tree = phylogeny.neighbor_joining_tree(seqs)
    assert tree is not None
    # number of terminals should equal number of input sequences
    assert len(tree.get_terminals()) == len(seqs)

    newick = phylogeny.to_newick(tree)
    assert isinstance(newick, str) and ";" in newick
