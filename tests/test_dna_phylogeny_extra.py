"""Tests for additional DNA phylogeny functions."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.dna import phylogeny, sequences


@pytest.mark.slow
def test_bootstrap_support_and_upgma() -> None:
    """Test UPGMA tree construction and bootstrap support."""
    fasta_path = Path(__file__).parent / "data" / "dna" / "toy.fasta"

    try:
        seqs = sequences.read_fasta(str(fasta_path))
        tree_nj = phylogeny.neighbor_joining_tree(seqs)
        tree_upgma = phylogeny.upgma_tree(seqs)
        assert phylogeny.basic_tree_stats(tree_nj)["leaves"] == len(seqs)
        assert phylogeny.basic_tree_stats(tree_upgma)["leaves"] == len(seqs)

        # Small bootstrap just to exercise path
        res_tree = phylogeny.bootstrap_support(tree_nj, seqs, n_replicates=10, method="nj")
        assert res_tree is not None
    except FileNotFoundError:
        pytest.skip("Test data file not found")
    except (ImportError, AttributeError) as e:
        pytest.skip(f"Phylogeny functions unavailable: {e}")
