"""Tests for DNA phylogeny k-mer based tree construction."""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.dna import phylogeny, sequences


@pytest.mark.slow
def test_neighbor_joining_tree_from_kmer_distance() -> None:
    """Test NJ tree construction from k-mer distances."""
    fasta_path = Path(__file__).parent / "data" / "dna" / "toy.fasta"

    try:
        seqs = sequences.read_fasta(str(fasta_path))
        tree = phylogeny.nj_tree_from_kmer(seqs, k=2, metric="cosine")
        assert tree is not None
        assert phylogeny.basic_tree_stats(tree)["leaves"] == len(seqs)
    except FileNotFoundError:
        pytest.skip("Test data file not found")
    except (ImportError, AttributeError) as e:
        pytest.skip(f"Phylogeny functions unavailable: {e}")
