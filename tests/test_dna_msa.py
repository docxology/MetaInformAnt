"""Tests for DNA multiple sequence alignment.

These tests use external alignment tools (muscle, mafft, clustalo) when available.
Tests will skip gracefully if tools are not installed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.dna import msa, sequences


@pytest.mark.external_tool
def test_align_msa_returns_equal_length_alignment() -> None:
    """Test MSA alignment returns equal-length sequences."""
    fasta_path = Path(__file__).parent / "data" / "dna" / "toy.fasta"

    try:
        seqs = sequences.read_fasta(str(fasta_path))
        aln = msa.align_msa(seqs, method="auto")
        assert set(aln.keys()) == set(seqs.keys())
        lengths = {len(s) for s in aln.values()}
        assert len(lengths) == 1  # all aligned sequences have same length
    except FileNotFoundError:
        pytest.skip("Test data file not found")
    except (ImportError, AttributeError) as e:
        pytest.skip(f"MSA functions unavailable: {e}")


@pytest.mark.external_tool
def test_align_msa_with_tmp_file(tmp_path: Path) -> None:
    """Test MSA alignment with temporary FASTA file."""
    fasta_content = """>seq1
ATCGATCGATCG
>seq2
ATCGATCGATCG
>seq3
TTTTTTTTTTTT
"""
    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(fasta_content)

    seqs = sequences.read_fasta(str(fasta_path))
    assert len(seqs) == 3

    try:
        aln = msa.align_msa(seqs)
        assert isinstance(aln, dict)
    except (ImportError, AttributeError) as e:
        pytest.skip(f"MSA functions unavailable: {e}")
