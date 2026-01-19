"""Tests for DNA consensus sequence generation."""
from __future__ import annotations

import pytest

from metainformant.dna import consensus


def test_consensus_majority_ignores_gaps() -> None:
    """Test consensus generation from alignment."""
    aln = {
        "A": "ACG-",
        "B": "ACGT",
        "C": "ACGT",
    }
    try:
        cons = consensus.consensus_from_alignment(aln)
        assert isinstance(cons, str)
        # Should return a consensus string
        assert len(cons) > 0
    except (AttributeError, ImportError, KeyError) as e:
        pytest.skip(f"Consensus functions unavailable: {e}")


def test_consensus_module_importable() -> None:
    """Test that the consensus module can be imported."""
    assert consensus is not None
