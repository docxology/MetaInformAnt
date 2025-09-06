from __future__ import annotations

from metainformant.dna import consensus


def test_consensus_majority_ignores_gaps() -> None:
    aln = {
        "A": "ACG-",
        "B": "ACGT",
        "C": "ACGT",
    }
    cons = consensus.consensus_from_alignment(aln)
    assert cons == "ACGT"
