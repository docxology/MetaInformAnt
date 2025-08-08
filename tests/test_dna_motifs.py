from __future__ import annotations

from metainformant.dna import motifs


def test_find_motif_positions_iupac() -> None:
    seq = "GGCATGATGCAATG"
    # ATN should match ATG and ATA etc
    pos = motifs.find_motif_positions(seq, "ATN")
    # positions of 'ATG' occur at index 3 and 6 and 10; ATN also includes index 11? (ATG at 3, 6, 10)
    assert 3 in pos and 6 in pos and 10 in pos


