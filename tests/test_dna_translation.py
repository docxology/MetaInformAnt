from __future__ import annotations

from metainformant.dna import translation


def test_translate_basic() -> None:
    assert translation.translate_dna("ATGGCC") == "MA"
    assert translation.translate_dna("ATGTAA", to_stop=True) == "M"


def test_find_orfs_simple() -> None:
    seq = "AAAATGAAATTTTGAATGCCCCCTAG"
    orfs = translation.find_orfs(seq, min_aa=2, include_reverse=False)
    # expect at least two ORFs (ATG...TGA and ATG...TAG)
    assert any(orf.protein.startswith("MK") for orf in orfs)
    assert any(orf.stop_codon in {"TGA", "TAG"} for orf in orfs)
