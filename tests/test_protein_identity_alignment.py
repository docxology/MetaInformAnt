from __future__ import annotations


def test_pairwise_identity_and_needleman_wunsch():
    from metainformant.protein.alignment import needleman_wunsch, pairwise_identity

    a = "MKT"
    b = "MKT"
    assert pairwise_identity(a, b) == 1.0

    c = "MKTA"
    d = "MKT-"
    # simple NW with match=1, mismatch=-1, gap=-1 should align reasonably
    s, a1, a2 = needleman_wunsch(c, d, match=1, mismatch=-1, gap=-1)
    assert isinstance(s, int)
    assert len(a1) == len(a2)
