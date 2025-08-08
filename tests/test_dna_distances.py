from __future__ import annotations

from metainformant.dna import distances


def test_p_distance_and_jc69() -> None:
    s1, s2 = "AAAA", "AAAT"
    p = distances.p_distance(s1, s2)
    assert abs(p - 0.25) < 1e-9
    d = distances.jc69_distance(s1, s2)
    assert abs(d - 0.304099) < 1e-3


