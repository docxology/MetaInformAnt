from __future__ import annotations

from metainformant.dna import distances


def test_kmer_distance_identical_zero_and_diff_positive() -> None:
    a = "ACGTACGT"
    b = "ACGTACGT"
    c = "TTTTTTTT"
    d_ab = distances.kmer_distance(a, b, k=2, metric="cosine")
    d_ac = distances.kmer_distance(a, c, k=2, metric="cosine")
    assert abs(d_ab) < 1e-12
    assert d_ac > 0
