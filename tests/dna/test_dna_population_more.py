from __future__ import annotations

from metainformant.dna import population


def test_segregating_sites_and_watterson_theta() -> None:
    seqs = [
        "ACGT",
        "ACGA",
        "ACGA",
    ]
    S = population.segregating_sites(seqs)
    assert S >= 1
    theta_w = population.wattersons_theta(seqs)
    assert theta_w > 0
