from __future__ import annotations

from metainformant.math import expected_pairwise_diversity, tajima_constants, tajimas_D


def test_expected_pairwise_diversity_and_tajima_constants():
    Ne = 1000.0
    mu = 1e-8
    pi = expected_pairwise_diversity(Ne, mu)
    assert abs(pi - (4.0 * Ne * mu)) < 1e-20

    const = tajima_constants(5)
    assert "a1" in const and abs(const["a1"] - (1 + 0.5 + 1 / 3 + 0.25)) < 1e-12


def test_tajimas_D_zero_when_pi_matches_S_over_a1():
    n = 10
    S = 12
    const = tajima_constants(n)
    pi = S / const["a1"] if const["a1"] > 0 else 0.0
    D = tajimas_D(S, pi, n)
    assert abs(D) < 1e-9
