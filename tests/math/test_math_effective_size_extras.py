from __future__ import annotations

from metainformant.math.population_genetics.statistics import effective_size_from_family_size_variance


def test_effective_size_from_family_size_variance():
    N = 1000.0
    Vk = 2.0
    Ne = effective_size_from_family_size_variance(N, Vk)
    assert abs(Ne - ((4.0 * N - 2.0) / (Vk + 2.0))) < 1e-12
