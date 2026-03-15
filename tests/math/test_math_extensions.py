from __future__ import annotations

from metainformant.math.population_genetics.coalescent import watterson_theta
from metainformant.math.population_genetics.ld import ld_decay_r2
from metainformant.math.quantitative_genetics.core import realized_heritability


def test_ld_decay_and_watterson_theta_and_realized_h2():
    r2t = ld_decay_r2(0.8, recombination_rate=0.1, generations=5)
    # (1-0.1)^(10) = 0.9^10 ≈ 0.34867844; 0.8 * that ≈ 0.27894275
    assert abs(r2t - 0.8 * (0.9**10)) < 1e-12

    theta = watterson_theta(num_segregating_sites=10, sample_size=5)
    # a1 = 1 + 1/2 + 1/3 + 1/4 = 2.083333...
    assert abs(theta - (10.0 / (1 + 0.5 + 1 / 3 + 0.25))) < 1e-12

    h2r = realized_heritability(response=0.6, selection_differential=1.2)
    assert abs(h2r - 0.5) < 1e-12
