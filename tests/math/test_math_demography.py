"""Tests for demographic models."""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.math.population_genetics.demography import (
    age_structure_model,
    bottleneck_effective_size,
    exponential_growth_effective_size,
    exponential_growth_model,
    island_model_update,
    logistic_growth_model,
    two_epoch_effective_size,
)


class TestExponentialGrowthEffectiveSize:
    """Wright's harmonic mean over the exponential_growth_model trajectory."""

    def test_no_growth(self):
        ne = exponential_growth_effective_size(1000, growth_rate=0.0, generations=10)
        assert ne == 1000.0

    def test_positive_growth_matches_geometric_series(self):
        # sum_{i<t} 1/(N0 (1+r)^i) = (1 - (1+r)^-t) / (N0 (1 - 1/(1+r)))
        n0, r, t = 10000, 0.23, 10
        expected = t * n0 * (1 - 1 / (1 + r)) / (1 - (1 + r) ** -t)
        ne = exponential_growth_effective_size(n0, growth_rate=r, generations=t)
        assert ne == pytest.approx(expected, rel=1e-12)
        # Above the founding size, but well below the arithmetic mean of the trajectory
        assert n0 < ne < sum(exponential_growth_model(n0, r, t - 1)) / t

    def test_negative_growth(self):
        ne = exponential_growth_effective_size(1000, growth_rate=-0.1, generations=10)
        assert ne == pytest.approx(10 / sum(1 / (1000 * 0.9**i) for i in range(10)), rel=1e-12)
        assert ne < 1000

    def test_zero_generations(self):
        ne = exponential_growth_effective_size(1000, growth_rate=0.1, generations=0)
        assert ne == 1000.0

    def test_rejects_rates_that_empty_the_population(self):
        with pytest.raises(ValueError, match="exceed -1"):
            exponential_growth_effective_size(1000, growth_rate=-1.0, generations=3)


class TestBottleneckEffectiveSize:
    """Harmonic mean over one generation at each size plus recovery generations."""

    def test_harmonic_mean_of_three_generations(self):
        ne = bottleneck_effective_size(10000, 100, 5)
        assert ne == pytest.approx(3 / (1 / 10000 + 1 / 100 + 1 / 5), rel=1e-12)
        # Dominated by the smallest generation
        assert 5 < ne < 15

    def test_mild_bottleneck(self):
        ne = bottleneck_effective_size(10000, 5000, 10000)
        assert ne == pytest.approx(3 / (1 / 10000 + 1 / 5000 + 1 / 10000), rel=1e-12)
        assert 5000 < ne < 10000

    def test_recovery_generations_raise_ne(self):
        ne_no_recovery = bottleneck_effective_size(10000, 100, 5000)
        ne_with_recovery = bottleneck_effective_size(10000, 100, 5000, recovery_generations=10)

        assert ne_with_recovery == pytest.approx(13 / (1 / 10000 + 1 / 100 + 11 / 5000), rel=1e-12)
        assert ne_with_recovery > ne_no_recovery

    def test_rejects_non_positive_sizes(self):
        with pytest.raises(ValueError, match="positive"):
            bottleneck_effective_size(10000, 100, 0)
        with pytest.raises(ValueError, match="negative"):
            bottleneck_effective_size(10000, 100, 50, recovery_generations=-1)


class TestTwoEpochEffectiveSize:
    """Coalescent effective size: the constant size with the same expected pairwise coalescence time."""

    def test_recent_expansion_still_reflects_the_ancient_size(self):
        ne = two_epoch_effective_size(1000, 10000, 50)
        assert ne == pytest.approx(10000 - 9000 * math.exp(-50 / 20000), rel=1e-12)
        assert 1000 < ne < 1100

    def test_recent_contraction_still_reflects_the_ancient_size(self):
        ne = two_epoch_effective_size(10000, 1000, 50)
        assert ne == pytest.approx(1000 + 9000 * math.exp(-50 / 2000), rel=1e-12)
        assert 9000 < ne < 10000

    def test_old_change_approaches_the_recent_size(self):
        assert two_epoch_effective_size(1000, 100, 10_000) == pytest.approx(100, rel=1e-6)

    def test_zero_time(self):
        assert two_epoch_effective_size(1000, 10000, 0) == 1000.0

    def test_matches_simulated_pairwise_coalescence(self):
        """E[T2] / 2 for a diploid pair: rate 1/(2 N_recent) until T, then 1/(2 N_ancient)."""
        ancient, recent, split = 1000, 100, 300
        rng = np.random.default_rng(20261026)
        first = rng.exponential(2 * recent, size=2_000_000)  # standard error ~0.15% of the mean
        coalescence = np.where(first < split, first, split + rng.exponential(2 * ancient, size=first.size))
        assert two_epoch_effective_size(ancient, recent, split) == pytest.approx(coalescence.mean() / 2, rel=0.01)


class TestExponentialGrowthModel:
    def test_single_generation(self):
        assert exponential_growth_model(10.0, 1.0, 1) == [10.0, 20.0]

    def test_compounding(self):
        sizes = exponential_growth_model(100.0, 0.1, 2)
        assert sizes[2] == pytest.approx(100.0 * 1.1**2)

    def test_zero_generations(self):
        assert exponential_growth_model(5.0, 0.2, 0) == [5.0]

    def test_validation(self):
        with pytest.raises(ValueError, match="positive"):
            exponential_growth_model(0.0, 0.1, 5)
        with pytest.raises(ValueError, match="negative"):
            exponential_growth_model(10.0, 0.1, -1)


class TestLogisticGrowthModel:
    def test_approaches_carrying_capacity(self):
        sizes = logistic_growth_model(10.0, 100.0, 0.3, 300)
        assert sizes[-1] < 100.0
        assert sizes[-1] > 99.0

    def test_starting_at_capacity_is_stationary(self):
        sizes = logistic_growth_model(100.0, 100.0, 0.5, 3)
        assert all(s == pytest.approx(100.0) for s in sizes)

    def test_validation(self):
        with pytest.raises(ValueError, match="Carrying capacity"):
            logistic_growth_model(10.0, 0.0, 0.1, 5)
        with pytest.raises(ValueError, match="negative"):
            logistic_growth_model(10.0, 100.0, 0.1, -1)


class TestAgeStructureModel:
    def test_single_generation_projection(self):
        fert = [0.0, 2.0]
        surv = [0.5, 0.0]
        res = age_structure_model(fert, surv, [10.0, 0.0], 1)
        # newborns = 10*fert[0] + 0*fert[1] = 0; age 1 = 10*surv[0] = 5
        assert res["population_history"][1] == [0.0, 5.0]
        assert res["total_populations"][1] == pytest.approx(5.0)

    def test_zero_generations_growth_rate(self):
        res = age_structure_model([1.0], [0.0], [10.0], 0)
        assert res["growth_rate"] == 0
        assert res["final_population"] == pytest.approx(10.0)

    def test_mismatched_lengths_raise(self):
        with pytest.raises(ValueError, match="same length"):
            age_structure_model([1.0, 2.0], [1.0, 1.0], [10.0], 1)


class TestIslandModelUpdate:
    def test_migration_pulls_toward_migrant_pool(self):
        assert island_model_update(0.2, 0.1, 0.4) == pytest.approx(0.22)

    def test_no_migration_is_identity(self):
        assert island_model_update(0.3, 0.0, 0.9) == pytest.approx(0.3)

    def test_validation(self):
        with pytest.raises(ValueError, match="between 0 and 1"):
            island_model_update(1.5, 0.1, 0.2)
        with pytest.raises(ValueError, match="between 0 and 1"):
            island_model_update(0.5, 1.2, 0.2)
        with pytest.raises(ValueError, match="between 0 and 1"):
            island_model_update(0.5, 0.1, -0.2)
