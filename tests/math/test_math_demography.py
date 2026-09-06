"""Tests for demographic models."""

from __future__ import annotations

import math

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
    """Test exponential growth effective size calculation."""

    def test_no_growth(self):
        """Test with zero growth rate."""
        ne = exponential_growth_effective_size(1000, growth_rate=0.0, generations=10)
        assert ne == 1000.0

    def test_positive_growth(self):
        """Test with positive growth rate."""
        ne = exponential_growth_effective_size(10000, growth_rate=0.23, generations=10)
        # Should be less than current size but greater than initial size
        assert ne < 10000
        assert ne > 1000  # Rough check

    def test_negative_growth(self):
        """Test with negative growth rate (decline)."""
        ne = exponential_growth_effective_size(1000, growth_rate=-0.1, generations=10)
        # Should be less than current size
        assert ne < 1000

    def test_zero_generations(self):
        """Test with zero generations."""
        ne = exponential_growth_effective_size(1000, growth_rate=0.1, generations=0)
        assert ne == 1000.0


class TestBottleneckEffectiveSize:
    """Test bottleneck effective size calculation."""

    def test_severe_bottleneck(self):
        """Test severe bottleneck."""
        ne = bottleneck_effective_size(10000, 100, 5)
        # Should be heavily weighted by bottleneck size
        assert ne < 1000  # Much less than pre-bottleneck
        assert ne > 50  # But more than just bottleneck (harmonic mean)

    def test_mild_bottleneck(self):
        """Test mild bottleneck."""
        ne = bottleneck_effective_size(10000, 5000, 5)
        # Should be intermediate
        assert ne < 10000
        assert ne > 1000

    def test_bottleneck_with_recovery(self):
        """Test bottleneck with recovery period."""
        ne_no_recovery = bottleneck_effective_size(10000, 100, 5)
        ne_with_recovery = bottleneck_effective_size(10000, 100, 5, recovery_generations=10)

        # With recovery, effective size should be larger
        assert ne_with_recovery > ne_no_recovery

    def test_zero_duration(self):
        """Test with zero bottleneck duration."""
        ne = bottleneck_effective_size(10000, 100, 0)
        assert ne == 10000.0


class TestTwoEpochEffectiveSize:
    """Test two-epoch effective size calculation."""

    def test_population_expansion(self):
        """Test population expansion."""
        ne = two_epoch_effective_size(1000, 10000, 50)
        # Should be intermediate (harmonic mean)
        assert ne > 1000
        assert ne < 10000

    def test_population_contraction(self):
        """Test population contraction."""
        ne = two_epoch_effective_size(10000, 1000, 50)
        # Should be intermediate (harmonic mean)
        assert ne > 1000
        assert ne < 10000

    def test_symmetric_epochs(self):
        """Test that expansion and contraction give same harmonic mean."""
        ne_expand = two_epoch_effective_size(1000, 10000, 50)
        ne_contract = two_epoch_effective_size(10000, 1000, 50)

        # Harmonic mean should be symmetric
        assert math.isclose(ne_expand, ne_contract, rel_tol=0.01)

    def test_zero_time(self):
        """Test with zero time since change."""
        ne = two_epoch_effective_size(1000, 10000, 0)
        assert ne == 10000.0


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
