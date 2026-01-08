"""Tests for demographic models."""

from __future__ import annotations

import math

import pytest

from metainformant.math.population_genetics.demography import (
    bottleneck_effective_size,
    exponential_growth_effective_size,
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

