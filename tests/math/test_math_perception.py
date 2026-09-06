"""Tests for perception: psychophysical laws and signal detection theory metrics."""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.math.perception.psychophysics import (
    fechner_law,
    michelson_contrast,
    stevens_power_law,
    weber_contrast,
)
from metainformant.math.perception.signal_detection import (
    criterion_c,
    d_prime,
    likelihood_ratio_beta,
    sdt_metrics,
)


class TestWeberContrast:
    def test_positive_contrast(self):
        assert weber_contrast(0.6, 0.4) == pytest.approx(0.5)

    def test_zero_background_nonzero_intensity_is_inf(self):
        assert weber_contrast(1.0, 0.0) == float("inf")

    def test_zero_background_zero_intensity(self):
        assert weber_contrast(0.0, 0.0) == 0.0

    def test_equal_intensities_is_zero(self):
        assert weber_contrast(0.5, 0.5) == pytest.approx(0.0)


class TestMichelsonContrast:
    def test_known_value(self):
        # (max - min) / (max + min)
        assert michelson_contrast(1.0, 0.5) == pytest.approx(1.0 / 3.0)

    def test_zero_denominator(self):
        assert michelson_contrast(0.0, 0.0) == 0.0

    def test_full_contrast(self):
        assert michelson_contrast(1.0, 0.0) == pytest.approx(1.0)


class TestFechnerLaw:
    def test_known_value(self):
        # S = k * ln(I / I0)
        assert fechner_law(math.e, threshold=1.0, k=2.0) == pytest.approx(2.0)

    def test_at_threshold_is_zero(self):
        assert fechner_law(1.0, threshold=1.0) == pytest.approx(0.0)

    def test_array_input(self):
        result = fechner_law(np.array([1.0, math.e]), threshold=1.0, k=1.0)
        assert isinstance(result, np.ndarray)
        assert result[1] == pytest.approx(1.0)


class TestStevensPowerLaw:
    def test_known_value(self):
        assert stevens_power_law(4.0, exponent=0.5, k=1.0) == pytest.approx(2.0)

    def test_scaling_constant(self):
        assert stevens_power_law(2.0, exponent=2.0, k=3.0) == pytest.approx(12.0)

    def test_array_input(self):
        result = stevens_power_law(np.array([1.0, 4.0]), exponent=0.5)
        assert result[1] == pytest.approx(2.0)


class TestDPrime:
    def test_equal_rates_zero_sensitivity(self):
        assert d_prime(0.5, 0.5) == pytest.approx(0.0)

    def test_positive_sensitivity(self):
        # Hits more frequent than false alarms -> d' > 0
        assert d_prime(0.8, 0.2) > 0.0

    def test_degenerate_rates_are_clipped_finite(self):
        # Rates of exactly 0/1 would give infinite z; clipping keeps d' finite
        assert math.isfinite(d_prime(1.0, 0.0))

    def test_symmetry(self):
        assert d_prime(0.8, 0.2) == pytest.approx(-d_prime(0.2, 0.8))


class TestCriterionC:
    def test_unbiased_at_half_rates(self):
        assert criterion_c(0.5, 0.5) == pytest.approx(0.0)

    def test_conservative_bias_positive(self):
        # Low false-alarm rate -> liberal/conservative sign convention: c > 0
        assert criterion_c(0.8, 0.1) > 0.0

    def test_degenerate_rates_are_clipped_finite(self):
        assert math.isfinite(criterion_c(1.0, 0.0))


class TestLikelihoodRatioBeta:
    def test_unbiased_beta_is_one(self):
        assert likelihood_ratio_beta(0.5, 0.5) == pytest.approx(1.0)

    def test_positive_for_valid_rates(self):
        assert likelihood_ratio_beta(0.8, 0.2) > 0.0


class TestSDTMetrics:
    def test_basic_counts(self):
        result = sdt_metrics(hits=80, misses=20, false_alarms=10, correct_rejections=90)
        assert result["hit_rate"] == pytest.approx(0.8)
        assert result["false_alarm_rate"] == pytest.approx(0.1)
        assert 0.0 <= result["d_prime"] <= 6.0
        assert result["beta"] > 0.0

    def test_no_signal_trials(self):
        # Empty noise trial class must not crash
        result = sdt_metrics(hits=5, misses=5, false_alarms=0, correct_rejections=0)
        assert result["false_alarm_rate"] == pytest.approx(0.0)

    def test_no_noise_trials(self):
        result = sdt_metrics(hits=0, misses=0, false_alarms=3, correct_rejections=7)
        assert result["hit_rate"] == pytest.approx(0.0)
