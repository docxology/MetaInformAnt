"""Tests for decision theory: DDM, island model.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.math.decision_theory.ddm import (
    ddm_analytic_accuracy,
    ddm_log_likelihood,
    ddm_mean_decision_time,
    fit_ddm_parameters,
    island_model_update,
)


# ---------------------------------------------------------------------------
# ddm_mean_decision_time
# ---------------------------------------------------------------------------


class TestDDMMeanDecisionTime:
    def test_basic_computation(self):
        t = ddm_mean_decision_time(drift_rate=1.0, boundary=2.0)
        assert t > 0

    def test_higher_drift_faster_decision(self):
        t_slow = ddm_mean_decision_time(drift_rate=0.5, boundary=2.0)
        t_fast = ddm_mean_decision_time(drift_rate=2.0, boundary=2.0)
        assert t_fast < t_slow

    def test_larger_boundary_slower_decision(self):
        t_small = ddm_mean_decision_time(drift_rate=1.0, boundary=1.0)
        t_large = ddm_mean_decision_time(drift_rate=1.0, boundary=4.0)
        assert t_large > t_small

    def test_zero_drift_returns_default(self):
        t = ddm_mean_decision_time(drift_rate=0.0, boundary=2.0)
        assert t == 1.0

    def test_negative_boundary_raises(self):
        with pytest.raises(ValueError):
            ddm_mean_decision_time(drift_rate=1.0, boundary=-1.0)


# ---------------------------------------------------------------------------
# ddm_analytic_accuracy
# ---------------------------------------------------------------------------


class TestDDMAnalyticAccuracy:
    def test_zero_drift_chance(self):
        acc = ddm_analytic_accuracy(drift_rate=0.0, boundary=2.0)
        assert acc == pytest.approx(0.5)

    def test_high_drift_high_accuracy(self):
        acc = ddm_analytic_accuracy(drift_rate=5.0, boundary=2.0)
        assert acc > 0.9

    def test_accuracy_between_0_and_1(self):
        acc = ddm_analytic_accuracy(drift_rate=1.0, boundary=1.0)
        assert 0.0 <= acc <= 1.0

    def test_negative_boundary_raises(self):
        with pytest.raises(ValueError):
            ddm_analytic_accuracy(drift_rate=1.0, boundary=-1.0)

    def test_increasing_drift_increases_accuracy(self):
        acc_low = ddm_analytic_accuracy(drift_rate=0.5, boundary=1.0)
        acc_high = ddm_analytic_accuracy(drift_rate=2.0, boundary=1.0)
        assert acc_high > acc_low


# ---------------------------------------------------------------------------
# ddm_log_likelihood
# ---------------------------------------------------------------------------


class TestDDMLogLikelihood:
    def test_valid_trial(self):
        ll = ddm_log_likelihood(rt=1.0, choice=1, drift_rate=1.0, boundary_separation=2.0)
        assert isinstance(ll, float)

    def test_rt_before_ndt_returns_neg_inf(self):
        ll = ddm_log_likelihood(rt=0.1, choice=1, drift_rate=1.0, boundary_separation=2.0, non_decision_time=0.5)
        assert ll == float("-inf")

    def test_invalid_choice_raises(self):
        with pytest.raises(ValueError, match="Choice must be 1 or 2"):
            ddm_log_likelihood(rt=1.0, choice=3, drift_rate=1.0, boundary_separation=2.0)

    def test_choice_2(self):
        ll = ddm_log_likelihood(rt=1.5, choice=2, drift_rate=1.0, boundary_separation=2.0)
        assert isinstance(ll, float)

    def test_different_starting_points(self):
        ll_half = ddm_log_likelihood(rt=1.0, choice=1, drift_rate=1.0, boundary_separation=2.0, starting_point=0.5)
        ll_high = ddm_log_likelihood(rt=1.0, choice=1, drift_rate=1.0, boundary_separation=2.0, starting_point=0.3)
        # Both should return finite or -inf values
        assert isinstance(ll_half, float)
        assert isinstance(ll_high, float)


# ---------------------------------------------------------------------------
# fit_ddm_parameters
# ---------------------------------------------------------------------------


class TestFitDDMParameters:
    def test_basic_fit(self):
        np.random.seed(42)
        rt_data = list(np.random.exponential(0.5, size=20) + 0.3)
        choice_data = [1 if np.random.random() > 0.3 else 2 for _ in range(20)]
        result = fit_ddm_parameters(rt_data, choice_data)
        assert "drift_rate" in result
        assert "boundary_separation" in result
        assert "starting_point" in result
        assert "non_decision_time" in result

    def test_mismatched_lengths_raises(self):
        with pytest.raises(ValueError, match="same length"):
            fit_ddm_parameters([1.0, 2.0], [1])

    def test_empty_data_raises(self):
        with pytest.raises(ValueError, match="No data"):
            fit_ddm_parameters([], [])

    def test_custom_bounds(self):
        np.random.seed(42)
        rt_data = list(np.random.exponential(0.5, size=15) + 0.3)
        choice_data = [1] * 10 + [2] * 5
        bounds = {
            "drift_rate": (0.1, 1.0),
            "boundary_separation": (0.5, 2.0),
            "starting_point": (0.3, 0.7),
            "non_decision_time": (0.0, 0.3),
        }
        result = fit_ddm_parameters(rt_data, choice_data, bounds=bounds)
        assert result["drift_rate"] >= 0.1
        assert result["boundary_separation"] >= 0.5


# ---------------------------------------------------------------------------
# island_model_update
# ---------------------------------------------------------------------------


class TestIslandModelUpdate:
    def test_no_migration(self):
        new_p = island_model_update(p=0.3, m=0.0, pm=0.8)
        assert new_p == pytest.approx(0.3)

    def test_full_migration(self):
        new_p = island_model_update(p=0.3, m=1.0, pm=0.8)
        assert new_p == pytest.approx(0.8)

    def test_partial_migration(self):
        new_p = island_model_update(p=0.2, m=0.5, pm=0.6)
        # (1-0.5)*0.2 + 0.5*0.6 = 0.1 + 0.3 = 0.4
        assert new_p == pytest.approx(0.4)

    def test_equal_frequencies_no_change(self):
        new_p = island_model_update(p=0.5, m=0.3, pm=0.5)
        assert new_p == pytest.approx(0.5)

    def test_frequency_bounds(self):
        new_p = island_model_update(p=0.0, m=0.1, pm=1.0)
        assert 0.0 <= new_p <= 1.0
