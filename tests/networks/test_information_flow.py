"""Behavior tests for Granger causality and the OLS helper in information_flow."""

from __future__ import annotations

import random

import pytest

from metainformant.information.network_info import information_flow
from metainformant.information.network_info.information_flow import _ols_rss


def _ar1(length: int, phi: float, noise: float, rng: random.Random) -> list[float]:
    series = [0.0]
    for _ in range(length - 1):
        series.append(phi * series[-1] + rng.gauss(0.0, noise))
    return series


def _driven_ar1(
    source: list[float],
    phi: float,
    gain: float,
    noise: float,
    rng: random.Random,
) -> list[float]:
    """AR(1) target driven one step behind by the source series."""
    target = [0.0]
    for t in range(1, len(source)):
        target.append(
            phi * target[t - 1] + gain * source[t - 1] + rng.gauss(0.0, noise)
        )
    return target


class TestGrangerCausality:
    def test_causal_series_detected(self) -> None:
        rng = random.Random(12345)
        n = 400
        source = _ar1(n, phi=0.6, noise=1.0, rng=rng)
        target = _driven_ar1(source, phi=0.5, gain=0.8, noise=1.0, rng=rng)
        result = information_flow.granger_causality(source, target, max_lag=3)
        assert result["is_causal"] is True
        assert result["p_value"] < 0.05
        assert result["rss_unrestricted"] < result["rss_restricted"]
        assert result["f_statistic"] > 0.0
        assert 1 <= result["optimal_lag"] <= 3

    def test_null_series_not_causal(self) -> None:
        rng = random.Random(999)
        n = 400
        source = _ar1(n, phi=0.6, noise=1.0, rng=rng)
        target = _ar1(n, phi=0.5, noise=1.0, rng=rng)
        result = information_flow.granger_causality(source, target, max_lag=3)
        assert result["is_causal"] is False
        assert result["p_value"] > 0.05

    def test_reverse_direction_not_causal(self) -> None:
        # y is driven by x; testing y -> x must not flag causality.
        # n = 1200: at n = 400 this seed draws a finite-sample false
        # positive (p = 0.0117 at lag 1) even though the reverse direction
        # is Granger-null in population (p = 0.903 at n = 1200).
        rng = random.Random(4242)
        n = 1200
        source = _ar1(n, phi=0.6, noise=1.0, rng=rng)
        target = _driven_ar1(source, phi=0.5, gain=0.8, noise=1.0, rng=rng)
        result = information_flow.granger_causality(target, source, max_lag=3)
        assert result["is_causal"] is False
        assert result["p_value"] > 0.05

    def test_lag_selection_prefers_true_lag(self) -> None:
        # y_t depends only on x_{t-1}: the penalised residual-SS criterion
        # must select lag 1 rather than the maximum lag.
        rng = random.Random(777)
        n = 400
        source = _ar1(n, phi=0.5, noise=1.0, rng=rng)
        target = _driven_ar1(source, phi=0.0, gain=0.9, noise=0.5, rng=rng)
        result = information_flow.granger_causality(source, target, max_lag=5)
        assert result["optimal_lag"] == 1
        assert result["is_causal"] is True

    def test_invalid_max_lag_raises(self) -> None:
        with pytest.raises(ValueError, match="at least 1"):
            information_flow.granger_causality(
                [0.0, 0.1, 0.2, 0.3], [0.0, 0.1, 0.2, 0.3], max_lag=0
            )

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            information_flow.granger_causality([0.0, 1.0], [0.0, 1.0, 2.0])

    def test_too_short_series_raise(self) -> None:
        with pytest.raises(ValueError, match="too short"):
            information_flow.granger_causality([0.0, 0.1, 0.2], [0.0, 0.1, 0.2])


class TestOlsRss:
    def test_perfect_line_has_zero_rss(self) -> None:
        rss = _ols_rss([1.0, 2.0, 3.0, 4.0], [[1.0], [2.0], [3.0], [4.0]])
        assert rss < 1e-18

    def test_hand_computed_residual_ss(self) -> None:
        # y = (1, 2, 2), x = (1, 2, 3): centered OLS gives
        # xbar=2, ybar=5/3, Sxx=2, Sxy=1, slope=1/2,
        # rss = Syy - Sxy^2 / Sxx = 2/3 - 1/2 = 1/6.
        rss = _ols_rss([1.0, 2.0, 2.0], [[1.0], [2.0], [3.0]])
        assert rss == pytest.approx(2.0 / 3.0 - 1.0 / 2.0)

    def test_multiple_regressors(self) -> None:
        # y = 1 + 2*x1 + 3*x2 exactly -> zero RSS with two regressors.
        y = [1.0, 6.0, 4.0, 9.0]
        x = [[0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 2.0]]
        assert _ols_rss(y, x) < 1e-18

    def test_degenerate_constant_series(self) -> None:
        # k == 0 guard: no regressors, pure sum of squares.
        assert _ols_rss([1.0, -1.0, 2.0], []) == 6.0
