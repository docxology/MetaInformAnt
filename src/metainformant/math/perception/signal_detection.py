"""Signal Detection Theory (SDT) metrics.

This module implements key metrics for analyzing decision making under uncertainty
using Signal Detection Theory (Sensitivity d', Bias).
"""

from __future__ import annotations

import math

import numpy as np


def _check_rate(name: str, rate: float) -> None:
    """Reject a rate that is not a proportion; only exact 0 and 1 are clipped, never an out-of-range value."""
    if not 0.0 <= rate <= 1.0:
        raise ValueError(f"{name} must be a proportion between 0 and 1, got {rate}")


def d_prime(hit_rate: float, false_alarm_rate: float, correction: float = 1e-6) -> float:
    """Calculate sensitivity index d' (d-prime).

    d' = Z(Hit Rate) - Z(False Alarm Rate)

    Args:
        hit_rate: Proportion of signal trials correctly identified (0 to 1)
        false_alarm_rate: Proportion of noise trials incorrectly identified as signal
        correction: Small epsilon to avoid infinite Z-scores at rates of exactly 0 or 1

    Returns:
        d' value. Higher values indicate better sensitivity.
        0 indicates inability to distinguish signal from noise.

    Raises:
        ValueError: If a rate is outside [0, 1] (NaN included)
    """
    from scipy.stats import norm

    _check_rate("hit_rate", hit_rate)
    _check_rate("false_alarm_rate", false_alarm_rate)

    # Clip rates to avoid infinity
    h = np.clip(hit_rate, correction, 1 - correction)
    fa = np.clip(false_alarm_rate, correction, 1 - correction)

    return float(norm.ppf(h)) - float(norm.ppf(fa))


def criterion_c(hit_rate: float, false_alarm_rate: float, correction: float = 1e-6) -> float:
    """Calculate decision criterion c (bias).

    c = -(Z(Hit Rate) + Z(False Alarm Rate)) / 2

    Args:
        hit_rate: Proportion of signal trials correctly identified
        false_alarm_rate: Proportion of noise trials incorrectly identified
        correction: Small epsilon

    Returns:
        c value.
        c = 0: Unbiased
        c > 0: Conservative (tendency to say "No")
        c < 0: Liberal (tendency to say "Yes")

    Raises:
        ValueError: If a rate is outside [0, 1] (NaN included)
    """
    from scipy.stats import norm

    _check_rate("hit_rate", hit_rate)
    _check_rate("false_alarm_rate", false_alarm_rate)

    h = np.clip(hit_rate, correction, 1 - correction)
    fa = np.clip(false_alarm_rate, correction, 1 - correction)

    return -(float(norm.ppf(h)) + float(norm.ppf(fa))) / 2.0


def likelihood_ratio_beta(hit_rate: float, false_alarm_rate: float, correction: float = 1e-6) -> float:
    """Calculate likelihood ratio β (beta).

    β = exp(d' * c)

    Args:
        hit_rate: Hit rate
        false_alarm_rate: False alarm rate
        correction: Small epsilon

    Returns:
        β value.
        β = 1: Unbiased
        β > 1: Conservative
        β < 1: Liberal
    """
    d = d_prime(hit_rate, false_alarm_rate, correction)
    c = criterion_c(hit_rate, false_alarm_rate, correction)

    return math.exp(d * c)


def sdt_metrics(hits: int, misses: int, false_alarms: int, correct_rejections: int) -> dict[str, float]:
    """Calculate standard SDT metrics from raw counts.

    Args:
        hits: Signal present, response "Yes"
        misses: Signal present, response "No"
        false_alarms: Signal absent, response "Yes"
        correct_rejections: Signal absent, response "No"

    Returns:
        Dictionary containing hit_rate, false_alarm_rate, d_prime, c, beta

    Raises:
        ValueError: If a count is negative, or there are no signal trials (hits + misses) or no noise
            trials (false alarms + correct rejections): the hit or false-alarm rate is then undefined,
            and no rate of 0 is invented for it
    """
    if min(hits, misses, false_alarms, correct_rejections) < 0:
        raise ValueError("Trial counts must be non-negative")
    n_signal = hits + misses
    n_noise = false_alarms + correct_rejections
    if n_signal == 0:
        raise ValueError("No signal trials (hits + misses = 0): the hit rate is undefined")
    if n_noise == 0:
        raise ValueError("No noise trials (false alarms + correct rejections = 0): the false-alarm rate is undefined")

    hit_rate = hits / n_signal
    fa_rate = false_alarms / n_noise

    return {
        "hit_rate": hit_rate,
        "false_alarm_rate": fa_rate,
        "d_prime": d_prime(hit_rate, fa_rate),
        "c": criterion_c(hit_rate, fa_rate),
        "beta": likelihood_ratio_beta(hit_rate, fa_rate),
    }
