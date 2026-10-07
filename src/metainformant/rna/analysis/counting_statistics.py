"""Counting intervals for descriptive campaign reports, not biological inference."""

from __future__ import annotations

import math
from numbers import Integral


def wilson_interval(successes: object, total: object, z: float = 1.959963984540054) -> tuple[float, float] | None:
    """Return a Wilson score interval for a proportion from counted evidence.

    This quantifies only the counting uncertainty of an observed denominator in
    one data-root snapshot. It is not an inferential or biological uncertainty.
    Returns ``None`` when the denominator is unavailable or zero.
    """

    if (
        isinstance(successes, bool)
        or isinstance(total, bool)
        or not isinstance(successes, Integral)
        or not isinstance(total, Integral)
    ):
        raise ValueError("Wilson interval requires integer counts")
    if not math.isfinite(z) or z <= 0:
        raise ValueError("Wilson interval requires a finite positive z score")
    successes, total = int(successes), int(total)
    if total <= 0 or successes < 0 or successes > total:
        return None
    p = successes / total
    z2 = z * z
    denom = 1.0 + z2 / total
    center = (p + z2 / (2.0 * total)) / denom
    half = z * math.sqrt(p * (1.0 - p) / total + z2 / (4.0 * total * total)) / denom
    return max(0.0, center - half), min(1.0, center + half)
