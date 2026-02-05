"""Morphometric profile analysis with allometry, indices, and size correction."""

from __future__ import annotations

import math
import statistics
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils.errors import ValidationError

from .measurement import Measurement


class MorphometricProfile:
    """A collection of morphometric measurements for a single specimen.

    Supports morphometric indices, allometric regressions,
    size correction, and coefficient of variation analysis.
    """

    def __init__(
        self,
        specimen_id: str,
        measurements: List[Measurement] | None = None,
        metadata: Dict[str, Any] | None = None,
    ) -> None:
        self.specimen_id = specimen_id
        self._measurements: Dict[str, Measurement] = {}
        if measurements:
            for m in measurements:
                self.add(m)
        self.metadata = metadata or {}

    @property
    def measurement_names(self) -> List[str]:
        return list(self._measurements.keys())

    @property
    def n_measurements(self) -> int:
        return len(self._measurements)

    def add(self, measurement: Measurement) -> None:
        """Add or replace a measurement."""
        self._measurements[measurement.name] = measurement

    def get(self, name: str) -> Optional[Measurement]:
        """Retrieve a measurement by name."""
        return self._measurements.get(name)

    def values_dict(self) -> Dict[str, float]:
        """Return all measurement values as a dict."""
        return {name: m.value for name, m in self._measurements.items()}

    def calculate_index(self, name: str, numerator: str, denominator: str) -> float:
        """Calculate a morphometric index (ratio * 100).

        Example: Cephalic Index (CI) = (Head Width / Head Length) * 100
        """
        m_num = self.get(numerator)
        m_den = self.get(denominator)

        if not m_num or not m_den:
            raise ValidationError(f"Missing measurements for index {name}: {numerator}, {denominator}")

        if m_num.unit != m_den.unit:
            try:
                m_num = m_num.convert(m_den.unit)
            except ValueError:
                raise ValidationError(f"Unit mismatch: {m_num.unit} vs {m_den.unit}")

        if m_den.value == 0:
            raise ValueError("Denominator is zero")

        return (m_num.value / m_den.value) * 100.0

    def geometric_mean_size(self) -> float:
        """Geometric mean of all measurements as a size proxy."""
        vals = [m.value for m in self._measurements.values() if m.value > 0]
        if not vals:
            return 0.0
        log_sum = sum(math.log(v) for v in vals)
        return math.exp(log_sum / len(vals))

    def size_corrected_values(self) -> Dict[str, float]:
        """Divide each measurement by the geometric mean size.

        Removes isometric size variation, leaving shape information.
        """
        gm = self.geometric_mean_size()
        if gm == 0:
            return {name: 0.0 for name in self._measurements}
        return {name: m.value / gm for name, m in self._measurements.items()}


def coefficient_of_variation(profiles: List[MorphometricProfile], measurement_name: str) -> Optional[float]:
    """CV (%) for a specific measurement across specimens.

    Args:
        profiles: List of MorphometricProfile objects.
        measurement_name: Name of the measurement to analyze.

    Returns:
        CV as percentage, or None if insufficient data.
    """
    values = []
    for p in profiles:
        m = p.get(measurement_name)
        if m is not None:
            values.append(m.value)

    if len(values) < 2:
        return None

    mean = statistics.mean(values)
    if mean == 0:
        return None

    sd = statistics.stdev(values)
    return (sd / mean) * 100.0


def allometric_regression(
    profiles: List[MorphometricProfile],
    x_measurement: str,
    y_measurement: str,
) -> Dict[str, float]:
    """Log-log linear regression for allometric scaling.

    Fits: log(y) = slope * log(x) + intercept
    Slope = 1.0 indicates isometry; >1 positive allometry; <1 negative allometry.

    Args:
        profiles: List of MorphometricProfile objects.
        x_measurement: Independent variable (e.g., body size).
        y_measurement: Dependent variable (e.g., head width).

    Returns:
        Dict with slope, intercept, r_squared, n, allometry_type.
    """
    log_x: List[float] = []
    log_y: List[float] = []

    for p in profiles:
        mx = p.get(x_measurement)
        my = p.get(y_measurement)
        if mx is not None and my is not None and mx.value > 0 and my.value > 0:
            log_x.append(math.log(mx.value))
            log_y.append(math.log(my.value))

    n = len(log_x)
    if n < 3:
        raise ValidationError(f"Need at least 3 valid pairs, got {n}")

    mean_x = sum(log_x) / n
    mean_y = sum(log_y) / n

    ss_xy = sum((log_x[i] - mean_x) * (log_y[i] - mean_y) for i in range(n))
    ss_xx = sum((log_x[i] - mean_x) ** 2 for i in range(n))
    ss_yy = sum((log_y[i] - mean_y) ** 2 for i in range(n))

    if ss_xx == 0:
        raise ValidationError("No variance in x measurement")

    slope = ss_xy / ss_xx
    intercept = mean_y - slope * mean_x

    r_squared = (ss_xy**2) / (ss_xx * ss_yy) if ss_yy > 0 else 0.0

    if abs(slope - 1.0) < 0.05:
        allometry_type = "isometric"
    elif slope > 1.0:
        allometry_type = "positive"
    else:
        allometry_type = "negative"

    return {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_squared,
        "n": n,
        "allometry_type": allometry_type,
    }


def compare_profiles(
    profile_a: MorphometricProfile,
    profile_b: MorphometricProfile,
) -> Dict[str, Any]:
    """Compare two morphometric profiles measurement by measurement.

    Returns:
        Dict with shared measurements, differences, and percent differences.
    """
    names_a = set(profile_a.measurement_names)
    names_b = set(profile_b.measurement_names)
    shared = names_a & names_b

    comparisons = {}
    for name in sorted(shared):
        ma = profile_a.get(name)
        mb = profile_b.get(name)
        if ma is not None and mb is not None:
            # Ensure same units
            if ma.unit != mb.unit:
                try:
                    mb = mb.convert(ma.unit)
                except ValueError:
                    continue

            diff = ma.value - mb.value
            mean_val = (ma.value + mb.value) / 2
            pct_diff = (diff / mean_val * 100) if mean_val != 0 else 0.0

            comparisons[name] = {
                "value_a": ma.value,
                "value_b": mb.value,
                "difference": diff,
                "percent_difference": pct_diff,
                "unit": ma.unit,
            }

    return {
        "specimen_a": profile_a.specimen_id,
        "specimen_b": profile_b.specimen_id,
        "shared_measurements": len(shared),
        "unique_to_a": sorted(names_a - names_b),
        "unique_to_b": sorted(names_b - names_a),
        "comparisons": comparisons,
    }


def summary_statistics(
    profiles: List[MorphometricProfile],
) -> Dict[str, Dict[str, float]]:
    """Summary statistics for each measurement across all profiles.

    Returns:
        Dict mapping measurement name to {mean, median, std, min, max, n, cv}.
    """
    # Collect all values per measurement
    measurement_values: Dict[str, List[float]] = {}
    for p in profiles:
        for name in p.measurement_names:
            m = p.get(name)
            if m is not None:
                measurement_values.setdefault(name, []).append(m.value)

    result = {}
    for name, vals in sorted(measurement_values.items()):
        if not vals:
            continue
        mean = statistics.mean(vals)
        std = statistics.stdev(vals) if len(vals) > 1 else 0.0
        result[name] = {
            "mean": mean,
            "median": statistics.median(vals),
            "std": std,
            "min": min(vals),
            "max": max(vals),
            "n": len(vals),
            "cv": (std / mean * 100) if mean != 0 else 0.0,
        }
    return result
