"""Peak areas must not depend on a NumPy 2-only function.

Regression tests for ``metainformant.metabolomics.visualization.plots``. ``detect_chromatographic_peaks`` integrated
each peak with ``np.trapezoid``, which exists only from NumPy 2.0, while ``pyproject.toml`` allows ``numpy>=1.24``:
on NumPy 1.x every call that found a peak raised ``AttributeError``. (``np.trapz``, the pre-2.0 name, is gone from
current NumPy, so neither name works everywhere.)

The behavioural tests check the areas against hand-computed and analytic values. The source guard exists because the
runtime failure cannot be reproduced on the NumPy that the test environment has: it fails on the old module whatever
NumPy is installed.

Real computation; no test doubles.
"""

from __future__ import annotations

import ast
import math
from pathlib import Path

import numpy as np
import pytest

from metainformant.metabolomics import visualization as visualization_package
from metainformant.metabolomics.visualization import plots
from metainformant.metabolomics.visualization.plots import detect_chromatographic_peaks

# attribute names that exist in only some of the NumPy releases the package supports (numpy>=1.24)
NUMPY_VERSION_SPECIFIC_NAMES = {"trapezoid"}


def _numpy_trapezoid_reference():
    """NumPy's own trapezoidal rule, under whichever name this NumPy has (an independent oracle)."""
    function = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
    if function is None:  # pragma: no cover - NumPy always has one of the two
        pytest.skip("this NumPy has neither np.trapezoid nor np.trapz")
    return function


class TestPeakAreas:
    def test_a_triangular_peak_has_the_hand_computed_area(self) -> None:
        x = np.arange(7, dtype=float)
        y = np.array([0, 10, 50, 100, 50, 10, 0], dtype=float)

        (peak,) = detect_chromatographic_peaks(x, y)

        # (0+10)/2 + (10+50)/2 + (50+100)/2 + (100+50)/2 + (50+10)/2 + (10+0)/2 = 5 + 30 + 75 + 75 + 30 + 5
        assert (peak.left_idx, peak.peak_idx, peak.right_idx) == (0, 3, 6)
        assert peak.area == 220.0

    def test_unevenly_spaced_retention_times_are_weighted_by_their_gaps(self) -> None:
        x = np.array([0.0, 1.0, 3.0, 6.0])
        y = np.array([1.0, 2.0, 4.0, 8.0])

        # (1+2)/2 * 1 + (2+4)/2 * 2 + (4+8)/2 * 3
        assert plots._trapezoid_area(y, x) == 25.5

    def test_a_gaussian_peak_has_the_analytic_area(self) -> None:
        height, centre, sigma = 1000.0, 50.0, 5.0
        x = np.linspace(0.0, 100.0, 2001)
        y = height * np.exp(-((x - centre) ** 2) / (2 * sigma**2))

        (peak,) = detect_chromatographic_peaks(x, y, min_intensity=1.0)

        assert peak.area == pytest.approx(height * sigma * math.sqrt(2 * math.pi), rel=1e-6)

    def test_the_helper_is_bitwise_equal_to_numpys_own_trapezoid(self) -> None:
        rng = np.random.default_rng(5)
        x = np.sort(rng.uniform(0, 100, size=57))
        y = rng.uniform(0, 1e4, size=57)

        assert plots._trapezoid_area(y, x) == float(_numpy_trapezoid_reference()(y, x))


class TestNumpyVersionCompatibility:
    def test_the_metabolomics_package_calls_no_numpy_function_that_is_missing_from_supported_versions(self) -> None:
        package_root = Path(visualization_package.__file__).resolve().parents[1]
        offences: list[str] = []
        for source in sorted(package_root.rglob("*.py")):
            for node in ast.walk(ast.parse(source.read_text(), filename=str(source))):
                if (
                    isinstance(node, ast.Attribute)
                    and node.attr in NUMPY_VERSION_SPECIFIC_NAMES
                    and isinstance(node.value, ast.Name)
                    and node.value.id in {"np", "numpy"}
                ):
                    offences.append(f"{source.relative_to(package_root)}:{node.lineno} np.{node.attr}")

        assert offences == [], f"NumPy-2-only calls while numpy>=1.24 is supported: {offences}"
