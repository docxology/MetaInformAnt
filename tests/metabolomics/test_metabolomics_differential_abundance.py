"""Tests for metabolomics differential abundance (Welch t-test statistics)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.metabolomics.analysis.identification import differential_abundance

scipy_stats = pytest.importorskip("scipy.stats")


def _welch_reference(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Independent scipy Welch reference for one metabolite row."""
    va = a.var(ddof=1) / a.size
    vb = b.var(ddof=1) / b.size
    se2 = va + vb
    t_ref = (a.mean() - b.mean()) / math.sqrt(se2)
    df_ref = se2**2 / (va**2 / (a.size - 1) + vb**2 / (b.size - 1))
    p_ref = 2.0 * scipy_stats.t.sf(abs(t_ref), df_ref)
    return t_ref, p_ref


def test_n5_vs_n4_matches_scipy_welch() -> None:
    """Unequal group sizes with unequal variances must follow Welch-Satterthwaite."""
    rng = np.random.default_rng(3)
    a = rng.normal(10.0, 2.0, (1, 5))
    b = rng.normal(11.0, 3.0, (1, 4))
    data = np.hstack([a, b])

    t_stats, p_values = differential_abundance(data, [0, 1, 2, 3, 4], [5, 6, 7, 8])

    t_ref, p_ref = _welch_reference(a[0], b[0])
    assert t_stats[0] == pytest.approx(t_ref, rel=1e-12, abs=1e-12)
    assert p_values[0] == pytest.approx(p_ref, rel=1e-12, abs=1e-12)


def test_n40_vs_n40_matches_scipy_welch() -> None:
    """Balanced large groups must match the scipy reference exactly."""
    rng = np.random.default_rng(4)
    a = rng.normal(0.0, 1.0, (1, 40))
    b = rng.normal(0.5, 1.0, (1, 40))
    data = np.hstack([a, b])

    t_stats, p_values = differential_abundance(
        data, list(range(40)), list(range(40, 80))
    )

    t_ref, p_ref = _welch_reference(a[0], b[0])
    assert t_stats[0] == pytest.approx(t_ref, rel=1e-12, abs=1e-12)
    assert p_values[0] == pytest.approx(p_ref, rel=1e-12, abs=1e-12)


def test_small_df_deviates_from_normal_approximation() -> None:
    """With tiny df the t distribution must give heavier tails than a normal."""
    rng = np.random.default_rng(5)
    a = rng.normal(0.0, 0.5, (1, 3))
    b = rng.normal(2.0, 2.0, (1, 3))
    data = np.hstack([a, b])

    t_stats, p_values = differential_abundance(data, [0, 1, 2], [3, 4, 5])

    t_ref, _ = _welch_reference(a[0], b[0])
    p_normal = 2.0 * scipy_stats.norm.sf(abs(t_ref))
    assert t_stats[0] == pytest.approx(t_ref, rel=1e-12, abs=1e-12)
    assert p_values[0] > p_normal + 1e-6


def test_zero_variance_groups() -> None:
    """Degenerate groups must not fabricate moderate p-values."""
    data = np.vstack([np.full(8, 5.0), np.full(8, 5.0)])
    t_stats, p_values = differential_abundance(data, [0, 1, 2, 3], [4, 5, 6, 7])
    assert t_stats[0] == 0.0
    assert p_values[0] == 1.0

    # One metabolite: four constant low samples vs four constant high samples.
    data_diff = np.array([[5.0] * 4 + [7.0] * 4])
    t_stats, p_values = differential_abundance(data_diff, [0, 1, 2, 3], [4, 5, 6, 7])
    assert abs(t_stats[0]) == float("inf")
    assert p_values[0] == 0.0


def test_identical_distributions_no_signal() -> None:
    """Halves of one Gaussian sample share a distribution: no signal may appear."""
    rng = np.random.default_rng(6)
    data = rng.normal(3.0, 1.0, (1, 10))

    t_stats, p_values = differential_abundance(data, [0, 1, 2, 3, 4], [5, 6, 7, 8, 9])

    assert abs(t_stats[0]) < 2.0
    assert p_values[0] > 0.05
