from __future__ import annotations

import pytest

from metainformant.math.population_genetics.effective_size import (
    effective_size_sex_ratio,
    harmonic_mean_effective_size,
)
from metainformant.math.population_genetics.statistics import effective_size_from_family_size_variance


def test_effective_size_from_family_size_variance():
    N = 1000.0
    Vk = 2.0
    Ne = effective_size_from_family_size_variance(N, Vk)
    assert abs(Ne - ((4.0 * N - 2.0) / (Vk + 2.0))) < 1e-12


def test_effective_size_from_family_size_list_mode():
    # mean k = 2, variance (population) = 0.5 -> Ne = 2 / (0.5 - 1) < 0? variance <= 1 -> inf
    assert effective_size_from_family_size_variance([2, 2, 2, 2]) == float("inf")
    # variance > 1: sizes [1, 5] -> mean 3, var 4 -> Ne = 3/3 = 1
    assert effective_size_from_family_size_variance([1, 5]) == pytest.approx(1.0)


def test_effective_size_from_family_size_empty_list():
    assert effective_size_from_family_size_variance([]) == 0.0


def test_sex_ratio_equal_reduces_to_wright_formula():
    assert effective_size_sex_ratio(100.0, 100.0) == pytest.approx(200.0)


def test_sex_ratio_unequal_uses_adjustment():
    val = effective_size_sex_ratio(100.0, 300.0, sex_ratio=0.25)
    harmonic = 2 * 100.0 * 300.0 / 400.0  # 150
    assert val == pytest.approx(harmonic / (0.25 * 0.75))


def test_sex_ratio_validation():
    with pytest.raises(ValueError, match="positive"):
        effective_size_sex_ratio(0.0, 100.0)
    with pytest.raises(ValueError, match="between 0 and 1"):
        effective_size_sex_ratio(100.0, 100.0, sex_ratio=1.0)


def test_harmonic_mean_effective_size_known_value():
    # harmonic mean of 100 and 200 = 2 / (1/100 + 1/200)
    val = harmonic_mean_effective_size([100.0, 200.0])
    assert val == pytest.approx(2.0 / (1.0 / 100.0 + 1.0 / 200.0))


def test_harmonic_mean_effective_size_validation():
    with pytest.raises(ValueError, match="empty"):
        harmonic_mean_effective_size([])
    with pytest.raises(ValueError, match="positive"):
        harmonic_mean_effective_size([100.0, 0.0])
