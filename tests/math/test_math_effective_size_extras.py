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
    # Ne = (4N - 2) / (Vk + 2) with N = len(list) and Vk the population variance of the counts.
    # [2, 2, 2, 2]: N = 4, Vk = 0 -> 14 / 2 = 7 (equal family sizes give Ne = 2N - 1)
    assert effective_size_from_family_size_variance([2, 2, 2, 2]) == pytest.approx(7.0)
    # [0, 1, 2, 3, 4]: N = 5, deviations -2..2 give a sum of squares of 10, Vk = 10 / 5 = 2 -> 18 / 4 = 4.5
    assert effective_size_from_family_size_variance([0, 1, 2, 3, 4]) == pytest.approx(4.5)
    # [0, 4]: N = 2, Vk = 4 -> 6 / 6 = 1
    assert effective_size_from_family_size_variance([0, 4]) == pytest.approx(1.0)


def test_effective_size_from_family_size_empty_list():
    assert effective_size_from_family_size_variance([]) == 0.0


def test_sex_ratio_equal_reduces_to_wright_formula():
    assert effective_size_sex_ratio(100.0, 100.0) == pytest.approx(200.0)


def test_sex_ratio_unequal_sexes_follow_wright_formula():
    # The breeding sex ratio (here 0.25) is fixed by the two counts, so Wright's formula
    # 4 Nm Nf / (Nm + Nf) applies as it stands: 4 * 100 * 300 / 400 = 300, below the census size 400.
    val = effective_size_sex_ratio(100.0, 300.0, sex_ratio=0.25)
    assert val == pytest.approx(4 * 100.0 * 300.0 / (100.0 + 300.0))
    assert val < 100.0 + 300.0


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
