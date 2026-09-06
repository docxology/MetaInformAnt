from __future__ import annotations

import pytest

from metainformant.math.population_genetics.ld import haldane_c_to_d, haldane_d_to_c, kosambi_c_to_d, kosambi_d_to_c
from metainformant.math.population_genetics.statistics import expected_r2_from_Ne_c


def test_haldane_and_kosambi_mapping_functions():
    # Haldane: c = 0.5(1 - exp(-2d)) and inverse d = -0.5 ln(1 - 2c)
    d = 0.01  # 1 cM = 0.01 Morgans
    c = haldane_d_to_c(d)
    d_back = haldane_c_to_d(c)
    assert abs(d - d_back) < 1e-12

    # Kosambi: d = 0.25 ln((1+2c)/(1-2c)) and inverse c = 0.5 tanh(2d)
    d = 0.02
    c = kosambi_d_to_c(d)
    d_back = kosambi_c_to_d(c)
    assert abs(d - d_back) < 1e-12


def test_expected_r2_from_Ne_c():
    Ne = 1000.0
    c = 0.01
    r2 = expected_r2_from_Ne_c(Ne, c)
    assert abs(r2 - (1.0 / (1.0 + 4.0 * Ne * c))) < 1e-18


def test_haldane_mapping_boundaries():
    assert haldane_c_to_d(0.0) == 0.0
    assert haldane_c_to_d(0.5) == float("inf")
    assert haldane_d_to_c(0.0) == 0.0
    with pytest.raises(ValueError, match="between 0 and 0.5"):
        haldane_c_to_d(0.6)


def test_kosambi_mapping_boundaries():
    assert kosambi_c_to_d(0.0) == 0.0
    assert kosambi_c_to_d(0.5) == float("inf")
    assert kosambi_d_to_c(0.0) == 0.0
    with pytest.raises(ValueError, match="between 0 and 0.5"):
        kosambi_c_to_d(-0.1)
    with pytest.raises(ValueError, match="negative"):
        kosambi_d_to_c(-1.0)


def test_expected_r2_from_Ne_c_rejects_missing_Ne():
    with pytest.raises(ValueError, match="Ne must be provided"):
        expected_r2_from_Ne_c(0.01, None)
