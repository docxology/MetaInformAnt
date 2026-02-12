from __future__ import annotations

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
