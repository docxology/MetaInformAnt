from __future__ import annotations

from metainformant.math import expected_segregating_sites


def test_expected_segregating_sites():
    theta = 0.01
    n = 5
    # a1 = 1 + 1/2 + 1/3 + 1/4
    a1 = 1 + 0.5 + 1.0 / 3.0 + 0.25
    E_S = expected_segregating_sites(theta, n)
    assert abs(E_S - theta * a1) < 1e-12
