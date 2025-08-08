from __future__ import annotations

from metainformant.math import expected_time_to_mrca, expected_total_branch_length


def test_expected_time_to_mrca_simple():
    # n=2: sum_{k=2}^2 1/(k(k-1)) = 1/2 => 4N * 1/2 = 2N
    assert abs(expected_time_to_mrca(2, 1000) - 2000.0) < 1e-9


def test_expected_total_branch_length_harmonic():
    # n=4: H_{3} = 1 + 1/2 + 1/3 = 1.833333...
    val = expected_total_branch_length(4, 100)
    assert abs(val - (4 * 100 * (1.0 + 0.5 + 1.0 / 3.0))) < 1e-9



