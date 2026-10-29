"""The k-mer enrichment p-value of ``motif_analysis`` is the exact binomial tail, not a normal approximation.

``_binomial_pvalue`` used the normal approximation ``0.5 * erfc(z / sqrt(2))``. For the rare events k-mer enrichment
tests (a 6-mer has an expected frequency of 1/4096) the normal tail is far too small: 8 occurrences among 10 000
k-mers got 1.9e-4 where the exact binomial tail is 3.7e-3, so k-mers were ranked, and reported in ``p_value``, with
values about twenty times too small. The oracle is ``scipy.stats.binom.sf``. No test doubles.
"""

from __future__ import annotations

import pytest

stats = pytest.importorskip("scipy.stats")

from metainformant.networks.regulatory.motif_analysis import _binomial_pvalue  # noqa: E402


@pytest.mark.parametrize(
    ("n", "p"),
    [(10000, 1 / 4096), (5000, 1 / 256), (500, 0.001), (50, 0.3), (100000, 1 / 1024), (12, 0.5)],
)
def test_the_tail_matches_scipy(n: int, p: float) -> None:
    for k in (1, 2, 3, 5, 8, 12, 20, 60, 400):
        if k > n:
            continue
        reference = float(stats.binom.sf(k - 1, n, p))
        assert _binomial_pvalue(k, n, p) == pytest.approx(reference, rel=1e-8, abs=1e-300), (k, n, p)


def test_the_rare_event_example() -> None:
    exact = float(stats.binom.sf(7, 10000, 1 / 4096))

    assert exact == pytest.approx(3.689e-3, rel=1e-3, abs=1e-300)
    assert _binomial_pvalue(8, 10000, 1 / 4096) == pytest.approx(exact, rel=1e-8, abs=1e-300)


def test_the_tail_is_a_survival_function() -> None:
    values = [_binomial_pvalue(k, 300, 0.02) for k in range(0, 40)]

    assert values[0] == 1.0
    assert all(later <= earlier + 1e-15 for earlier, later in zip(values, values[1:]))
    assert all(0.0 <= value <= 1.0 for value in values)


def test_degenerate_arguments_keep_their_conventions() -> None:
    assert _binomial_pvalue(3, 0, 0.1) == 1.0  # nothing to test
    assert _binomial_pvalue(3, 10, 0.0) == 1.0 and _binomial_pvalue(3, 10, 1.0) == 1.0  # no usable expectation
    assert _binomial_pvalue(0, 10, 0.2) == 1.0
    assert _binomial_pvalue(11, 10, 0.2) == 0.0
