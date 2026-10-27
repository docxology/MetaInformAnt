"""Regression tests: the ``ld_matrix`` of ``compute_credible_set`` is a correlation matrix, and says so.

Defect: the docstring described ``ld_matrix`` as an "LD (r-squared) matrix", but the adjustment squares the
entries (``1 - r**2``), so a caller who followed the docstring and passed r-squared values was penalised with
r**4. The rest of the package (``conditional_analysis``, ``_apply_ld_adjustment``) treats the matrix as
correlations, and so do the existing LD tests; the docstring is what was wrong.

The numerical expectations are taken from the documented contract and computed here by hand: the log Bayes factor
of a variant is lowered by ``log(1 - r**2)``, with r its largest absolute correlation with any variant whose raw
Bayes factor is larger, and the factor ``1 - r**2`` is floored at 0.01.
"""

from __future__ import annotations

import inspect
import math

import pytest

from metainformant.gwas.finemapping.credible_sets import compute_credible_set

SE = 1.0  # standard errors of 1 put the Bayes factors on the z scale and keep the ranking by |z|


def _raw_log_bfs(z: list[float]) -> list[float]:
    result = compute_credible_set(z, standard_errors=SE)
    assert result["status"] == "success"
    return result["log_bayes_factors"]


def _adjusted_log_bfs(z: list[float], ld: list[list[float]]) -> list[float]:
    result = compute_credible_set(z, ld_matrix=ld, standard_errors=SE)
    assert result["status"] == "success"
    return result["log_bayes_factors"]


def _deltas(z: list[float], ld: list[list[float]]) -> list[float]:
    """Adjusted minus raw log Bayes factor per variant: the log of the LD penalty factor."""
    return [after - before for before, after in zip(_raw_log_bfs(z), _adjusted_log_bfs(z, ld))]


def _two_variant_ld(r: float) -> list[list[float]]:
    return [[1.0, r], [r, 1.0]]


def test_the_ld_penalty_is_one_minus_r_squared_of_the_correlation() -> None:
    deltas = _deltas([6.0, 5.0], _two_variant_ld(0.6))

    assert deltas[0] == pytest.approx(0.0, abs=1e-12)  # the strongest variant is never penalised
    assert deltas[1] == pytest.approx(math.log(1 - 0.6**2), abs=1e-12)


def test_an_r_squared_reading_of_the_matrix_is_not_what_the_code_does() -> None:
    """If entries were r-squared the factor would be 1 - 0.6 = 0.4; the code uses 1 - 0.6**2 = 0.64."""
    delta = _deltas([6.0, 5.0], _two_variant_ld(0.6))[1]

    assert delta != pytest.approx(math.log(1 - 0.6), abs=1e-3)
    assert math.exp(delta) == pytest.approx(0.64, abs=1e-12)


def test_the_sign_of_the_correlation_does_not_matter() -> None:
    assert _deltas([6.0, 5.0], _two_variant_ld(-0.6))[1] == pytest.approx(math.log(0.64), abs=1e-12)


def test_the_penalty_factor_is_floored_at_one_percent() -> None:
    delta = _deltas([6.0, 5.0], _two_variant_ld(0.999))[1]  # 1 - 0.999**2 = 0.002, below the floor

    assert delta == pytest.approx(math.log(0.01), abs=1e-12)


def test_the_strongest_correlation_with_a_larger_bayes_factor_decides() -> None:
    """Variant 2 correlates 0.3 with variant 0 and 0.8 with variant 1: the larger one, 0.8, is used.

    Variant 1 correlates 0.5 with variant 0. The correlation of variant 0 with the weaker variants is irrelevant.
    """
    ld = [
        [1.0, 0.5, 0.3],
        [0.5, 1.0, 0.8],
        [0.3, 0.8, 1.0],
    ]
    deltas = _deltas([6.0, 5.0, 4.0], ld)

    assert deltas[0] == pytest.approx(0.0, abs=1e-12)
    assert deltas[1] == pytest.approx(math.log(1 - 0.5**2), abs=1e-12)
    assert deltas[2] == pytest.approx(math.log(1 - 0.8**2), abs=1e-12)


def test_uncorrelated_variants_keep_their_bayes_factors() -> None:
    identity = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]

    assert _deltas([6.0, 5.0, 4.0], identity) == pytest.approx([0.0, 0.0, 0.0], abs=1e-12)


def test_pips_stay_a_distribution_after_the_adjustment() -> None:
    result = compute_credible_set(
        [6.0, 5.0, 4.0], ld_matrix=[[1, 0.5, 0.3], [0.5, 1, 0.8], [0.3, 0.8, 1]], standard_errors=SE
    )

    assert sum(result["pips"]) == pytest.approx(1.0, abs=1e-12)


def _ld_matrix_entry(function: object) -> str:
    """The ``ld_matrix`` entry of the Args section of a docstring (up to ``coverage:``), whitespace-normalised."""
    docstring = inspect.getdoc(function) or ""
    assert "ld_matrix:" in docstring
    return " ".join(docstring.split("ld_matrix:", 1)[1].split("coverage:", 1)[0].split())


def test_the_docstring_describes_a_correlation_matrix_not_r_squared() -> None:
    entry = _ld_matrix_entry(compute_credible_set).lower()

    assert "correlation" in entry
    assert "1 - r" in entry.replace("``", "")  # the documented adjustment is 1 - r**2
    # The only mention of r-squared is the warning that the matrix is not one.
    assert "r-squared" not in entry.replace("not an r-squared matrix", "")
