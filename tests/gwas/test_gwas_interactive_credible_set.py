"""The plotting credible set uses the Wakefield implementation, not ABF = 1/p.

Zero-mocks: ``finemapping.credible_sets.compute_credible_set`` on the same
Z-scores and standard errors is the oracle.
"""

from __future__ import annotations

import pytest
from scipy import stats

from metainformant.gwas.finemapping import credible_sets
from metainformant.gwas.visualization.interactive.finemapping import compute_credible_set


def _results(betas, ses):
    return [
        {"beta": b, "se": s, "p_value": float(2 * stats.norm.sf(abs(b / s)))} for b, s in zip(betas, ses, strict=True)
    ]


def test_beta_and_se_use_the_wakefield_variance_model() -> None:
    betas, ses = [0.45, 0.02, -0.01, 0.03], [0.06, 0.05, 0.05, 0.07]
    result = compute_credible_set(_results(betas, ses))
    oracle = credible_sets.compute_credible_set([b / s for b, s in zip(betas, ses)], standard_errors=ses)

    assert result["pips"] == pytest.approx(oracle["pips"], rel=1e-12)
    assert result["credible_set_indices"][0] == 0
    assert result["credible_set_size"] == oracle["n_snps"]


def test_p_values_alone_fall_back_to_the_z_scale() -> None:
    p_values = [1e-12, 0.4, 0.7, 0.2]
    result = compute_credible_set([{"p_value": p} for p in p_values])
    z_scores = [float(stats.norm.isf(p / 2)) for p in p_values]
    oracle = credible_sets.compute_credible_set(z_scores)

    assert result["pips"] == pytest.approx(oracle["pips"], rel=1e-12)
    assert result["credible_set_indices"] == [0]


def test_extreme_z_scores_keep_nonzero_p_values() -> None:
    assert credible_sets._z_to_p(10.0) == pytest.approx(2 * stats.norm.sf(10.0), rel=1e-12)
    assert credible_sets._z_to_p(30.0) > 0
