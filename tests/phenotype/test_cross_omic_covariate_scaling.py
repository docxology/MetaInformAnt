"""Covariate units must not change the genotype standard error or significance."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats

from metainformant.phenotype.integration.cross_omic import _covariate_adjusted_fit


@pytest.mark.parametrize("scale", [1e-8, 1e8])
def test_covariate_units_preserve_genotype_inference(scale: float) -> None:
    rng = np.random.default_rng(421)
    n = 80
    genotype = rng.integers(0, 3, n).astype(float)
    covariate = rng.normal(size=n)
    target = 0.15 * genotype + 0.8 * covariate + rng.normal(size=n)
    design = np.column_stack([np.ones(n), covariate, genotype])
    q, r = np.linalg.qr(design, mode="reduced")
    coefficients = np.linalg.solve(r, q.T @ target)
    residual = target - design @ coefficients
    residual_sd = np.sqrt(float(residual @ residual) / (n - 3))
    inverse_r = np.linalg.solve(r, np.eye(3))
    expected_se = residual_sd * np.linalg.norm(inverse_r[-1])
    expected_p = 2.0 * stats.t.sf(abs(coefficients[-1] / expected_se), n - 3)

    fitted = _covariate_adjusted_fit(target.tolist(), genotype.tolist(), (scale * covariate).reshape(-1, 1).tolist())
    assert fitted is not None
    assert fitted[0] == pytest.approx(coefficients[-1], rel=1e-7, abs=1e-10)
    assert fitted[4] == pytest.approx(expected_se, rel=1e-7, abs=1e-12)
    assert fitted[3] == pytest.approx(expected_p, rel=1e-7, abs=1e-12)
