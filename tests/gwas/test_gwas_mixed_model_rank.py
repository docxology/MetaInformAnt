"""A monomorphic SNP fails the mixed-model test whether or not the solver notices the singularity.

Zero-mocks: the real EMMA per-SNP test runs on a rotated design built here. A rounding-level perturbation
leaves the SNP column collinear with the intercept but the system only nearly singular, which is what
OpenBLAS solved silently on Linux CI.
"""

from __future__ import annotations

import math

import numpy as np

from metainformant.gwas.analysis.mixed_model import _emma_test_snp


def test_a_numerically_monomorphic_snp_is_reported_singular() -> None:
    rng = np.random.default_rng(11)
    n = 12
    rotation, _ = np.linalg.qr(rng.normal(size=(n, n)))  # orthogonal, like the kinship eigenvectors
    intercept = np.ones(n)
    snp = np.ones(n) + 1e-15 * rng.normal(size=n)  # heterozygous everywhere, up to rounding error
    phenotype = rng.normal(size=n)

    result = _emma_test_snp(
        snp_rot=rotation.T @ snp,
        y_rot=rotation.T @ phenotype,
        X_rot=(rotation.T @ intercept).reshape(-1, 1),
        eigenvalues=np.ones(n),
        delta=1.0,
    )

    assert "singular" in result["error"]
    assert all(math.isnan(result[key]) for key in ("beta", "se", "t_stat", "p_value"))


def test_a_polymorphic_snp_is_still_tested() -> None:
    rng = np.random.default_rng(12)
    n = 12
    rotation, _ = np.linalg.qr(rng.normal(size=(n, n)))
    snp = rng.integers(0, 3, size=n).astype(float)
    phenotype = 0.8 * snp + rng.normal(size=n)

    result = _emma_test_snp(
        snp_rot=rotation.T @ snp,
        y_rot=rotation.T @ phenotype,
        X_rot=(rotation.T @ np.ones(n)).reshape(-1, 1),
        eigenvalues=np.ones(n),
        delta=1.0,
    )

    assert "error" not in result
    assert 0.0 < result["p_value"] < 1.0
