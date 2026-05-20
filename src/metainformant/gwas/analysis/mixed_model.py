"""Mixed linear model (MLM) for GWAS using the EMMA algorithm.

Implements the Efficient Mixed-Model Association (EMMA) algorithm for
genome-wide association testing with population structure correction.

Reference: Kang et al. (2008) Genetics 178:1709-1723.

The core idea:
    y = Xb + Zu + e
    where u ~ N(0, sigma_g^2 * K), e ~ N(0, sigma_e^2 * I)

EMMA eigendecomposes K once, rotates the model, then tests each SNP in O(n) time.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore


def association_test_mixed(
    genotypes: List[int],
    phenotypes: List[float],
    kinship_matrix: List[List[float]],
    covariates: Optional[List[List[float]]] = None,
) -> Dict[str, Any]:
    """Perform mixed model association test for a single SNP.

    Args:
        genotypes: Genotype values (0, 1, 2) for each sample
        phenotypes: Phenotype values for each sample
        kinship_matrix: Kinship matrix (n_samples x n_samples)
        covariates: Optional covariate matrix (covariates x samples)

    Returns:
        Dictionary with association test results
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required for mixed model association testing")

    n = len(genotypes)
    if n != len(phenotypes):
        raise ValueError("Genotypes and phenotypes must have same length")
    if len(kinship_matrix) != n or any(len(row) != n for row in kinship_matrix):
        raise ValueError("Kinship matrix dimensions must match number of samples")

    logger.debug(f"Running mixed model association test on {n} samples")

    try:
        y = np.array(phenotypes, dtype=float)
        K = np.array(kinship_matrix, dtype=float)
        snp = np.array(genotypes, dtype=float)

        # Build covariate matrix: intercept + covariates
        X_cols = [np.ones(n)]
        if covariates:
            for cov in covariates:
                X_cols.append(np.array(cov, dtype=float))
        X = np.column_stack(X_cols)

        # Eigendecompose K
        eigenvalues, eigenvectors = _emma_eigendecompose(K)

        # Rotate y, X
        Ut = eigenvectors.T
        y_rot = Ut @ y
        X_rot = Ut @ X

        # Estimate variance components via REML
        sigma_g, sigma_e = _emma_reml(y_rot, X_rot, eigenvalues)

        # Compute delta = sigma_e / sigma_g (or handle edge cases)
        if sigma_g > 1e-10:
            delta = sigma_e / sigma_g
        else:
            delta = 1e10  # essentially no genetic variance

        heritability = sigma_g / (sigma_g + sigma_e) if (sigma_g + sigma_e) > 1e-10 else 0.0

        # Test this SNP
        snp_rot = Ut @ snp
        result = _emma_test_snp(snp_rot, y_rot, X_rot, eigenvalues, delta)

        return {
            "status": "success",
            "beta": result["beta"],
            "se": result["se"],
            "t_stat": result["t_stat"],
            "p_value": result["p_value"],
            "sigma_g": float(sigma_g),
            "sigma_e": float(sigma_e),
            "heritability": float(heritability),
            "n_samples": n,
            "test_type": "mixed",
            "converged": True,
        }

    except Exception as e:
        logger.warning(f"Mixed model association test failed: {e}")
        return {
            "status": "error",
            "beta": 0.0,
            "se": 0.0,
            "t_stat": 0.0,
            "p_value": 1.0,
            "sigma_g": 0.0,
            "sigma_e": 0.0,
            "heritability": 0.0,
            "n_samples": n,
            "test_type": "mixed",
            "converged": False,
            "error": str(e),
        }


def run_mixed_model_gwas(
    genotype_matrix: List[List[int]],
    phenotypes: List[float],
    kinship_matrix: List[List[float]],
    variant_info: Optional[List[Dict[str, Any]]] = None,
    covariates: Optional[List[List[float]]] = None,
) -> List[Dict[str, Any]]:
    """Run mixed model GWAS across all variants.

    Eigendecomposes K once, then tests each SNP efficiently.

    Args:
        genotype_matrix: Genotype matrix (variants x samples)
        phenotypes: Phenotype values
        kinship_matrix: Kinship matrix (n_samples x n_samples)
        variant_info: Optional variant metadata (chrom, pos, id, etc.)
        covariates: Optional covariate matrix (covariates x samples)

    Returns:
        List of result dictionaries, one per variant
    """
    if not HAS_NUMPY:
        raise ImportError("numpy is required for mixed model GWAS")

    n_variants = len(genotype_matrix)
    n_samples = len(phenotypes)

    logger.info(f"Running mixed model GWAS: {n_variants} variants x {n_samples} samples")

    y = np.array(phenotypes, dtype=float)
    K = np.array(kinship_matrix, dtype=float)

    # Build covariate matrix
    X_cols = [np.ones(n_samples)]
    if covariates:
        for cov in covariates:
            X_cols.append(np.array(cov, dtype=float))
    X = np.column_stack(X_cols)

    # Eigendecompose K once
    eigenvalues, eigenvectors = _emma_eigendecompose(K)
    Ut = eigenvectors.T

    # Rotate y and X once
    y_rot = Ut @ y
    X_rot = Ut @ X

    # Estimate variance components via REML (null model)
    sigma_g, sigma_e = _emma_reml(y_rot, X_rot, eigenvalues)
    delta = sigma_e / sigma_g if sigma_g > 1e-10 else 1e10
    heritability = sigma_g / (sigma_g + sigma_e) if (sigma_g + sigma_e) > 1e-10 else 0.0

    logger.info(f"REML estimates: sigma_g={sigma_g:.4f}, sigma_e={sigma_e:.4f}, h2={heritability:.4f}")

    # Test each SNP
    results = []
    for i in range(n_variants):
        snp = np.array(genotype_matrix[i], dtype=float)
        snp_rot = Ut @ snp
        snp_result = _emma_test_snp(snp_rot, y_rot, X_rot, eigenvalues, delta)

        result: Dict[str, Any] = {
            "variant_index": i,
            "beta": snp_result["beta"],
            "se": snp_result["se"],
            "t_stat": snp_result["t_stat"],
            "p_value": snp_result["p_value"],
            "test_type": "mixed",
        }

        if variant_info and i < len(variant_info):
            result["variant_id"] = variant_info[i].get("id", f"variant_{i}")
            result["chrom"] = variant_info[i].get("chrom", "")
            result["pos"] = variant_info[i].get("pos", 0)

        results.append(result)

    logger.info(f"Mixed model GWAS complete: {len(results)} variants tested")
    return results


def _emma_eigendecompose(K: Any) -> Tuple[Any, Any]:
    """Eigendecompose the kinship matrix K.

    Args:
        K: Kinship matrix (n x n numpy array)

    Returns:
        Tuple of (eigenvalues, eigenvectors) sorted descending
    """
    # Use symmetric eigendecomposition for better numerical stability
    eigenvalues, eigenvectors = np.linalg.eigh(K)

    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Clamp small negative eigenvalues to zero
    eigenvalues = np.maximum(eigenvalues, 0.0)

    return eigenvalues, eigenvectors


def _emma_reml(
    y_rot: Any,
    X_rot: Any,
    eigenvalues: Any,
    n_grid: int = 100,
) -> Tuple[float, float]:
    """Estimate variance components using REML via grid search + optimization.

    In the rotated model:
        y* = X*b + e*
        Var(e*) = sigma_g * diag(lambda + delta)
        delta = sigma_e / sigma_g

    We profile out b and sigma_g, optimize over delta.

    Args:
        y_rot: Rotated phenotypes (U' y)
        X_rot: Rotated covariates (U' X)
        eigenvalues: Eigenvalues of K
        n_grid: Number of grid points for initial search

    Returns:
        Tuple of (sigma_g, sigma_e)
    """
    n = len(y_rot)

    def reml_log_likelihood(log_delta: float) -> float:
        """Compute REML log-likelihood for a given log(delta)."""
        delta = math.exp(log_delta)
        d = eigenvalues + delta  # lambda_i + delta

        # Avoid division by zero
        d = np.maximum(d, 1e-10)

        # Weighted least squares: solve (X' D^-1 X) b = X' D^-1 y
        d_inv = 1.0 / d
        XtDinvX = X_rot.T @ np.diag(d_inv) @ X_rot
        XtDinvy = X_rot.T @ (d_inv * y_rot)

        try:
            beta = np.linalg.solve(XtDinvX, XtDinvy)
        except np.linalg.LinAlgError:
            return -1e30

        residuals = y_rot - X_rot @ beta
        p = X_rot.shape[1]

        # sigma_g (profiled)
        sigma_g = float(np.sum(d_inv * residuals**2) / (n - p))
        if sigma_g <= 0:
            return -1e30

        # REML log-likelihood
        ll = -0.5 * (n - p) * math.log(2 * math.pi * sigma_g)
        ll -= 0.5 * np.sum(np.log(d))
        ll -= 0.5 * (n - p)

        # REML correction for fixed effects
        sign, logdet = np.linalg.slogdet(XtDinvX)
        if sign > 0:
            ll -= 0.5 * logdet

        return float(ll)

    # Grid search over log(delta)
    log_delta_grid = np.linspace(-10, 10, n_grid)
    best_ll = -1e30
    best_log_delta = 0.0

    for ld in log_delta_grid:
        ll = reml_log_likelihood(ld)
        if ll > best_ll:
            best_ll = ll
            best_log_delta = ld

    # Refine with Brent-style golden section search
    best_log_delta = _golden_section_search(reml_log_likelihood, best_log_delta - 1.0, best_log_delta + 1.0)

    # Extract final estimates
    delta = math.exp(best_log_delta)
    d = eigenvalues + delta
    d = np.maximum(d, 1e-10)
    d_inv = 1.0 / d
    p = X_rot.shape[1]

    XtDinvX = X_rot.T @ np.diag(d_inv) @ X_rot
    XtDinvy = X_rot.T @ (d_inv * y_rot)

    try:
        beta = np.linalg.solve(XtDinvX, XtDinvy)
    except np.linalg.LinAlgError:
        return 0.0, 1.0

    residuals = y_rot - X_rot @ beta
    sigma_g = float(np.sum(d_inv * residuals**2) / (n - p))
    sigma_e = sigma_g * delta

    # Ensure non-negative
    sigma_g = max(sigma_g, 0.0)
    sigma_e = max(sigma_e, 0.0)

    return sigma_g, sigma_e


def _golden_section_search(
    func: Any,
    a: float,
    b: float,
    tol: float = 1e-5,
    max_iter: int = 100,
) -> float:
    """Golden section search for maximum of a unimodal function.

    Args:
        func: Function to maximize
        a: Lower bound
        b: Upper bound
        tol: Convergence tolerance
        max_iter: Maximum iterations

    Returns:
        x value at maximum
    """
    gr = (math.sqrt(5) + 1) / 2  # golden ratio

    c = b - (b - a) / gr
    d = a + (b - a) / gr

    for _ in range(max_iter):
        if abs(b - a) < tol:
            break
        if func(c) < func(d):
            a = c
        else:
            b = d
        c = b - (b - a) / gr
        d = a + (b - a) / gr

    return (a + b) / 2


def _emma_test_snp(
    snp_rot: Any,
    y_rot: Any,
    X_rot: Any,
    eigenvalues: Any,
    delta: float,
) -> Dict[str, float]:
    """Test a single SNP using pre-computed EMMA parameters.

    Args:
        snp_rot: Rotated SNP genotypes (U' x)
        y_rot: Rotated phenotypes (U' y)
        X_rot: Rotated covariates (U' X)
        eigenvalues: Eigenvalues of K
        delta: sigma_e / sigma_g ratio

    Returns:
        Dictionary with beta, se, t_stat, p_value
    """
    n = len(y_rot)
    d = eigenvalues + delta
    d = np.maximum(d, 1e-10)
    d_inv = 1.0 / d

    # Add SNP to design matrix
    X_full = np.column_stack([X_rot, snp_rot.reshape(-1, 1)])
    p = X_full.shape[1]

    # Weighted least squares
    XtDinvX = X_full.T @ np.diag(d_inv) @ X_full
    XtDinvy = X_full.T @ (d_inv * y_rot)

    try:
        beta_vec = np.linalg.solve(XtDinvX, XtDinvy)
    except np.linalg.LinAlgError:
        return {"beta": 0.0, "se": 0.0, "t_stat": 0.0, "p_value": 1.0}

    residuals = y_rot - X_full @ beta_vec
    dof = n - p

    if dof <= 0:
        return {"beta": 0.0, "se": 0.0, "t_stat": 0.0, "p_value": 1.0}

    sigma2 = float(np.sum(d_inv * residuals**2) / dof)

    # Variance of beta
    try:
        cov_beta = sigma2 * np.linalg.inv(XtDinvX)
    except np.linalg.LinAlgError:
        return {"beta": 0.0, "se": 0.0, "t_stat": 0.0, "p_value": 1.0}

    # SNP coefficient is the last one
    snp_idx = p - 1
    beta = float(beta_vec[snp_idx])
    se = float(math.sqrt(max(cov_beta[snp_idx, snp_idx], 0.0)))
    t_stat = beta / se if se > 1e-10 else 0.0

    # p-value from t-distribution
    try:
        from scipy import stats as scipy_stats

        p_value = float(2 * scipy_stats.t.sf(abs(t_stat), dof))
    except ImportError:
        # Normal approximation fallback
        p_value = 2 * _normal_cdf_complement(abs(t_stat))

    return {"beta": beta, "se": se, "t_stat": t_stat, "p_value": p_value}


def _normal_cdf_complement(x: float) -> float:
    """Compute 1 - Phi(x) using Abramowitz & Stegun approximation."""
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    x_scaled = abs(x) / math.sqrt(2.0)
    t = 1.0 / (1.0 + p * x_scaled)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-x_scaled * x_scaled)
    cdf = 0.5 * (1 + y) if x >= 0 else 0.5 * (1 - y)
    return 1.0 - cdf
