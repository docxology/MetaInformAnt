"""SNP heritability estimation and partitioning.

Estimates narrow-sense SNP heritability (h2_SNP) using REML variance component
estimation on genomic relationship matrices. Supports per-chromosome partitioning
and visualization of heritability contributions.

The core model:
    y = mu + g + e
    where Var(g) = sigma_g^2 * K, Var(e) = sigma_e^2 * I
    h2 = sigma_g^2 / (sigma_g^2 + sigma_e^2)

Reference: Yang et al. (2011) Nature Genetics 43:519-525 (GCTA-GREML).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore


def estimate_heritability(
    kinship_matrix: Any,
    phenotypes: List[float],
    method: str = "reml",
) -> Dict[str, Any]:
    """Estimate SNP heritability (h2_SNP) using REML variance component estimation.

    Eigendecomposes the kinship matrix K, rotates phenotypes into the
    eigenspace, then performs a grid search over h2 in [0, 1] to maximize
    the restricted log-likelihood.

    Args:
        kinship_matrix: Kinship/GRM matrix (n x n), numpy array or list of lists.
        phenotypes: Phenotype values for each sample.
        method: Estimation method. Currently only "reml" is supported.

    Returns:
        Dictionary with status, h2 estimate, standard error, variance components,
        log-likelihood, sample size, and method used.
    """
    if not HAS_NUMPY:
        return {"status": "error", "message": "numpy is required for heritability estimation"}

    n = len(phenotypes)
    if n < 3:
        return {
            "status": "error",
            "message": f"Need at least 3 samples for heritability estimation, got {n}",
        }

    try:
        y = np.array(phenotypes, dtype=float)
        K = np.asarray(kinship_matrix, dtype=float)
    except (ValueError, TypeError) as e:
        return {"status": "error", "message": f"Failed to convert inputs to arrays: {e}"}

    if K.shape != (n, n):
        return {
            "status": "error",
            "message": f"Kinship matrix shape {K.shape} does not match {n} samples",
        }

    # Check for zero-variance phenotypes
    if np.var(y) < 1e-12:
        return {
            "status": "error",
            "message": "Phenotype has zero or near-zero variance",
        }

    logger.debug(f"Estimating heritability with {n} samples using method={method}")

    try:
        # Eigendecompose K
        eigenvalues, eigenvectors = np.linalg.eigh(K)
        # Sort descending
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        # Clamp small negative eigenvalues
        eigenvalues = np.maximum(eigenvalues, 0.0)

        # Rotate phenotypes
        Ut = eigenvectors.T
        y_rot = Ut @ y

        # Intercept in rotated space
        ones_rot = Ut @ np.ones(n)

        # Grid search over h2 in [0, 1]
        n_grid = 100
        h2_grid = np.linspace(0.001, 0.999, n_grid)
        best_ll = -np.inf
        best_h2 = 0.5

        for h2_candidate in h2_grid:
            ll = _reml_log_likelihood(y_rot, ones_rot, eigenvalues, h2_candidate, n)
            if ll > best_ll:
                best_ll = ll
                best_h2 = h2_candidate

        # Refine with golden section search around the best grid point
        step = 1.0 / n_grid
        lo = max(0.001, best_h2 - step)
        hi = min(0.999, best_h2 + step)
        best_h2 = _golden_section_max(
            lambda h2: _reml_log_likelihood(y_rot, ones_rot, eigenvalues, h2, n),
            lo,
            hi,
        )
        best_ll = _reml_log_likelihood(y_rot, ones_rot, eigenvalues, best_h2, n)

        # Derive variance components from h2
        total_var = float(np.var(y))
        sigma_g = best_h2 * total_var
        sigma_e = (1.0 - best_h2) * total_var

        # Approximate standard error via second derivative of log-likelihood
        h2_se = _approximate_h2_se(y_rot, ones_rot, eigenvalues, best_h2, n)

        return {
            "status": "success",
            "h2": float(best_h2),
            "h2_se": float(h2_se),
            "sigma_g": float(sigma_g),
            "sigma_e": float(sigma_e),
            "log_likelihood": float(best_ll),
            "n_samples": n,
            "method": method,
        }

    except Exception as e:
        logger.warning(f"Heritability estimation failed: {e}")
        return {"status": "error", "message": str(e)}


def partition_heritability_by_chromosome(
    kinship_matrices: Dict[int, Any],
    phenotypes: List[float],
) -> Dict[str, Any]:
    """Partition SNP heritability by chromosome.

    Estimates per-chromosome h2 by fitting each chromosome's kinship matrix
    independently (one-at-a-time approach). The total h2 is the sum of
    per-chromosome estimates.

    Args:
        kinship_matrices: Maps chromosome number to kinship matrix for
            variants on that chromosome.
        phenotypes: Phenotype values for each sample.

    Returns:
        Dictionary with per-chromosome h2 estimates, total h2, and
        number of chromosomes analyzed.
    """
    if not HAS_NUMPY:
        return {"status": "error", "message": "numpy is required for heritability partitioning"}

    n = len(phenotypes)
    if n < 3:
        return {
            "status": "error",
            "message": f"Need at least 3 samples, got {n}",
        }

    if not kinship_matrices:
        return {"status": "error", "message": "No kinship matrices provided"}

    logger.info(f"Partitioning heritability across {len(kinship_matrices)} chromosomes")

    per_chromosome: Dict[str, Dict[str, float]] = {}
    total_h2 = 0.0

    for chrom, K_chr in sorted(kinship_matrices.items()):
        result = estimate_heritability(K_chr, phenotypes)

        if result["status"] == "success":
            chr_h2 = result["h2"]
            chr_se = result["h2_se"]
            per_chromosome[str(chrom)] = {"h2": chr_h2, "h2_se": chr_se}
            total_h2 += chr_h2
        else:
            logger.warning(f"Chromosome {chrom} estimation failed: {result.get('message', 'unknown')}")
            per_chromosome[str(chrom)] = {"h2": 0.0, "h2_se": 0.0}

    # Cap total h2 at 1.0 (sum of independent estimates can exceed 1)
    total_h2 = min(total_h2, 1.0)

    return {
        "status": "success",
        "per_chromosome": per_chromosome,
        "total_h2": float(total_h2),
        "n_chromosomes": len(kinship_matrices),
    }


def heritability_bar_chart(
    h2_data: Dict[str, Any],
    output_file: Optional[Union[str, Path]] = None,
    title: str = "SNP Heritability by Chromosome",
) -> Dict[str, Any]:
    """Create a bar chart of per-chromosome heritability estimates.

    Plots h2 per chromosome with error bars (h2_se), draws a horizontal
    line at the total h2, and colors bars by relative contribution.

    Args:
        h2_data: Output from partition_heritability_by_chromosome containing
            per_chromosome and total_h2 keys.
        output_file: Path to save the figure. If None, the figure is not saved.
        title: Chart title.

    Returns:
        Dictionary with status and output path (if saved).
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return {"status": "skipped", "output_path": None, "message": "matplotlib not available"}

    per_chromosome = h2_data.get("per_chromosome", {})
    total_h2 = h2_data.get("total_h2", 0.0)

    if not per_chromosome:
        return {"status": "failed", "output_path": None, "message": "No per-chromosome data"}

    # Sort chromosomes numerically where possible
    def _chr_sort_key(c: str) -> Tuple[int, str]:
        try:
            return (0, str(int(c)).zfill(5))
        except ValueError:
            return (1, c)

    sorted_chroms = sorted(per_chromosome.keys(), key=_chr_sort_key)
    h2_values = [per_chromosome[c]["h2"] for c in sorted_chroms]
    h2_se_values = [per_chromosome[c]["h2_se"] for c in sorted_chroms]

    try:
        fig, ax = plt.subplots(figsize=(max(8, len(sorted_chroms) * 0.6), 5))

        # Color bars by relative contribution
        max_h2 = max(h2_values) if max(h2_values) > 0 else 1.0
        colors = plt.cm.YlOrRd([v / max_h2 for v in h2_values])  # type: ignore[attr-defined]

        x_positions = range(len(sorted_chroms))
        bars = ax.bar(
            x_positions,
            h2_values,
            yerr=h2_se_values,
            capsize=3,
            color=colors,
            edgecolor="gray",
            linewidth=0.5,
        )

        # Horizontal line at total h2
        ax.axhline(y=total_h2, color="steelblue", linestyle="--", linewidth=1.5, label=f"Total h2 = {total_h2:.3f}")

        ax.set_xlabel("Chromosome")
        ax.set_ylabel("Heritability (h2)")
        ax.set_title(title)
        ax.set_xticks(list(x_positions))
        ax.set_xticklabels(sorted_chroms, rotation=45 if len(sorted_chroms) > 10 else 0)
        ax.legend(loc="upper right")
        ax.set_ylim(bottom=0)

        plt.tight_layout()

        output_path_str: Optional[str] = None
        if output_file is not None:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=150, bbox_inches="tight")
            output_path_str = str(output_path)
            logger.info(f"Heritability bar chart saved to {output_path}")

        plt.close(fig)

        return {"status": "success", "output_path": output_path_str}

    except Exception as e:
        logger.warning(f"Heritability bar chart failed: {e}")
        return {"status": "failed", "output_path": None, "message": str(e)}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _reml_log_likelihood(
    y_rot: Any,
    ones_rot: Any,
    eigenvalues: Any,
    h2: float,
    n: int,
) -> float:
    """Compute the restricted log-likelihood for a given h2.

    In the rotated space with V = h2 * diag(lambda) + (1-h2) * I:

    Args:
        y_rot: Rotated phenotypes (U' y).
        ones_rot: Rotated intercept vector (U' 1).
        eigenvalues: Eigenvalues of K.
        h2: Candidate heritability value in (0, 1).
        n: Number of samples.

    Returns:
        REML log-likelihood value.
    """
    # Diagonal of V in eigenspace: h2 * lambda_i + (1 - h2)
    d = h2 * eigenvalues + (1.0 - h2)
    d = np.maximum(d, 1e-10)
    d_inv = 1.0 / d

    # Weighted least squares for intercept: beta = (1'V^-1 1)^-1 1'V^-1 y
    ot_dinv_o = float(np.sum(d_inv * ones_rot**2))
    if ot_dinv_o < 1e-30:
        return -1e30
    ot_dinv_y = float(np.sum(d_inv * ones_rot * y_rot))
    beta_hat = ot_dinv_y / ot_dinv_o

    residuals = y_rot - beta_hat * ones_rot
    p = 1  # one fixed effect (intercept)

    # Profiled sigma2
    sigma2 = float(np.sum(d_inv * residuals**2)) / (n - p)
    if sigma2 <= 0:
        return -1e30

    # REML log-likelihood
    ll = -0.5 * (n - p) * math.log(2 * math.pi * sigma2)
    ll -= 0.5 * float(np.sum(np.log(d)))
    ll -= 0.5 * (n - p)
    ll -= 0.5 * math.log(ot_dinv_o)

    return float(ll)


def _approximate_h2_se(
    y_rot: Any,
    ones_rot: Any,
    eigenvalues: Any,
    h2: float,
    n: int,
    delta: float = 0.005,
) -> float:
    """Approximate standard error of h2 using the curvature of the log-likelihood.

    Computes SE = 1 / sqrt(-d2L/dh2^2) using finite differences.

    Args:
        y_rot: Rotated phenotypes.
        ones_rot: Rotated intercept.
        eigenvalues: Eigenvalues of K.
        h2: Estimated h2.
        n: Number of samples.
        delta: Step size for finite difference.

    Returns:
        Approximate standard error of h2.
    """
    h2_lo = max(0.001, h2 - delta)
    h2_hi = min(0.999, h2 + delta)
    actual_delta = (h2_hi - h2_lo) / 2.0

    ll_lo = _reml_log_likelihood(y_rot, ones_rot, eigenvalues, h2_lo, n)
    ll_mid = _reml_log_likelihood(y_rot, ones_rot, eigenvalues, h2, n)
    ll_hi = _reml_log_likelihood(y_rot, ones_rot, eigenvalues, h2_hi, n)

    # Second derivative via central finite difference
    d2ll = (ll_hi - 2.0 * ll_lo + ll_lo) / (actual_delta**2)
    # Correct formula: (ll_hi - 2*ll_mid + ll_lo) / delta^2
    d2ll = (ll_hi - 2.0 * ll_mid + ll_lo) / (actual_delta**2)

    if d2ll >= 0:
        # Likelihood surface is not concave here; return a large SE
        return 0.5

    fisher_info = -d2ll
    se = 1.0 / math.sqrt(fisher_info)

    # Cap SE to reasonable range
    return min(se, 0.5)


def _golden_section_max(
    func: Any,
    a: float,
    b: float,
    tol: float = 1e-6,
    max_iter: int = 100,
) -> float:
    """Golden section search for the maximum of a unimodal function on [a, b].

    Args:
        func: Function to maximize (float -> float).
        a: Lower bound.
        b: Upper bound.
        tol: Convergence tolerance.
        max_iter: Maximum iterations.

    Returns:
        x value at maximum.
    """
    gr = (math.sqrt(5) + 1) / 2

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

    return (a + b) / 2.0
