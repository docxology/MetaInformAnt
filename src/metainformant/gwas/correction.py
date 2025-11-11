"""Multiple testing correction for GWAS."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from ..core.logging import get_logger

logger = get_logger(__name__)

# Try importing scipy for statistical functions
try:
    from scipy.stats import chi2

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    chi2 = None


def bonferroni_correction(
    pvalues: list[float],
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Apply Bonferroni correction for multiple testing.

    Args:
        pvalues: List of p-values
        alpha: Significance level

    Returns:
        Dictionary with:
        - corrected_alpha: Bonferroni-corrected significance threshold
        - significant_count: Number of significant results
        - significant_indices: Indices of significant results
    """
    if not pvalues:
        return {
            "status": "failed",
            "error": "No p-values provided",
        }

    n_tests = len(pvalues)
    corrected_alpha = alpha / n_tests if n_tests > 0 else alpha

    significant_indices = [i for i, p in enumerate(pvalues) if p < corrected_alpha]
    significant_count = len(significant_indices)

    return {
        "status": "success",
        "method": "bonferroni",
        "n_tests": n_tests,
        "alpha": alpha,
        "corrected_alpha": corrected_alpha,
        "significant_count": significant_count,
        "significant_indices": significant_indices,
    }


def fdr_correction(
    pvalues: list[float],
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Apply Benjamini-Hochberg FDR correction.

    Args:
        pvalues: List of p-values
        alpha: FDR level

    Returns:
        Dictionary with:
        - corrected_pvalues: FDR-adjusted p-values
        - significant_count: Number of significant results
        - significant_indices: Indices of significant results
    """
    if not pvalues:
        return {
            "status": "failed",
            "error": "No p-values provided",
        }

    n_tests = len(pvalues)
    pvalues_arr = np.array(pvalues)

    # Sort p-values and get original indices
    sorted_indices = np.argsort(pvalues_arr)
    sorted_pvalues = pvalues_arr[sorted_indices]

    # Compute FDR-adjusted p-values (BH procedure)
    # p_adj[i] = min(p_adj[i+1], p[i] * n / (i+1))
    adjusted_pvalues = np.zeros(n_tests)
    adjusted_pvalues[n_tests - 1] = sorted_pvalues[n_tests - 1]

    for i in range(n_tests - 2, -1, -1):
        adjusted_pvalues[i] = min(
            adjusted_pvalues[i + 1],
            sorted_pvalues[i] * n_tests / (i + 1),
        )

    # Clamp to 1.0
    adjusted_pvalues = np.minimum(adjusted_pvalues, 1.0)

    # Map back to original order
    corrected_pvalues = np.zeros(n_tests)
    for i in range(n_tests):
        corrected_pvalues[sorted_indices[i]] = adjusted_pvalues[i]

    # Find significant results
    significant_indices = [i for i, p in enumerate(corrected_pvalues) if p < alpha]
    significant_count = len(significant_indices)

    return {
        "status": "success",
        "method": "fdr_bh",
        "n_tests": n_tests,
        "alpha": alpha,
        "corrected_pvalues": corrected_pvalues.tolist(),
        "significant_count": significant_count,
        "significant_indices": significant_indices,
    }


def genomic_control(
    chi2_stats: list[float] | None = None,
    pvalues: list[float] | None = None,
) -> dict[str, Any]:
    """Calculate genomic inflation factor (lambda_GC).

    Args:
        chi2_stats: Optional list of chi-square statistics
        pvalues: Optional list of p-values (converted to chi-square if chi2_stats not provided)

    Returns:
        Dictionary with:
        - lambda_gc: Genomic inflation factor
        - median_chi2: Median chi-square statistic
        - median_chi2_expected: Expected median chi-square under null (≈ 0.456)
    """
    if chi2_stats:
        chi2_arr = np.array(chi2_stats)
    elif pvalues:
        # Convert p-values to chi-square statistics
        # chi2 = -2 * ln(p) for df=1 (two-tailed test)
        pvalues_arr = np.array(pvalues)
        chi2_arr = -2.0 * np.log(np.maximum(pvalues_arr, 1e-300))
    else:
        return {
            "status": "failed",
            "error": "Either chi2_stats or pvalues must be provided",
        }

    # Remove invalid values
    valid_mask = np.isfinite(chi2_arr) & (chi2_arr > 0)
    if not np.any(valid_mask):
        return {
            "status": "failed",
            "error": "No valid chi-square statistics",
        }

    chi2_valid = chi2_arr[valid_mask]
    median_chi2 = np.median(chi2_valid)

    # Expected median for chi-square with 1 df is approximately 0.456
    median_chi2_expected = 0.456
    lambda_gc = median_chi2 / median_chi2_expected

    return {
        "status": "success",
        "lambda_gc": float(lambda_gc),
        "median_chi2": float(median_chi2),
        "median_chi2_expected": median_chi2_expected,
        "n_tests": int(np.sum(valid_mask)),
    }


def permutation_test(
    genotypes: list[list[int]],
    phenotypes: list[float],
    covariates: list[list[float]] | None = None,
    n_perm: int = 1000,
    seed: int | None = None,
) -> dict[str, Any]:
    """Perform permutation-based p-value calculation.

    Args:
        genotypes: Genotype matrix (samples x variants)
        phenotypes: Phenotype values
        covariates: Optional covariates
        n_perm: Number of permutations
        seed: Random seed

    Returns:
        Dictionary with permutation-based p-values

    Note:
        This is computationally expensive and may take a long time for large datasets.
    """
    logger.info(f"permutation_test: Running {n_perm} permutations")

    if seed is not None:
        np.random.seed(seed)

    num_variants = len(genotypes[0]) if genotypes else 0
    num_samples = len(genotypes) if genotypes else 0

    if num_variants == 0 or num_samples == 0:
        return {
            "status": "failed",
            "error": "Invalid genotype matrix",
        }

    # Calculate observed statistics (simplified - would need full association test)
    # For demonstration, we'll just use correlation
    observed_stats: list[float] = []

    for var_idx in range(num_variants):
        var_genotypes = [genotypes[sample_idx][var_idx] for sample_idx in range(num_samples)]
        # Simple correlation as statistic
        valid_indices = [i for i in range(num_samples) if var_genotypes[i] != -1]
        if len(valid_indices) < 3:
            observed_stats.append(0.0)
            continue

        genos = np.array([var_genotypes[i] for i in valid_indices])
        phenos = np.array([phenotypes[i] for i in valid_indices])

        if np.std(genos) > 0 and np.std(phenos) > 0:
            corr = np.corrcoef(genos, phenos)[0, 1]
            observed_stats.append(abs(corr))
        else:
            observed_stats.append(0.0)

    # Permute phenotypes and compute null distribution
    perm_pvalues: list[float] = []

    for var_idx in range(num_variants):
        observed_stat = observed_stats[var_idx]
        if observed_stat == 0:
            perm_pvalues.append(1.0)
            continue

        null_stats: list[float] = []
        perm_phenos = np.array(phenotypes.copy())

        for perm in range(n_perm):
            # Permute phenotypes
            np.random.shuffle(perm_phenos)

            var_genotypes = [genotypes[sample_idx][var_idx] for sample_idx in range(num_samples)]
            valid_indices = [i for i in range(num_samples) if var_genotypes[i] != -1]
            if len(valid_indices) < 3:
                continue

            genos = np.array([var_genotypes[i] for i in valid_indices])
            phenos_perm = perm_phenos[valid_indices]

            if np.std(genos) > 0 and np.std(phenos_perm) > 0:
                corr = np.corrcoef(genos, phenos_perm)[0, 1]
                null_stats.append(abs(corr))

        if null_stats:
            # P-value = proportion of null stats >= observed
            p_value = sum(1 for ns in null_stats if ns >= observed_stat) / len(null_stats)
            perm_pvalues.append(p_value)
        else:
            perm_pvalues.append(1.0)

    logger.info(f"permutation_test: Completed permutations for {num_variants} variants")

    return {
        "status": "success",
        "method": "permutation",
        "n_perm": n_perm,
        "perm_pvalues": perm_pvalues,
    }

