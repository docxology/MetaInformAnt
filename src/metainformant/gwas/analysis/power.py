"""GWAS statistical power estimation, convergence, and saturation analysis.

Provides functions for:
- NCP-based power calculation (non-central chi-squared)
- Power curves across sample sizes
- SNP/sample subsampling convergence analysis
- Block-jackknife standard error estimation
- Saturation curve fitting (exponential model)

Reference: Sham & Purcell (2014) Nature Reviews Genetics 15:335-346.
"""

from __future__ import annotations

import math
import random
from typing import Any, Dict, List, Optional

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore

# SciPy is needed for non-central chi-squared CDF
try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None  # type: ignore


# ---------------------------------------------------------------------------
# 1. Power calculation
# ---------------------------------------------------------------------------


def compute_power(
    n_samples: int,
    maf: float,
    beta: float,
    alpha: float = 5e-8,
    h2: Optional[float] = None,
) -> Dict[str, Any]:
    """Compute statistical power for a single-SNP GWAS test.

    Uses the non-centrality parameter (NCP) of a chi-squared distribution
    with 1 degree of freedom:

        NCP = N × 2 × MAF × (1 - MAF) × β²
        Power = 1 - F_{χ²(1, NCP)}(χ²_{1,1-α})

    Args:
        n_samples: Number of samples.
        maf: Minor allele frequency (0 < MAF ≤ 0.5).
        beta: Effect size (regression coefficient).
        alpha: Significance threshold (default 5e-8 for genome-wide).
        h2: Optional heritability. If provided, residual variance is adjusted
            to (1 - h2) rather than 1.0.

    Returns:
        Dictionary with power, ncp, threshold, and input parameters.
    """
    if n_samples < 1:
        return {"status": "error", "message": "n_samples must be >= 1"}
    if not (0 < maf <= 0.5):
        return {"status": "error", "message": "MAF must be in (0, 0.5]"}
    if alpha <= 0 or alpha >= 1:
        return {"status": "error", "message": "alpha must be in (0, 1)"}

    # Residual variance (standardized trait)
    sigma_e_sq = (1.0 - h2) if h2 is not None and 0 < h2 < 1 else 1.0

    # Non-centrality parameter
    ncp = n_samples * 2.0 * maf * (1.0 - maf) * (beta**2) / sigma_e_sq

    # Chi-squared threshold for significance
    if HAS_SCIPY:
        threshold = scipy_stats.chi2.ppf(1.0 - alpha, df=1)
        power = 1.0 - scipy_stats.ncx2.cdf(threshold, df=1, nc=ncp)
    else:
        # Fallback: approximate using normal distribution
        z_alpha = _norm_ppf(1.0 - alpha / 2.0)
        z_power = math.sqrt(ncp) - z_alpha
        power = _norm_cdf(z_power)
        threshold = z_alpha**2

    logger.debug(
        f"Power calc: N={n_samples}, MAF={maf:.3f}, β={beta:.3f}, "
        f"α={alpha:.2e} → NCP={ncp:.2f}, power={power:.4f}"
    )

    return {
        "status": "success",
        "power": float(power),
        "ncp": float(ncp),
        "threshold": float(threshold),
        "n_samples": n_samples,
        "maf": float(maf),
        "beta": float(beta),
        "alpha": float(alpha),
        "h2": h2,
    }


def power_curve(
    sample_sizes: List[int],
    maf: float,
    beta: float,
    alpha: float = 5e-8,
    h2: Optional[float] = None,
) -> Dict[str, Any]:
    """Compute power across a range of sample sizes.

    Args:
        sample_sizes: List of sample sizes to evaluate.
        maf: Minor allele frequency.
        beta: Effect size.
        alpha: Significance threshold.
        h2: Optional heritability.

    Returns:
        Dictionary with sample_sizes and corresponding power values.
    """
    powers = []
    for n in sample_sizes:
        result = compute_power(n, maf, beta, alpha, h2)
        powers.append(
            result.get("power", 0.0) if result.get("status") == "success" else 0.0
        )

    logger.info(
        f"Power curve: {len(sample_sizes)} sample sizes, MAF={maf:.3f}, β={beta:.3f}, "
        f"power range [{min(powers):.4f}, {max(powers):.4f}]"
    )

    return {
        "status": "success",
        "sample_sizes": sample_sizes,
        "powers": powers,
        "maf": float(maf),
        "beta": float(beta),
        "alpha": float(alpha),
        "h2": h2,
    }


# ---------------------------------------------------------------------------
# 2. Subsampling convergence
# ---------------------------------------------------------------------------


def subsample_convergence(
    genotype_matrix: List[List[int]],
    traits: List[float],
    variants_info: Optional[List[Dict[str, Any]]] = None,
    fractions: Optional[List[float]] = None,
    n_replicates: int = 5,
    alpha: float = 5e-8,
    config: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Assess GWAS result stability by subsampling variants at different fractions.

    For each fraction f ∈ fractions, randomly selects f×M variants (M total),
    runs association testing, and tracks:
    - Genomic inflation factor (λ_GC)
    - Number of significant hits
    - Mean effect size (|β|)
    - Mean -log10(p)

    Args:
        genotype_matrix: Variant × sample genotype matrix (0/1/2/-1).
        traits: Phenotype values per sample.
        variants_info: Optional variant metadata (chrom, pos, etc.).
        fractions: Data fractions to evaluate (default 10% to 100%).
        n_replicates: Number of random replicates per fraction.
        alpha: Significance threshold.
        config: Optional GWAS config dict.

    Returns:
        Dictionary with fractions and per-fraction mean/std of each metric.
    """
    from metainformant.gwas.analysis.association import association_test_linear
    from metainformant.gwas.analysis.correction import genomic_control

    if fractions is None:
        fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    n_variants = len(genotype_matrix)
    n_samples = len(genotype_matrix[0]) if genotype_matrix else 0

    if n_variants == 0 or n_samples == 0:
        return {"status": "error", "message": "Empty genotype matrix"}

    logger.info(
        f"Convergence analysis: {n_variants} variants, {len(fractions)} fractions, "
        f"{n_replicates} replicates each"
    )

    metrics: Dict[str, List[Dict[str, float]]] = {
        "lambda_gc": [],
        "n_significant": [],
        "mean_abs_beta": [],
        "mean_neg_log10_p": [],
    }

    actual_traits = traits[:n_samples]

    for frac in fractions:
        m_sub = max(1, int(frac * n_variants))
        frac_lambdas = []
        frac_hits = []
        frac_betas = []
        frac_logps = []

        for rep in range(n_replicates):
            # Random subsample of variant indices
            if frac >= 1.0:
                sub_indices = list(range(n_variants))
            else:
                sub_indices = sorted(random.sample(range(n_variants), m_sub))

            # Run association on subsampled variants
            p_values = []
            betas = []
            for idx in sub_indices:
                result = association_test_linear(genotype_matrix[idx], actual_traits)
                pv = result.get("p_value", 1.0)
                b = result.get("beta", 0.0)
                p_values.append(pv)
                betas.append(abs(b))

            # Compute metrics
            n_sig = sum(1 for p in p_values if p < alpha)
            mean_beta = sum(betas) / len(betas) if betas else 0.0
            neg_log10_p = [-math.log10(max(p, 1e-300)) for p in p_values]
            mean_nlp = sum(neg_log10_p) / len(neg_log10_p) if neg_log10_p else 0.0

            # Genomic control lambda
            gc_result = genomic_control(p_values=p_values)
            lambda_val = (
                gc_result.get("lambda_gc", 1.0) if isinstance(gc_result, dict) else 1.0
            )

            frac_lambdas.append(lambda_val)
            frac_hits.append(n_sig)
            frac_betas.append(mean_beta)
            frac_logps.append(mean_nlp)

        # Aggregate over replicates
        metrics["lambda_gc"].append(
            {
                "mean": _mean(frac_lambdas),
                "std": _std(frac_lambdas),
                "values": frac_lambdas,
            }
        )
        metrics["n_significant"].append(
            {
                "mean": _mean(frac_hits),
                "std": _std(frac_hits),
                "values": frac_hits,
            }
        )
        metrics["mean_abs_beta"].append(
            {
                "mean": _mean(frac_betas),
                "std": _std(frac_betas),
                "values": frac_betas,
            }
        )
        metrics["mean_neg_log10_p"].append(
            {
                "mean": _mean(frac_logps),
                "std": _std(frac_logps),
                "values": frac_logps,
            }
        )

        logger.info(
            f"  Fraction {frac:.0%}: λ_GC={_mean(frac_lambdas):.4f}±{_std(frac_lambdas):.4f}, "
            f"hits={_mean(frac_hits):.1f}±{_std(frac_hits):.1f}"
        )

    return {
        "status": "success",
        "fractions": fractions,
        "n_replicates": n_replicates,
        "n_variants_total": n_variants,
        "metrics": metrics,
    }


# ---------------------------------------------------------------------------
# 3. Jackknife SE estimation
# ---------------------------------------------------------------------------


def jackknife_se(
    genotype_matrix: List[List[int]],
    traits: List[float],
    n_blocks: int = 10,
    statistic: str = "lambda_gc",
    alpha: float = 5e-8,
) -> Dict[str, Any]:
    """Estimate standard error of a GWAS summary statistic via block-jackknife.

    Divides variants into n_blocks contiguous blocks, leaving each out in turn,
    recomputing the statistic, then applying the jackknife SE formula:

        SE = sqrt((J-1)/J × Σ(θ_{-j} - θ̄)²)

    Args:
        genotype_matrix: Variant × sample genotype matrix.
        traits: Phenotype values.
        n_blocks: Number of jackknife blocks.
        statistic: Which statistic to estimate SE for.
            Options: "lambda_gc", "n_significant", "mean_abs_beta".
        alpha: Significance threshold.

    Returns:
        Dictionary with full estimate, SE, 95% CI, and block estimates.
    """

    n_variants = len(genotype_matrix)
    n_samples = len(genotype_matrix[0]) if genotype_matrix else 0

    if n_variants < n_blocks:
        return {
            "status": "error",
            "message": f"Need at least {n_blocks} variants for jackknife",
        }

    actual_traits = traits[:n_samples]

    # Divide into blocks
    block_size = n_variants // n_blocks
    blocks = []
    for b in range(n_blocks):
        start = b * block_size
        end = start + block_size if b < n_blocks - 1 else n_variants
        blocks.append(list(range(start, end)))

    logger.info(f"Jackknife SE: {n_blocks} blocks, statistic={statistic}")

    # Full-dataset estimate
    full_stat = _compute_statistic(
        genotype_matrix, actual_traits, list(range(n_variants)), statistic, alpha
    )

    # Leave-one-block-out estimates
    block_estimates = []
    for b in range(n_blocks):
        # All indices except block b
        keep_indices = []
        for bb in range(n_blocks):
            if bb != b:
                keep_indices.extend(blocks[bb])

        theta_minus_j = _compute_statistic(
            genotype_matrix, actual_traits, keep_indices, statistic, alpha
        )
        block_estimates.append(theta_minus_j)

    # Jackknife SE
    theta_bar = sum(block_estimates) / n_blocks
    variance = (
        (n_blocks - 1.0) / n_blocks * sum((t - theta_bar) ** 2 for t in block_estimates)
    )
    se = math.sqrt(variance) if variance > 0 else 0.0

    # 95% CI
    z = 1.96
    ci_lower = full_stat - z * se
    ci_upper = full_stat + z * se

    logger.info(
        f"Jackknife {statistic}: estimate={full_stat:.4f}, SE={se:.4f}, "
        f"95% CI=[{ci_lower:.4f}, {ci_upper:.4f}]"
    )

    return {
        "status": "success",
        "statistic": statistic,
        "estimate": float(full_stat),
        "se": float(se),
        "ci_lower": float(ci_lower),
        "ci_upper": float(ci_upper),
        "n_blocks": n_blocks,
        "block_estimates": [float(t) for t in block_estimates],
    }


# ---------------------------------------------------------------------------
# 4. Saturation analysis
# ---------------------------------------------------------------------------


def saturation_analysis(
    convergence_data: Dict[str, Any],
    metric: str = "n_significant",
    convergence_threshold: float = 0.05,
) -> Dict[str, Any]:
    """Fit saturation curve to convergence data.

    Fits the exponential model K(f) = K_∞ × (1 - e^{-κf}) to the observed
    metric values and determines whether results have saturated.

    Saturation is declared when the marginal gain from the last 10% of data
    is below the convergence_threshold fraction of K_∞.

    Args:
        convergence_data: Output from subsample_convergence().
        metric: Which metric to analyze (key in convergence_data["metrics"]).
        convergence_threshold: Fraction of K_∞ below which marginal gain
            indicates saturation.

    Returns:
        Dictionary with fitted parameters, R², saturation flag.
    """
    fractions = convergence_data.get("fractions", [])
    metric_data = convergence_data.get("metrics", {}).get(metric, [])

    if not fractions or not metric_data:
        return {"status": "error", "message": f"No data for metric '{metric}'"}

    observed = [m["mean"] for m in metric_data]

    if all(v == 0 for v in observed):
        return {
            "status": "success",
            "k_inf": 0.0,
            "kappa": 0.0,
            "r_squared": 1.0,
            "is_saturated": True,
            "marginal_gain_at_100pct": 0.0,
            "metric": metric,
            "message": "All values are zero; trivially saturated",
        }

    # Fit K_∞ and κ via grid search (no scipy.optimize dependency)
    k_inf_est = max(observed) * 1.2 if max(observed) > 0 else 1.0
    best_r2 = -1e30
    best_k_inf = k_inf_est
    best_kappa = 1.0

    for k_inf_mult in [0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0]:
        k_inf_try = max(observed) * k_inf_mult if max(observed) > 0 else 1.0
        for kappa_try in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0]:
            predicted = [
                k_inf_try * (1.0 - math.exp(-kappa_try * f)) for f in fractions
            ]
            r2 = _r_squared(observed, predicted)
            if r2 > best_r2:
                best_r2 = r2
                best_k_inf = k_inf_try
                best_kappa = kappa_try

    # Marginal gain from 0.9 → 1.0
    val_90 = best_k_inf * (1.0 - math.exp(-best_kappa * 0.9))
    val_100 = best_k_inf * (1.0 - math.exp(-best_kappa * 1.0))
    marginal_gain = (val_100 - val_90) / best_k_inf if best_k_inf > 0 else 0.0

    is_saturated = marginal_gain < convergence_threshold

    logger.info(
        f"Saturation analysis ({metric}): K_∞={best_k_inf:.4f}, κ={best_kappa:.2f}, "
        f"R²={best_r2:.4f}, saturated={is_saturated}"
    )

    return {
        "status": "success",
        "k_inf": float(best_k_inf),
        "kappa": float(best_kappa),
        "r_squared": float(best_r2),
        "is_saturated": is_saturated,
        "marginal_gain_at_100pct": float(marginal_gain),
        "metric": metric,
        "fractions": fractions,
        "observed": observed,
        "predicted": [
            best_k_inf * (1.0 - math.exp(-best_kappa * f)) for f in fractions
        ],
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _compute_statistic(
    genotype_matrix: List[List[int]],
    traits: List[float],
    variant_indices: List[int],
    statistic: str,
    alpha: float,
) -> float:
    """Compute a single summary statistic on a subset of variants."""
    from metainformant.gwas.analysis.association import association_test_linear
    from metainformant.gwas.analysis.correction import genomic_control

    p_values = []
    betas = []
    for idx in variant_indices:
        result = association_test_linear(genotype_matrix[idx], traits)
        p_values.append(result.get("p_value", 1.0))
        betas.append(abs(result.get("beta", 0.0)))

    if statistic == "lambda_gc":
        gc_result = genomic_control(p_values=p_values)
        return gc_result.get("lambda_gc", 1.0) if isinstance(gc_result, dict) else 1.0
    elif statistic == "n_significant":
        return float(sum(1 for p in p_values if p < alpha))
    elif statistic == "mean_abs_beta":
        return sum(betas) / len(betas) if betas else 0.0
    else:
        return 0.0


def _mean(values: List[float]) -> float:
    """Compute mean of a list of values."""
    if not values:
        return 0.0
    return sum(values) / len(values)


def _std(values: List[float]) -> float:
    """Compute sample standard deviation."""
    if len(values) < 2:
        return 0.0
    m = _mean(values)
    variance = sum((x - m) ** 2 for x in values) / (len(values) - 1)
    return math.sqrt(variance)


def _r_squared(observed: List[float], predicted: List[float]) -> float:
    """Compute R² (coefficient of determination)."""
    if len(observed) != len(predicted) or not observed:
        return 0.0
    mean_obs = sum(observed) / len(observed)
    ss_res = sum((o - p) ** 2 for o, p in zip(observed, predicted))
    ss_tot = sum((o - mean_obs) ** 2 for o in observed)
    if ss_tot == 0:
        return 1.0 if ss_res == 0 else 0.0
    return 1.0 - ss_res / ss_tot


def _norm_cdf(x: float) -> float:
    """Standard normal CDF (no scipy fallback)."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _norm_ppf(p: float) -> float:
    """Approximate standard normal quantile function (Abramowitz & Stegun)."""
    if p <= 0:
        return -10.0
    if p >= 1:
        return 10.0
    if p == 0.5:
        return 0.0

    # Rational approximation
    if p < 0.5:
        t = math.sqrt(-2.0 * math.log(p))
    else:
        t = math.sqrt(-2.0 * math.log(1.0 - p))

    c0 = 2.515517
    c1 = 0.802853
    c2 = 0.010328
    d1 = 1.432788
    d2 = 0.189269
    d3 = 0.001308

    result = t - (c0 + c1 * t + c2 * t * t) / (
        1.0 + d1 * t + d2 * t * t + d3 * t * t * t
    )

    return result if p > 0.5 else -result
