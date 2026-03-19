"""LD Score Regression (LDSR) — SNP-heritability and intercept estimation.

Implements a simplified version of the LDSC method (Bulik-Sullivan et al.
2015 Nat Genet). Requires only summary statistics and per-SNP LD scores
(sum of r² in a cis-window), which are computed from the genotype matrix.

Method:
  E[χ²_j] = (N h²_g / M) × ℓ_j + Na + 1
  where ℓ_j = sum_k r²_jk (LD score for variant j)
        h²_g = SNP heritability
        N    = sample size
        M    = number of variants
        a    = confounding factor (from intercept)

The regression of χ² on LD scores estimates:
  Slope → (N h²_g / M)  →  h²_g = slope × M / N
  Intercept → Na + 1    →  confounding λ

References:
  Bulik-Sullivan et al. (2015) LD Score Regression distinguishes
    confounding from polygenicity. Nat Genet 47:291-295.
  Finucane et al. (2015) Partitioning heritability by functional
    annotation using GWAS summary statistics. Nat Genet 47:1228-1235.
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

try:
    from scipy import stats as scipy_stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    import statsmodels.api as sm
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def compute_ld_scores(
    genotype_matrix: List[List[float]],
    positions: Optional[List[int]] = None,
    window_kb: int = 1000,
    max_variants: int = 2000,
) -> List[float]:
    """Compute per-SNP LD scores: ℓ_j = Σ_k r²_{jk} within window.

    Args:
        genotype_matrix: Variant-major dosage matrix (variants × samples).
        positions: Genomic positions in bp for each variant.
        window_kb: Half-window in kb for LD score computation (default 1000 kb).
        max_variants: Cap on variants to process (for performance).

    Returns:
        List of LD scores, one per variant. Returns [1.0, ...] as fallback.
    """
    if not HAS_NUMPY:
        logger.warning("numpy required for LD score computation")
        return [1.0] * len(genotype_matrix)

    n_var = min(len(genotype_matrix), max_variants)
    if n_var < 2:
        return [1.0] * n_var

    geno = np.array(genotype_matrix[:n_var], dtype=float)
    # Center and scale
    means = geno.mean(axis=1, keepdims=True)
    geno -= means
    stds = geno.std(axis=1, keepdims=True)
    stds[stds < 1e-10] = 1.0
    geno /= stds

    n_samples = geno.shape[1]
    window_bp = window_kb * 1000

    ld_scores = []
    for i in range(n_var):
        total_r2 = 0.0
        for j in range(n_var):
            if i == j:
                total_r2 += 1.0
                continue
            # Check window
            if positions is not None:
                if abs(positions[i] - positions[j]) > window_bp:
                    continue
            # r² = (corr)²
            cov = float(np.dot(geno[i], geno[j])) / n_samples
            r2 = min(1.0, cov ** 2)
            total_r2 += r2
        ld_scores.append(total_r2)

    logger.info(f"LD scores: {n_var} variants, mean l_score={sum(ld_scores)/len(ld_scores):.2f}, "
                f"max={max(ld_scores):.2f}")
    return ld_scores


def ldsr_regression(
    chi2_stats: List[float],
    ld_scores: List[float],
    n_samples: int,
    n_variants: int,
    intercept_constraint: Optional[float] = None,
) -> Dict[str, Any]:
    """LD score regression of χ² statistics on LD scores.

    Fits: χ²_j = slope × ℓ_j + intercept + error
    Recovers:
      h²_snp = slope × M / N   (SNP heritability)
      intercept ≈ 1 + confounding

    Args:
        chi2_stats: Per-variant chi-square statistics (β²/se²).
        ld_scores: Per-variant LD scores from compute_ld_scores().
        n_samples: GWAS sample size N.
        n_variants: Total marker count M.
        intercept_constraint: If set, constrain intercept to this value
            (useful for meta-analysis where confounding is known).

    Returns:
        Dict with:
            h2_snp: float — SNP-based heritability estimate
            h2_snp_se: float — standard error on h2_snp
            intercept: float — estimated intercept (should be ~1.0)
            intercept_se: float
            lambda_gc_explained: float — fraction of λ_GC from polygenicity
            mean_chi2: float — empirical mean χ²
            attenuation_ratio: float — (intercept-1)/(mean_chi2-1)
            method: str
    """
    n = min(len(chi2_stats), len(ld_scores))
    if n < 5:
        return {"h2_snp": 0.0, "intercept": 1.0, "method": "insufficient_data"}

    chi2_arr = [max(c, 0.0) for c in chi2_stats[:n]]
    ld_arr = ld_scores[:n]
    mean_chi2 = sum(chi2_arr) / n

    if HAS_STATSMODELS and HAS_NUMPY:
        X = np.column_stack([ld_arr, [1.0] * n])
        y = np.array(chi2_arr)

        if intercept_constraint is not None:
            # Constrained: regress (chi2 - intercept) on LD scores
            y_adj = y - intercept_constraint
            X_c = np.array(ld_arr).reshape(-1, 1)
            ols = sm.OLS(y_adj, X_c).fit()
            slope = float(ols.params[0])
            slope_se = float(ols.bse[0])
            intercept = intercept_constraint
            intercept_se = 0.0
        else:
            ols = sm.OLS(y, X).fit()
            slope = float(ols.params[0])
            slope_se = float(ols.bse[0])
            intercept = float(ols.params[1])
            intercept_se = float(ols.bse[1])

        method = "statsmodels_OLS"
    elif HAS_SCIPY and HAS_NUMPY:
        slope, intercept, r_val, p_val, se = scipy_stats.linregress(ld_arr, chi2_arr)
        slope_se = se
        intercept_se = se * float(np.std(ld_arr))
        method = "scipy_linregress"
    else:
        # Fallback: manual OLS
        n_f = float(n)
        mean_ld = sum(ld_arr) / n_f
        mean_chi2_v = mean_chi2
        slope_num = sum((ld_arr[i] - mean_ld) * (chi2_arr[i] - mean_chi2_v) for i in range(n))
        slope_den = sum((v - mean_ld) ** 2 for v in ld_arr)
        slope = slope_num / slope_den if slope_den > 1e-10 else 0.0
        intercept = mean_chi2_v - slope * mean_ld
        slope_se = intercept_se = float("nan")
        method = "manual_OLS"

    # SNP heritability
    h2_snp = slope * n_variants / n_samples
    h2_snp_se = slope_se * n_variants / n_samples if not math.isnan(slope_se) else float("nan")
    h2_snp = max(0.0, min(1.0, h2_snp))

    # Attenuation ratio: fraction of χ² inflation from confounding
    # (intercept - 1) / (mean_chi2 - 1)  — should be ~0 for clean GWAS
    if mean_chi2 > 1.001:
        attenuation_ratio = (intercept - 1.0) / (mean_chi2 - 1.0)
    else:
        attenuation_ratio = 0.0

    # Fraction of λ_GC explained by polygenicity (vs. confounding)
    # λ_GC ≈ median chi2 / 0.4549
    median_chi2 = sorted(chi2_arr)[n // 2]
    lambda_gc = median_chi2 / 0.4549
    polygenicity_fraction = 1.0 - max(0.0, (intercept - 1.0)) / max(0.001, lambda_gc - 1.0)

    logger.info(
        f"LDSR: h²_snp={h2_snp:.4f}±{h2_snp_se:.4f}, "
        f"intercept={intercept:.4f}±{intercept_se:.4f}, "
        f"attenuation={attenuation_ratio:.3f} ({method})"
    )

    return {
        "h2_snp": round(h2_snp, 6),
        "h2_snp_se": round(h2_snp_se, 6) if not math.isnan(h2_snp_se) else None,
        "intercept": round(intercept, 4),
        "intercept_se": round(intercept_se, 4) if not math.isnan(intercept_se) else None,
        "slope": round(slope, 6),
        "slope_se": round(slope_se, 6) if not math.isnan(slope_se) else None,
        "mean_chi2": round(mean_chi2, 4),
        "lambda_gc_median": round(lambda_gc, 4),
        "attenuation_ratio": round(attenuation_ratio, 4),
        "polygenicity_fraction": round(polygenicity_fraction, 4),
        "n_snps_used": n,
        "method": method,
    }


def ldsr_plot(
    chi2_stats: List[float],
    ld_scores: List[float],
    ldsr_result: Dict[str, Any],
    output_path: Optional[str] = None,
) -> Any:
    """Scatter plot of χ² vs LD score with fitted LDSR line.

    Args:
        chi2_stats: Per-variant chi-square values.
        ld_scores: Per-variant LD scores.
        ldsr_result: Output of ldsr_regression().
        output_path: File path to save.

    Returns:
        matplotlib Figure or None.
    """
    if not HAS_MATPLOTLIB or not HAS_NUMPY:
        return None

    n = min(len(chi2_stats), len(ld_scores))
    chi2_arr = chi2_stats[:n]
    ld_arr = ld_scores[:n]

    slope = ldsr_result.get("slope", 0.0)
    intercept = ldsr_result.get("intercept", 1.0)
    h2_snp = ldsr_result.get("h2_snp", 0.0)
    attenuation = ldsr_result.get("attenuation_ratio", 0.0)

    fig, ax = plt.subplots(figsize=(8, 6))

    # Scatter (subsample for display)
    if n > 3000:
        import random as _rnd
        idx = _rnd.sample(range(n), 3000)
        xs = [ld_arr[i] for i in idx]
        ys = [chi2_arr[i] for i in idx]
    else:
        xs, ys = ld_arr, chi2_arr

    ax.scatter(xs, ys, alpha=0.15, s=4, color="#4C72B0", label=f"SNPs (n={n:,})")

    # Regression line
    if xs:
        x_lo, x_hi = min(xs), max(xs)
        ax.plot([x_lo, x_hi], [slope * x_lo + intercept, slope * x_hi + intercept],
                "r-", linewidth=2, label=f"LDSR fit (slope={slope:.4f})")

    # Annotations
    ax.axhline(1.0, color="#888", linewidth=1, linestyle=":", label="χ²=1 (null)")
    ax.set_xlabel("LD score (l_j)", fontsize=12)
    ax.set_ylabel("χ² statistic", fontsize=12)
    ax.set_title(
        f"LD Score Regression\nh²_SNP={h2_snp:.4f}, "
        f"intercept={intercept:.3f}, attenuation={attenuation:.3f}",
        fontsize=11,
    )
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Clip y-axis for clarity (top hit may have huge chi2)
    y_clips = sorted(ys)
    if len(y_clips) > 10:
        ax.set_ylim(0, min(y_clips[int(len(y_clips) * 0.99)], 50))

    plt.tight_layout()
    if output_path:
        from pathlib import Path
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"LDSR plot saved to {output_path}")
    return fig
