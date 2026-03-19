"""Polygenic Risk Score (PRS) construction and validation.

Implements the Clumping + Thresholding (C+T) method, the most widely
validated non-Bayesian PRS approach:

1. Clumping: Keep the most significant variant in each LD block (r²<threshold)
2. Thresholding: Apply p-value cutoffs [5e-8, 1e-4, 0.01, 0.05, 0.5]
3. Score: PRS_i = Σ (β_j × dosage_{i,j}) for selected variants
4. Validate: Compute R² (linear) or AUC (binary) on phenotypes

References:
  - Purcell et al. (2009) Nat Genet — original C+T PRS
  - Choi et al. (2020) Nat Protoc — PRS tutorial and best practices
  - Dudbridge (2013) PloS Genet — statistical theory of C+T
"""
from __future__ import annotations

import math
import statistics as _stats
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def clump_variants(
    association_results: List[Dict[str, Any]],
    ld_matrix: List[List[float]],
    r2_threshold: float = 0.1,
    p_threshold: float = 1.0,
) -> List[int]:
    """Greedy LD clumping: keep most significant variant in each LD block.

    Algorithm (mirrors PLINK --clump):
    1. Sort variants by p-value (ascending)
    2. For each variant i in order:
       a. If already pruned, skip
       b. Else add i to clump representatives
       c. Prune all j with r²(i,j) > r2_threshold

    Args:
        association_results: List of result dicts (must have 'p_value' key).
        ld_matrix: Symmetric r² matrix (n_variants × n_variants), 0-indexed
            matching the order of association_results.
        r2_threshold: LD r² above which a variant is pruned (default 0.1).
        p_threshold: Pre-filter: only consider variants with p < this.

    Returns:
        List of integer indices (into association_results) of clump representatives.
    """
    n = len(association_results)
    if n == 0:
        return []

    # Pre-filter by p-value threshold
    eligible = [
        i for i, r in enumerate(association_results)
        if r.get("p_value", 1.0) < p_threshold
    ]
    if not eligible:
        return []

    # Sort by p-value
    eligible_sorted = sorted(eligible, key=lambda i: association_results[i]["p_value"])

    pruned = set()
    representatives = []
    n_ld = len(ld_matrix)

    for idx in eligible_sorted:
        if idx in pruned:
            continue
        representatives.append(idx)
        # Prune all variants in LD with this one
        if idx < n_ld:
            for j in eligible_sorted:
                if j != idx and j not in pruned and j < n_ld and idx < len(ld_matrix[j]):
                    if ld_matrix[idx][j] > r2_threshold:
                        pruned.add(j)

    logger.info(f"Clumping (r²<{r2_threshold}): {len(representatives)}/{len(eligible)} "
                f"representative variants selected from {n} total")
    return representatives


def compute_prs(
    genotype_matrix: List[List[float]],
    betas: List[float],
    selected_indices: List[int],
    n_samples: int,
) -> List[float]:
    """Compute polygenic risk scores for all samples.

    PRS_i = Σ_{j in selected} (β_j × dosage_{i,j}) / n_variants_used

    Args:
        genotype_matrix: Variant-major (variants × samples).
        betas: Effect sizes for each variant, matching genotype_matrix rows.
        selected_indices: Which variants (rows) to include.
        n_samples: Number of samples.

    Returns:
        List of PRS values, one per sample. Normalized by number of variants.
    """
    if not selected_indices:
        return [0.0] * n_samples

    prs = [0.0] * n_samples
    n_used = 0
    for var_idx in selected_indices:
        if var_idx >= len(genotype_matrix) or var_idx >= len(betas):
            continue
        dosages = genotype_matrix[var_idx]
        beta = betas[var_idx]
        for s in range(min(n_samples, len(dosages))):
            d = dosages[s]
            if not math.isnan(d):
                prs[s] += beta * d
        n_used += 1

    # Normalize by variant count (mean score, not cumulative)
    if n_used > 0:
        prs = [v / n_used for v in prs]

    logger.info(f"PRS computed: {n_used} variants, {n_samples} samples")
    return prs


def prs_r2(prs: List[float], phenotypes: List[float]) -> float:
    """Compute Pearson R² between PRS and continuous phenotype.

    Args:
        prs: Polygenic risk scores (one per sample).
        phenotypes: Phenotype values (matched order).

    Returns:
        R² value ∈ [0, 1]. Returns 0.0 if inputs are degenerate.
    """
    n = min(len(prs), len(phenotypes))
    if n < 3:
        return 0.0
    prs = prs[:n]
    y = phenotypes[:n]
    mean_prs = sum(prs) / n
    mean_y = sum(y) / n
    num = sum((prs[i] - mean_prs) * (y[i] - mean_y) for i in range(n))
    den_prs = math.sqrt(sum((v - mean_prs) ** 2 for v in prs))
    den_y = math.sqrt(sum((v - mean_y) ** 2 for v in y))
    if den_prs < 1e-10 or den_y < 1e-10:
        return 0.0
    r = num / (den_prs * den_y)
    return min(1.0, r ** 2)


def prs_full_analysis(
    association_results: List[Dict[str, Any]],
    genotype_matrix: List[List[float]],
    phenotypes: List[float],
    ld_matrix: List[List[float]],
    p_thresholds: Optional[List[float]] = None,
    r2_clump: float = 0.1,
) -> Dict[str, Any]:
    """Full C+T PRS analysis across multiple p-value thresholds.

    Args:
        association_results: GWAS results (must have 'p_value', 'beta').
        genotype_matrix: Variant-major dosage matrix.
        phenotypes: Phenotype values per sample.
        ld_matrix: Pairwise r² matrix.
        p_thresholds: P-value thresholds to test. Default: [5e-8, 1e-4, 0.01, 0.05, 0.5]
        r2_clump: LD r² threshold for clumping.

    Returns:
        Dict with:
            'thresholds': list of dicts {p_threshold, n_variants, r2_phenotype, prs_scores}
            'best_threshold': p-value threshold with highest R²
            'best_r2': highest R²
    """
    if p_thresholds is None:
        p_thresholds = [5e-8, 1e-4, 1e-2, 5e-2, 0.5]

    betas = [r.get("beta", 0.0) for r in association_results]
    n_samples = len(phenotypes)
    threshold_results = []

    for p_thresh in p_thresholds:
        clumped_idx = clump_variants(
            association_results, ld_matrix,
            r2_threshold=r2_clump, p_threshold=p_thresh,
        )
        prs_scores = compute_prs(genotype_matrix, betas, clumped_idx, n_samples)
        r2 = prs_r2(prs_scores, phenotypes)
        threshold_results.append({
            "p_threshold": p_thresh,
            "n_variants": len(clumped_idx),
            "r2_phenotype": round(r2, 6),
            "prs_scores": prs_scores,
        })
        logger.info(f"  PRS (p<{p_thresh:.0e}, n={len(clumped_idx)} variants): R²={r2:.4f}")

    # Best threshold
    best = max(threshold_results, key=lambda x: x["r2_phenotype"])
    logger.info(f"Best PRS: p<{best['p_threshold']:.0e} → R²={best['r2_phenotype']:.4f} "
                f"({best['n_variants']} variants)")

    return {
        "thresholds": threshold_results,
        "best_threshold": best["p_threshold"],
        "best_r2": best["r2_phenotype"],
        "best_n_variants": best["n_variants"],
    }


def prs_distribution_plot(
    prs_result: Dict[str, Any],
    phenotypes: List[float],
    output_path: Optional[Union[str, Path]] = None,
) -> Any:
    """Two-panel PRS visualization.

    Left:  R² vs p-value threshold (N selected variants as secondary bar)
    Right: PRS distribution histogram colored by phenotype tertile (low/mid/high)

    Args:
        prs_result: Output of prs_full_analysis().
        phenotypes: Continuous phenotype values.
        output_path: Where to save the figure.

    Returns:
        matplotlib Figure or None.
    """
    if not HAS_MATPLOTLIB:
        logger.warning("matplotlib not available for PRS plot")
        return None

    thresholds = prs_result.get("thresholds", [])
    best_thresh = prs_result.get("best_threshold")
    best_prs = next((t["prs_scores"] for t in thresholds
                     if t["p_threshold"] == best_thresh), [])

    if not thresholds:
        return None

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f"Polygenic Risk Score — C+T Analysis\nBest R²={prs_result.get('best_r2', 0):.4f} "
        f"(p<{best_thresh:.0e}, n={prs_result.get('best_n_variants', 0)} variants)",
        fontsize=12, fontweight="bold",
    )

    # Left: R² profile across p thresholds
    p_labels = [f"{t['p_threshold']:.0e}" for t in thresholds]
    r2_vals = [t["r2_phenotype"] for t in thresholds]
    n_var_vals = [t["n_variants"] for t in thresholds]
    colors = ["#C44E52" if t["p_threshold"] == best_thresh else "#4C72B0" for t in thresholds]

    bars = ax1.bar(p_labels, r2_vals, color=colors, edgecolor="white", linewidth=0.5, alpha=0.9)
    ax1_twin = ax1.twinx()
    ax1_twin.plot(p_labels, n_var_vals, "ko--", markersize=5, label="N variants", zorder=5)
    ax1_twin.set_ylabel("# variants selected", fontsize=9, color="#333")
    ax1_twin.tick_params(axis="y", labelcolor="#333", labelsize=8)
    ax1.set_xlabel("P-value threshold", fontsize=11)
    ax1.set_ylabel("PRS R² with phenotype", fontsize=11)
    ax1.set_title("PRS Performance vs Threshold", fontsize=10)

    for bar, r2 in zip(bars, r2_vals):
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.0005,
                 f"{r2:.4f}", ha="center", va="bottom", fontsize=7)

    # Right: PRS distribution by phenotype tertile (best threshold)
    if best_prs and phenotypes:
        n = min(len(best_prs), len(phenotypes))
        pheno = list(phenotypes[:n])
        prs_v = list(best_prs[:n])

        # Split by tertile
        sorted_pheno = sorted(pheno)
        q33 = sorted_pheno[n // 3]
        q66 = sorted_pheno[2 * n // 3]

        prs_low = [prs_v[i] for i in range(n) if pheno[i] <= q33]
        prs_mid = [prs_v[i] for i in range(n) if q33 < pheno[i] <= q66]
        prs_high = [prs_v[i] for i in range(n) if pheno[i] > q66]

        n_bins = max(8, n // 5)
        all_vals = prs_v
        min_v, max_v = min(all_vals), max(all_vals)
        bin_edges = [min_v + (max_v - min_v) * i / n_bins for i in range(n_bins + 1)]

        for vals, label, color in [
            (prs_low, "Low phenotype (T1)", "#4C72B0"),
            (prs_mid, "Mid phenotype (T2)", "#55A868"),
            (prs_high, "High phenotype (T3)", "#C44E52"),
        ]:
            if vals:
                ax2.hist(vals, bins=bin_edges, alpha=0.55, label=label, color=color, edgecolor="white")

        ax2.set_xlabel("Polygenic Risk Score", fontsize=11)
        ax2.set_ylabel("Count", fontsize=11)
        ax2.set_title("PRS Distribution by Phenotype Tertile", fontsize=10)
        ax2.legend(fontsize=9)

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"PRS distribution plot saved to {output_path}")
    return fig
