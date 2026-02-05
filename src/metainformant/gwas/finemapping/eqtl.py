"""eQTL analysis methods for gene expression-variant associations.

This module provides comprehensive eQTL (expression Quantitative Trait Loci)
analysis including cis-eQTL scanning, trans-eQTL analysis, conditional analysis,
and effect size estimation.

Example:
    >>> from metainformant.gwas.finemapping.eqtl import cis_eqtl_scan
    >>> results = cis_eqtl_scan(
    ...     expression_matrix=expr_df,
    ...     genotype_matrix=geno_df,
    ...     gene_positions=gene_pos,
    ...     variant_positions=var_pos,
    ... )
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# cis-eQTL Analysis
# ---------------------------------------------------------------------------


def cis_eqtl_scan(
    expression_matrix: pd.DataFrame,
    genotype_matrix: pd.DataFrame,
    gene_positions: pd.DataFrame,
    variant_positions: pd.DataFrame,
    cis_window: int = 1_000_000,
    maf_threshold: float = 0.05,
) -> pd.DataFrame:
    """Scan for cis-eQTLs within a window around each gene.

    For each gene, tests all variants within the cis_window (default 1Mb)
    for association with gene expression levels.

    Args:
        expression_matrix: Gene expression matrix (genes x samples).
        genotype_matrix: Genotype dosage matrix (variants x samples).
        gene_positions: DataFrame with columns [gene_id, chrom, tss_position].
        variant_positions: DataFrame with columns [variant_id, chrom, position].
        cis_window: Window size (bp) around TSS to test (default 1Mb).
        maf_threshold: Minimum minor allele frequency threshold.

    Returns:
        DataFrame with columns:
            - gene_id: Gene identifier
            - variant_id: Variant identifier
            - beta: Effect size
            - se: Standard error
            - pvalue: P-value from association test
            - distance: Distance from variant to TSS
    """
    logger.info(
        f"Starting cis-eQTL scan: {len(gene_positions)} genes, "
        f"{len(variant_positions)} variants, window={cis_window:,}bp"
    )

    results = []

    # Ensure sample alignment
    common_samples = list(
        set(expression_matrix.columns) & set(genotype_matrix.columns)
    )
    if len(common_samples) < 10:
        logger.warning(f"Only {len(common_samples)} common samples found")
        return pd.DataFrame()

    expr = expression_matrix[common_samples]
    geno = genotype_matrix[common_samples]

    for _, gene_row in gene_positions.iterrows():
        gene_id = gene_row["gene_id"]
        gene_chrom = str(gene_row["chrom"])
        tss = gene_row["tss_position"]

        if gene_id not in expr.index:
            continue

        # Find variants in cis window
        cis_variants = variant_positions[
            (variant_positions["chrom"].astype(str) == gene_chrom)
            & (abs(variant_positions["position"] - tss) <= cis_window)
        ]

        if len(cis_variants) == 0:
            continue

        gene_expr = expr.loc[gene_id].values

        for _, var_row in cis_variants.iterrows():
            var_id = var_row["variant_id"]
            if var_id not in geno.index:
                continue

            var_geno = geno.loc[var_id].values

            # Filter by MAF
            maf = _compute_maf(var_geno)
            if maf < maf_threshold:
                continue

            # Run association
            beta, se, pval = _linear_regression(gene_expr, var_geno)
            distance = int(var_row["position"] - tss)

            results.append({
                "gene_id": gene_id,
                "variant_id": var_id,
                "beta": beta,
                "se": se,
                "pvalue": pval,
                "distance": distance,
                "maf": maf,
            })

    result_df = pd.DataFrame(results)
    logger.info(f"cis-eQTL scan complete: {len(result_df)} tests")
    return result_df


def trans_eqtl_scan(
    expression_matrix: pd.DataFrame,
    genotype_matrix: pd.DataFrame,
    gene_positions: pd.DataFrame,
    variant_positions: pd.DataFrame,
    cis_window: int = 1_000_000,
    maf_threshold: float = 0.05,
    pvalue_threshold: float = 1e-5,
) -> pd.DataFrame:
    """Scan for trans-eQTLs (variants outside cis window).

    Tests variants that are on different chromosomes or beyond the cis window.

    Args:
        expression_matrix: Gene expression matrix (genes x samples).
        genotype_matrix: Genotype dosage matrix (variants x samples).
        gene_positions: DataFrame with columns [gene_id, chrom, tss_position].
        variant_positions: DataFrame with columns [variant_id, chrom, position].
        cis_window: Window to exclude as cis (default 1Mb).
        maf_threshold: Minimum minor allele frequency.
        pvalue_threshold: Only report results below this threshold.

    Returns:
        DataFrame with trans-eQTL results.
    """
    logger.info("Starting trans-eQTL scan")

    results = []
    common_samples = list(
        set(expression_matrix.columns) & set(genotype_matrix.columns)
    )
    expr = expression_matrix[common_samples]
    geno = genotype_matrix[common_samples]

    for _, gene_row in gene_positions.iterrows():
        gene_id = gene_row["gene_id"]
        gene_chrom = str(gene_row["chrom"])
        tss = gene_row["tss_position"]

        if gene_id not in expr.index:
            continue

        gene_expr = expr.loc[gene_id].values

        # Trans = different chromosome OR beyond cis window
        trans_variants = variant_positions[
            (variant_positions["chrom"].astype(str) != gene_chrom)
            | (abs(variant_positions["position"] - tss) > cis_window)
        ]

        for _, var_row in trans_variants.iterrows():
            var_id = var_row["variant_id"]
            if var_id not in geno.index:
                continue

            var_geno = geno.loc[var_id].values
            maf = _compute_maf(var_geno)
            if maf < maf_threshold:
                continue

            beta, se, pval = _linear_regression(gene_expr, var_geno)

            if pval < pvalue_threshold:
                results.append({
                    "gene_id": gene_id,
                    "variant_id": var_id,
                    "beta": beta,
                    "se": se,
                    "pvalue": pval,
                    "is_trans": True,
                    "maf": maf,
                })

    result_df = pd.DataFrame(results)
    logger.info(f"trans-eQTL scan complete: {len(result_df)} significant hits")
    return result_df


def conditional_eqtl(
    expression_vector: np.ndarray,
    genotype_matrix: pd.DataFrame,
    lead_variant_id: str,
    test_variant_ids: list[str] | None = None,
) -> pd.DataFrame:
    """Conditional eQTL analysis controlling for lead variant.

    Tests remaining variants after conditioning on the lead eQTL variant
    to identify independent secondary signals.

    Args:
        expression_vector: Gene expression values.
        genotype_matrix: Genotype dosages (variants x samples).
        lead_variant_id: Lead variant to condition on.
        test_variant_ids: Variants to test (default: all except lead).

    Returns:
        DataFrame with conditional association results.
    """
    logger.info(f"Conditional eQTL analysis, conditioning on {lead_variant_id}")

    if lead_variant_id not in genotype_matrix.index:
        return pd.DataFrame()

    lead_geno = genotype_matrix.loc[lead_variant_id].values

    if test_variant_ids is None:
        test_variant_ids = [
            v for v in genotype_matrix.index if v != lead_variant_id
        ]

    results = []
    for var_id in test_variant_ids:
        if var_id not in genotype_matrix.index:
            continue

        var_geno = genotype_matrix.loc[var_id].values

        # Residualize expression on lead variant
        _, _, lead_resid = _linear_regression_residuals(
            expression_vector, lead_geno
        )

        # Test residuals vs this variant
        beta, se, pval = _linear_regression(lead_resid, var_geno)

        results.append({
            "variant_id": var_id,
            "conditioned_on": lead_variant_id,
            "beta": beta,
            "se": se,
            "pvalue": pval,
        })

    return pd.DataFrame(results)


def eqtl_effect_sizes(
    expression_matrix: pd.DataFrame,
    genotype_matrix: pd.DataFrame,
    associations: pd.DataFrame,
) -> pd.DataFrame:
    """Compute detailed effect sizes for eQTL associations.

    Args:
        expression_matrix: Gene expression (genes x samples).
        genotype_matrix: Genotype dosages (variants x samples).
        associations: DataFrame with gene_id, variant_id columns.

    Returns:
        DataFrame with effect size details:
            - beta: Effect size
            - se: Standard error
            - r_squared: Variance explained
            - mean_expr_by_genotype: Mean expression per genotype class
    """
    common_samples = list(
        set(expression_matrix.columns) & set(genotype_matrix.columns)
    )
    expr = expression_matrix[common_samples]
    geno = genotype_matrix[common_samples]

    results = []
    for _, assoc in associations.iterrows():
        gene_id = assoc["gene_id"]
        var_id = assoc["variant_id"]

        if gene_id not in expr.index or var_id not in geno.index:
            continue

        gene_expr = expr.loc[gene_id].values
        var_geno = geno.loc[var_id].values

        beta, se, pval = _linear_regression(gene_expr, var_geno)
        r2 = _compute_r_squared(gene_expr, var_geno)

        # Mean expression by genotype (0, 1, 2)
        means = {}
        for gt in [0, 1, 2]:
            mask = np.round(var_geno) == gt
            if mask.sum() > 0:
                means[f"mean_expr_gt{gt}"] = np.mean(gene_expr[mask])

        results.append({
            "gene_id": gene_id,
            "variant_id": var_id,
            "beta": beta,
            "se": se,
            "pvalue": pval,
            "r_squared": r2,
            **means,
        })

    return pd.DataFrame(results)


def eqtl_summary_stats(
    cis_results: pd.DataFrame,
    fdr_threshold: float = 0.05,
) -> dict[str, Any]:
    """Generate eQTL summary statistics.

    Args:
        cis_results: Results from cis_eqtl_scan.
        fdr_threshold: FDR threshold for significant eQTLs.

    Returns:
        Dictionary with summary statistics.
    """
    if len(cis_results) == 0:
        return {"n_tests": 0, "n_egenes": 0, "n_eqtls": 0}

    # FDR correction
    from scipy import stats
    pvals = cis_results["pvalue"].values
    _, fdr_pvals = stats.false_discovery_control(pvals, method="bh"), None
    
    # Simple Benjamini-Hochberg
    n = len(pvals)
    sorted_idx = np.argsort(pvals)
    sorted_pvals = pvals[sorted_idx]
    fdr = np.zeros(n)
    for i, p in enumerate(sorted_pvals):
        fdr[sorted_idx[i]] = p * n / (i + 1)
    fdr = np.minimum.accumulate(fdr[::-1])[::-1]
    
    significant = fdr < fdr_threshold

    n_egenes = cis_results.loc[significant, "gene_id"].nunique()
    n_eqtls = significant.sum()

    return {
        "n_tests": len(cis_results),
        "n_egenes": n_egenes,
        "n_eqtls": int(n_eqtls),
        "fdr_threshold": fdr_threshold,
        "mean_effect_size": float(cis_results["beta"].abs().mean()),
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _compute_maf(genotypes: np.ndarray) -> float:
    """Compute minor allele frequency from dosages."""
    valid = ~np.isnan(genotypes)
    if valid.sum() == 0:
        return 0.0
    af = np.mean(genotypes[valid]) / 2
    return min(af, 1 - af)


def _linear_regression(
    y: np.ndarray, x: np.ndarray
) -> tuple[float, float, float]:
    """Simple linear regression returning beta, SE, p-value."""
    valid = ~(np.isnan(y) | np.isnan(x))
    if valid.sum() < 3:
        return np.nan, np.nan, 1.0

    y_v, x_v = y[valid], x[valid]
    n = len(y_v)

    # Center variables
    x_c = x_v - np.mean(x_v)
    y_c = y_v - np.mean(y_v)

    # Calculate beta
    ss_x = np.sum(x_c**2)
    if ss_x == 0:
        return np.nan, np.nan, 1.0

    beta = np.sum(x_c * y_c) / ss_x

    # Residuals and SE
    y_pred = np.mean(y_v) + beta * x_c
    resid = y_v - y_pred
    mse = np.sum(resid**2) / (n - 2) if n > 2 else 1.0
    se = np.sqrt(mse / ss_x) if ss_x > 0 else np.nan

    # T-statistic and p-value
    if se > 0 and not np.isnan(se):
        t_stat = beta / se
        from scipy import stats
        pval = 2 * stats.t.sf(abs(t_stat), n - 2)
    else:
        pval = 1.0

    return float(beta), float(se), float(pval)


def _linear_regression_residuals(
    y: np.ndarray, x: np.ndarray
) -> tuple[float, float, np.ndarray]:
    """Linear regression returning beta, SE, and residuals."""
    valid = ~(np.isnan(y) | np.isnan(x))
    y_v, x_v = y.copy(), x.copy()
    y_v[~valid] = np.nan
    x_v[~valid] = np.nan

    beta, se, _ = _linear_regression(y, x)

    if np.isnan(beta):
        return beta, se, y

    residuals = y - (np.nanmean(y) + beta * (x - np.nanmean(x)))
    return beta, se, residuals


def _compute_r_squared(y: np.ndarray, x: np.ndarray) -> float:
    """Compute R-squared between response and predictor."""
    valid = ~(np.isnan(y) | np.isnan(x))
    if valid.sum() < 3:
        return 0.0

    y_v, x_v = y[valid], x[valid]
    correlation = np.corrcoef(y_v, x_v)[0, 1]
    return float(correlation**2) if not np.isnan(correlation) else 0.0
