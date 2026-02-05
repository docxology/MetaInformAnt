"""Differential expression analysis for single-cell data.

Provides statistical testing for identifying genes that are differentially
expressed between groups of cells. Supports Wilcoxon rank-sum test, Welch's
t-test, and pseudobulk aggregation for multi-sample experimental designs.
Also includes fold change computation, volcano plot data preparation, and
gene set activity scoring.

Statistical methods:
    - Wilcoxon rank-sum: Non-parametric test comparing expression ranks
      between two groups. Robust to non-normality and outliers.
    - t-test: Welch's t-test for unequal variances. Assumes approximate
      normality (reasonable for log-normalized data).
    - Pseudobulk: Aggregates single cells by sample, then performs
      sample-level testing. More appropriate for multi-donor designs
      where cells within a donor are not independent.
"""

from __future__ import annotations

import math
import random
from collections import Counter, defaultdict
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional scientific dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None  # type: ignore[assignment]


def differential_expression(
    expression_matrix: Any,
    groups: list[int],
    gene_names: list[str],
    method: str = "wilcoxon",
    min_cells: int = 3,
    min_log2fc: float = 0.0,
) -> list[dict]:
    """Test for differential expression between two cell groups.

    Performs gene-by-gene statistical testing to identify genes whose
    expression differs significantly between cells in group 0 and group 1.
    Returns results sorted by adjusted p-value.

    Args:
        expression_matrix: Expression matrix (cells x genes). Accepts lists
            of lists, numpy arrays, or scipy sparse matrices.
        groups: Group assignment for each cell. Must contain exactly two
            unique values (typically 0 and 1).
        gene_names: Gene names corresponding to columns in the expression
            matrix.
        method: Statistical method. One of "wilcoxon" (Wilcoxon rank-sum
            test), "t_test" (Welch's t-test), or "pseudobulk" (not valid
            here; use pseudobulk_de instead).
        min_cells: Minimum number of cells per group required to test a
            gene. Genes where either group has fewer expressing cells are
            skipped.
        min_log2fc: Minimum absolute log2 fold change to include in results.

    Returns:
        Sorted list of dictionaries, each containing:
            - gene: str gene name.
            - log2fc: float log2 fold change (group 1 vs group 0).
            - p_value: float raw p-value from the statistical test.
            - adjusted_p: float Benjamini-Hochberg adjusted p-value.
            - pct_group1: float percentage of group-1 cells expressing the gene.
            - pct_group2: float percentage of group-0 cells expressing the gene.
            - mean_group1: float mean expression in group 1.
            - mean_group2: float mean expression in group 0.

    Raises:
        ValueError: If groups does not contain exactly two unique values,
            or dimensions are inconsistent.
    """
    valid_methods = ("wilcoxon", "t_test")
    if method not in valid_methods:
        raise ValueError(f"Invalid method '{method}'. Must be one of {valid_methods}")

    matrix = _to_list_matrix(expression_matrix)
    n_cells = len(matrix)
    n_genes = len(matrix[0]) if n_cells > 0 else 0

    if len(groups) != n_cells:
        raise ValueError(f"groups length ({len(groups)}) must match expression matrix " f"rows ({n_cells})")
    if len(gene_names) != n_genes:
        raise ValueError(f"gene_names length ({len(gene_names)}) must match expression " f"matrix columns ({n_genes})")

    unique_groups = sorted(set(groups))
    if len(unique_groups) != 2:
        raise ValueError(f"groups must contain exactly 2 unique values, got {len(unique_groups)}: " f"{unique_groups}")

    g1_val, g2_val = unique_groups
    g1_indices = [i for i, g in enumerate(groups) if g == g1_val]
    g2_indices = [i for i, g in enumerate(groups) if g == g2_val]

    logger.info(f"Running DE analysis ({method}): {len(g1_indices)} vs {len(g2_indices)} cells, " f"{n_genes} genes")

    results: list[dict] = []
    raw_p_values: list[float] = []

    for gene_idx in range(n_genes):
        expr_g1 = [matrix[i][gene_idx] for i in g1_indices]
        expr_g2 = [matrix[i][gene_idx] for i in g2_indices]

        # Skip genes with too few expressing cells
        n_expr_g1 = sum(1 for v in expr_g1 if v > 0)
        n_expr_g2 = sum(1 for v in expr_g2 if v > 0)
        if n_expr_g1 < min_cells and n_expr_g2 < min_cells:
            continue

        mean_g1 = sum(expr_g1) / len(expr_g1) if expr_g1 else 0.0
        mean_g2 = sum(expr_g2) / len(expr_g2) if expr_g2 else 0.0

        log2fc = compute_log_fold_change(mean_g1, mean_g2)

        if abs(log2fc) < min_log2fc:
            continue

        # Statistical test
        if method == "wilcoxon":
            p_value = _wilcoxon_rank_sum(expr_g1, expr_g2)
        else:  # t_test
            p_value = _welch_t_test(expr_g1, expr_g2)

        pct_g1 = (n_expr_g1 / len(expr_g1) * 100) if expr_g1 else 0.0
        pct_g2 = (n_expr_g2 / len(expr_g2) * 100) if expr_g2 else 0.0

        results.append(
            {
                "gene": gene_names[gene_idx],
                "log2fc": log2fc,
                "p_value": p_value,
                "adjusted_p": 0.0,  # placeholder, filled below
                "pct_group1": pct_g1,
                "pct_group2": pct_g2,
                "mean_group1": mean_g1,
                "mean_group2": mean_g2,
            }
        )
        raw_p_values.append(p_value)

    # Benjamini-Hochberg correction
    adjusted = _benjamini_hochberg(raw_p_values)
    for i, res in enumerate(results):
        res["adjusted_p"] = adjusted[i]

    # Sort by adjusted p-value
    results.sort(key=lambda x: x["adjusted_p"])

    logger.info(
        f"DE analysis complete: {len(results)} genes tested, "
        f"{sum(1 for r in results if r['adjusted_p'] < 0.05)} significant (FDR < 0.05)"
    )

    return results


def pseudobulk_de(
    expression_matrix: Any,
    cell_labels: list[str],
    sample_labels: list[str],
    groups: list[int],
    gene_names: list[str] | None = None,
    min_cells_per_sample: int = 5,
) -> list[dict]:
    """Pseudobulk differential expression for multi-sample designs.

    Aggregates single-cell expression by sample (summing counts per gene
    within each sample), then performs sample-level statistical testing
    between groups. This approach accounts for the hierarchical structure
    of multi-donor/multi-sample experiments where cells within a sample
    are not independent.

    Args:
        expression_matrix: Expression matrix (cells x genes).
        cell_labels: Cell type labels for each cell. Only cells of the
            same type are aggregated (pass a single type or filter first).
        sample_labels: Sample identifier for each cell.
        groups: Group assignment for each sample (must be consistent
            within a sample). Exactly two unique values required.
        gene_names: Gene names for columns. Auto-generated if None.
        min_cells_per_sample: Minimum cells per sample required for
            inclusion in the pseudobulk analysis.

    Returns:
        Sorted list of dictionaries with DE results (same format as
        differential_expression).

    Raises:
        ValueError: If input dimensions are inconsistent or groups
            do not contain exactly two unique values.
    """
    matrix = _to_list_matrix(expression_matrix)
    n_cells = len(matrix)
    n_genes = len(matrix[0]) if n_cells > 0 else 0

    if len(cell_labels) != n_cells:
        raise ValueError(f"cell_labels length ({len(cell_labels)}) must match rows ({n_cells})")
    if len(sample_labels) != n_cells:
        raise ValueError(f"sample_labels length ({len(sample_labels)}) must match rows ({n_cells})")
    if len(groups) != n_cells:
        raise ValueError(f"groups length ({len(groups)}) must match rows ({n_cells})")

    if gene_names is None:
        gene_names = [f"gene_{i}" for i in range(n_genes)]
    elif len(gene_names) != n_genes:
        raise ValueError(f"gene_names length ({len(gene_names)}) must match columns ({n_genes})")

    logger.info(f"Running pseudobulk DE: {n_cells} cells, {n_genes} genes")

    # Aggregate expression per sample
    samples = sorted(set(sample_labels))
    sample_sums: dict[str, list[float]] = {}
    sample_counts: dict[str, int] = {}
    sample_groups: dict[str, int] = {}

    for sample in samples:
        sample_indices = [i for i, s in enumerate(sample_labels) if s == sample]
        if len(sample_indices) < min_cells_per_sample:
            logger.debug(
                f"Skipping sample '{sample}': only {len(sample_indices)} cells " f"(min={min_cells_per_sample})"
            )
            continue

        # Sum expression across cells in this sample
        gene_sums = [0.0] * n_genes
        for idx in sample_indices:
            for j in range(n_genes):
                gene_sums[j] += matrix[idx][j]

        sample_sums[sample] = gene_sums
        sample_counts[sample] = len(sample_indices)

        # Determine group for this sample (majority vote)
        sample_group_vals = [groups[i] for i in sample_indices]
        group_counts = Counter(sample_group_vals)
        sample_groups[sample] = group_counts.most_common(1)[0][0]

    # Check we have two groups
    unique_groups = sorted(set(sample_groups.values()))
    if len(unique_groups) != 2:
        raise ValueError(f"After aggregation, need exactly 2 groups, got {len(unique_groups)}: " f"{unique_groups}")

    g1_val, g2_val = unique_groups
    g1_samples = [s for s, g in sample_groups.items() if g == g1_val]
    g2_samples = [s for s, g in sample_groups.items() if g == g2_val]

    logger.info(
        f"Pseudobulk aggregated: {len(g1_samples)} samples in group 1, " f"{len(g2_samples)} samples in group 2"
    )

    # Run DE on pseudobulk expression
    results: list[dict] = []
    raw_p_values: list[float] = []

    for gene_idx in range(n_genes):
        vals_g1 = [sample_sums[s][gene_idx] for s in g1_samples]
        vals_g2 = [sample_sums[s][gene_idx] for s in g2_samples]

        mean_g1 = sum(vals_g1) / len(vals_g1) if vals_g1 else 0.0
        mean_g2 = sum(vals_g2) / len(vals_g2) if vals_g2 else 0.0

        log2fc = compute_log_fold_change(mean_g1, mean_g2)

        # Welch's t-test on sample-level aggregates
        p_value = _welch_t_test(vals_g1, vals_g2)

        pct_g1 = (sum(1 for v in vals_g1 if v > 0) / len(vals_g1) * 100) if vals_g1 else 0.0
        pct_g2 = (sum(1 for v in vals_g2 if v > 0) / len(vals_g2) * 100) if vals_g2 else 0.0

        results.append(
            {
                "gene": gene_names[gene_idx],
                "log2fc": log2fc,
                "p_value": p_value,
                "adjusted_p": 0.0,
                "pct_group1": pct_g1,
                "pct_group2": pct_g2,
                "mean_group1": mean_g1,
                "mean_group2": mean_g2,
            }
        )
        raw_p_values.append(p_value)

    # Benjamini-Hochberg correction
    adjusted = _benjamini_hochberg(raw_p_values)
    for i, res in enumerate(results):
        res["adjusted_p"] = adjusted[i]

    results.sort(key=lambda x: x["adjusted_p"])

    logger.info(
        f"Pseudobulk DE complete: {len(results)} genes, "
        f"{sum(1 for r in results if r['adjusted_p'] < 0.05)} significant (FDR < 0.05)"
    )

    return results


def compute_log_fold_change(
    mean_a: float,
    mean_b: float,
    pseudocount: float = 1.0,
) -> float:
    """Compute log2 fold change between two mean expression values.

    Uses a pseudocount to avoid division by zero and log of zero.

    Args:
        mean_a: Mean expression in group A (numerator).
        mean_b: Mean expression in group B (denominator).
        pseudocount: Added to both values before computing the ratio.
            Default 1.0.

    Returns:
        Log2 fold change: log2((mean_a + pseudocount) / (mean_b + pseudocount)).
    """
    return math.log2((mean_a + pseudocount) / (mean_b + pseudocount))


def volcano_data(
    de_results: list[dict],
    fc_threshold: float = 1.0,
    p_threshold: float = 0.05,
) -> dict:
    """Prepare volcano plot data with significance classification.

    Classifies each gene as "up" (significantly upregulated), "down"
    (significantly downregulated), or "ns" (not significant) based on
    fold change and p-value thresholds.

    Args:
        de_results: List of DE result dictionaries from
            differential_expression or pseudobulk_de.
        fc_threshold: Minimum absolute log2 fold change for significance.
        p_threshold: Maximum adjusted p-value for significance.

    Returns:
        Dictionary with keys:
            - genes: list[str] of gene names.
            - log2fc: list[float] of log2 fold change values (x-axis).
            - neg_log10_p: list[float] of -log10(adjusted_p) values (y-axis).
            - classification: list[str] of "up", "down", or "ns" per gene.
            - n_up: int count of upregulated genes.
            - n_down: int count of downregulated genes.
            - n_ns: int count of non-significant genes.
    """
    genes: list[str] = []
    log2fc_values: list[float] = []
    neg_log10_p: list[float] = []
    classification: list[str] = []

    for res in de_results:
        gene = res.get("gene", "unknown")
        fc = res.get("log2fc", 0.0)
        adj_p = res.get("adjusted_p", 1.0)

        genes.append(gene)
        log2fc_values.append(fc)

        # Compute -log10(p), clamping to avoid log(0)
        safe_p = max(adj_p, 1e-300)
        neg_log10 = -math.log10(safe_p)
        neg_log10_p.append(neg_log10)

        # Classify
        if adj_p < p_threshold and fc > fc_threshold:
            classification.append("up")
        elif adj_p < p_threshold and fc < -fc_threshold:
            classification.append("down")
        else:
            classification.append("ns")

    n_up = classification.count("up")
    n_down = classification.count("down")
    n_ns = classification.count("ns")

    logger.info(f"Volcano data prepared: {n_up} up, {n_down} down, {n_ns} not significant")

    return {
        "genes": genes,
        "log2fc": log2fc_values,
        "neg_log10_p": neg_log10_p,
        "classification": classification,
        "n_up": n_up,
        "n_down": n_down,
        "n_ns": n_ns,
    }


def gene_set_scoring(
    expression_matrix: Any,
    gene_sets: dict[str, list[str]],
    gene_names: list[str],
    method: str = "mean",
    n_background: int = 50,
    seed: int | None = None,
) -> dict:
    """Score gene set activity per cell.

    Computes a score for each gene set in each cell, reflecting how
    strongly the cell expresses that gene set relative to background.

    Args:
        expression_matrix: Expression matrix (cells x genes).
        gene_sets: Dictionary mapping gene set names to lists of gene names.
        gene_names: Gene names for columns in expression_matrix.
        method: Scoring method. "mean" computes mean expression of gene set
            genes minus background. "sum" uses the total expression.
        n_background: Number of random background genes for normalization
            (used only with "mean" method).
        seed: Random seed for reproducible background selection.

    Returns:
        Dictionary with keys:
            - scores: dict mapping gene set names to lists of per-cell
              scores (one float per cell).
            - n_cells: int number of cells scored.
            - n_gene_sets: int number of gene sets.
            - gene_set_sizes: dict mapping gene set names to the number
              of genes found in the expression matrix.
    """
    valid_methods = ("mean", "sum")
    if method not in valid_methods:
        raise ValueError(f"Invalid method '{method}'. Must be one of {valid_methods}")

    matrix = _to_list_matrix(expression_matrix)
    n_cells = len(matrix)
    n_genes = len(matrix[0]) if n_cells > 0 else 0

    if len(gene_names) != n_genes:
        raise ValueError(f"gene_names length ({len(gene_names)}) must match " f"columns ({n_genes})")

    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    logger.info(f"Scoring {len(gene_sets)} gene sets across {n_cells} cells (method={method})")

    rng = random.Random(seed)
    all_indices = list(range(n_genes))

    scores: dict[str, list[float]] = {}
    gene_set_sizes: dict[str, int] = {}

    for gs_name, gs_genes in gene_sets.items():
        gs_indices = [gene_to_idx[g] for g in gs_genes if g in gene_to_idx]
        gene_set_sizes[gs_name] = len(gs_indices)

        if not gs_indices:
            scores[gs_name] = [0.0] * n_cells
            continue

        # Build background gene set (excluding current gene set)
        non_gs_indices = [i for i in all_indices if i not in set(gs_indices)]
        n_bg = min(n_background, len(non_gs_indices))
        bg_indices = rng.sample(non_gs_indices, n_bg) if n_bg > 0 else []

        cell_scores: list[float] = []
        for cell_idx in range(n_cells):
            row = matrix[cell_idx]

            if method == "mean":
                gs_mean = sum(row[j] for j in gs_indices) / len(gs_indices)
                bg_mean = sum(row[j] for j in bg_indices) / len(bg_indices) if bg_indices else 0.0
                cell_scores.append(gs_mean - bg_mean)
            else:  # sum
                gs_sum = sum(row[j] for j in gs_indices)
                cell_scores.append(gs_sum)

        scores[gs_name] = cell_scores

    logger.info(f"Gene set scoring complete: {len(gene_sets)} sets scored")

    return {
        "scores": scores,
        "n_cells": n_cells,
        "n_gene_sets": len(gene_sets),
        "gene_set_sizes": gene_set_sizes,
    }


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _to_list_matrix(data: Any) -> list[list[float]]:
    """Convert various matrix types to list of lists of floats.

    Args:
        data: Expression matrix. Supports list[list[float]], numpy arrays,
            scipy sparse matrices, or objects with a .toarray() method.

    Returns:
        List of lists of floats.
    """
    if isinstance(data, list):
        return data

    if HAS_NUMPY and isinstance(data, np.ndarray):
        return data.tolist()

    if hasattr(data, "toarray"):
        return data.toarray().tolist()

    if hasattr(data, "values"):
        return data.values.tolist()

    raise TypeError(f"Unsupported expression matrix type: {type(data)}")


def _wilcoxon_rank_sum(a: list[float], b: list[float]) -> float:
    """Perform Wilcoxon rank-sum test (Mann-Whitney U test).

    Computes the U statistic and derives a p-value using normal
    approximation for sample sizes > 20, otherwise returns an
    approximate p-value.

    Args:
        a: Values from group A.
        b: Values from group B.

    Returns:
        Two-sided p-value.
    """
    if HAS_SCIPY:
        try:
            stat, p_val = scipy_stats.ranksums(a, b)
            return float(p_val)
        except (ValueError, TypeError):
            return 1.0

    # Pure Python fallback: Mann-Whitney U with normal approximation
    n_a = len(a)
    n_b = len(b)
    if n_a == 0 or n_b == 0:
        return 1.0

    # Combine and rank
    combined = [(v, 0) for v in a] + [(v, 1) for v in b]
    combined.sort(key=lambda x: x[0])

    # Assign ranks (handle ties with average rank)
    ranks: list[float] = [0.0] * len(combined)
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg_rank = (i + j + 1) / 2.0  # 1-indexed average
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j

    # Sum of ranks for group A
    rank_sum_a = sum(ranks[i] for i in range(len(combined)) if combined[i][1] == 0)

    # U statistic
    u_a = rank_sum_a - n_a * (n_a + 1) / 2.0

    # Normal approximation
    mu = n_a * n_b / 2.0
    sigma = math.sqrt(n_a * n_b * (n_a + n_b + 1) / 12.0)

    if sigma == 0:
        return 1.0

    z = (u_a - mu) / sigma

    # Two-sided p-value from standard normal
    p_value = 2.0 * _standard_normal_cdf(-abs(z))
    return p_value


def _welch_t_test(a: list[float], b: list[float]) -> float:
    """Perform Welch's t-test for unequal variances.

    Args:
        a: Values from group A.
        b: Values from group B.

    Returns:
        Two-sided p-value.
    """
    if HAS_SCIPY:
        try:
            stat, p_val = scipy_stats.ttest_ind(a, b, equal_var=False)
            return float(p_val) if not math.isnan(p_val) else 1.0
        except (ValueError, TypeError):
            return 1.0

    # Pure Python fallback
    n_a = len(a)
    n_b = len(b)
    if n_a < 2 or n_b < 2:
        return 1.0

    mean_a = sum(a) / n_a
    mean_b = sum(b) / n_b
    var_a = sum((x - mean_a) ** 2 for x in a) / (n_a - 1)
    var_b = sum((x - mean_b) ** 2 for x in b) / (n_b - 1)

    se = math.sqrt(var_a / n_a + var_b / n_b)
    if se == 0:
        return 1.0

    t_stat = (mean_a - mean_b) / se

    # Welch-Satterthwaite degrees of freedom
    num = (var_a / n_a + var_b / n_b) ** 2
    denom = (var_a / n_a) ** 2 / (n_a - 1) + (var_b / n_b) ** 2 / (n_b - 1)
    df = num / denom if denom > 0 else 1.0

    # Approximate p-value using normal distribution for large df
    p_value = 2.0 * _standard_normal_cdf(-abs(t_stat))
    return p_value


def _standard_normal_cdf(x: float) -> float:
    """Approximate CDF of the standard normal distribution.

    Uses the Abramowitz and Stegun approximation (formula 7.1.26).

    Args:
        x: Value at which to evaluate the CDF.

    Returns:
        Approximate P(Z <= x) for Z ~ N(0,1).
    """
    # Error function approximation
    sign = 1 if x >= 0 else -1
    x = abs(x) / math.sqrt(2.0)

    # Rational approximation of erf
    t = 1.0 / (1.0 + 0.3275911 * x)
    a1, a2, a3, a4, a5 = 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429
    erf = 1.0 - (a1 * t + a2 * t**2 + a3 * t**3 + a4 * t**4 + a5 * t**5) * math.exp(-x * x)

    return 0.5 * (1.0 + sign * erf)


def _benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Apply Benjamini-Hochberg FDR correction to a list of p-values.

    Args:
        p_values: Raw p-values.

    Returns:
        Adjusted p-values (same length and order as input).
    """
    n = len(p_values)
    if n == 0:
        return []

    # Sort p-values and track original indices
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])

    adjusted = [0.0] * n
    min_so_far = 1.0

    for rank_idx in range(n - 1, -1, -1):
        orig_idx, p_val = indexed[rank_idx]
        rank = rank_idx + 1  # 1-indexed rank
        adj = p_val * n / rank
        adj = min(adj, min_so_far)
        adj = min(adj, 1.0)
        adjusted[orig_idx] = adj
        min_so_far = adj

    return adjusted


__all__ = [
    "compute_log_fold_change",
    "differential_expression",
    "gene_set_scoring",
    "pseudobulk_de",
    "volcano_data",
]
