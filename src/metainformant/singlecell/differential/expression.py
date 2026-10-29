"""Differential expression analysis for single-cell data.

Provides statistical testing for identifying genes that are differentially
expressed between groups of cells. Supports Wilcoxon rank-sum test, Welch's
t-test, and pseudobulk aggregation for multi-sample experimental designs.
Also includes fold change computation, volcano plot data preparation, and
gene set activity scoring.

Statistical methods:
    - Wilcoxon rank-sum: Non-parametric test comparing expression ranks
      between two groups. Tie-aware (exact or tie-corrected asymptotic
      Mann-Whitney U), which matters for massively tied single-cell
      counts. Robust to non-normality and outliers.
    - t-test: Welch's t-test for unequal variances. Assumes approximate
      normality (reasonable for log-normalized data).
    - Pseudobulk: Aggregates single cells by sample, CPM-normalizes the
      per-sample sums, then performs sample-level testing. More
      appropriate for multi-donor designs where cells within a donor are
      not independent.

Undefined tests: a gene whose test cannot be computed (an empty group, fewer than two values for a t-test,
constant values, non-finite values) gets ``p_value = NaN`` and ``adjusted_p = NaN``, and its
``test_status`` says why. No p-value is ever substituted for a test that did not run: such genes take no
part in the Benjamini-Hochberg correction (it runs over the genes that have a p-value) and sort after
every gene that has one.
"""

from __future__ import annotations

import math
import random
from collections.abc import Sequence
from typing import Any, NamedTuple

from metainformant.core.utils.errors import ValidationError
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional scientific dependencies
np: Any  # None when numpy is unavailable (optional dependency)
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None

scipy_stats: Any  # None when scipy is unavailable (optional dependency)
try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None


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
    expression differs significantly between the cells of the two groups.
    Returns results sorted by adjusted p-value.

    The two group labels are ordered: "group 1" below is the cells carrying the smaller label (0 for the
    usual 0/1 labels) and "group 2" the cells carrying the larger one. ``log2fc`` is group 1 over group 2,
    so a gene that is higher in the cells labelled 1 has a NEGATIVE ``log2fc``.

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
        min_cells: Minimum number of expressing cells required to test a
            gene. A gene is skipped only when BOTH groups have fewer
            expressing cells than this.
        min_log2fc: Minimum absolute log2 fold change to include in results.
            Applied AFTER multiple-testing correction: Benjamini-Hochberg
            adjusted p-values are computed over the full tested family
            (every gene passing min_cells that has a p-value), so this
            effect-size filter never shrinks the BH family or inflates
            adjusted significance.

    Returns:
        List of dictionaries sorted by adjusted p-value (genes without a p-value last), each containing:
            - gene: str gene name.
            - log2fc: float log2 fold change, log2((mean_group1 + 1) / (mean_group2 + 1)).
            - p_value: float raw two-sided p-value; NaN when the test is undefined for this gene.
            - adjusted_p: float Benjamini-Hochberg adjusted p-value; NaN when ``p_value`` is NaN.
            - pct_group1: float percentage of group-1 cells expressing the gene.
            - pct_group2: float percentage of group-2 cells expressing the gene.
            - mean_group1: float mean expression in group 1.
            - mean_group2: float mean expression in group 2.
            - test_status: ``"ok"``, or why the test is undefined (for example ``"zero_variance: ..."``).

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

        # Statistical test. A gene whose test is undefined keeps NaN and a status; it is not part of the BH family.
        if method == "wilcoxon":
            outcome = _wilcoxon_rank_sum(expr_g1, expr_g2)
        else:  # t_test
            outcome = _welch_t_test(expr_g1, expr_g2)

        pct_g1 = (n_expr_g1 / len(expr_g1) * 100) if expr_g1 else 0.0
        pct_g2 = (n_expr_g2 / len(expr_g2) * 100) if expr_g2 else 0.0

        results.append(
            {
                "gene": gene_names[gene_idx],
                "log2fc": log2fc,
                "p_value": outcome.p_value,
                "adjusted_p": math.nan,  # filled below
                "pct_group1": pct_g1,
                "pct_group2": pct_g2,
                "mean_group1": mean_g1,
                "mean_group2": mean_g2,
                "test_status": outcome.status,
            }
        )
        raw_p_values.append(outcome.p_value)

    # Benjamini-Hochberg correction over the FULL family of genes that have a p-value: the
    # min_log2fc effect-size filter is applied only afterwards, so it can
    # never shrink the BH family or inflate adjusted significance.
    adjusted = _benjamini_hochberg(raw_p_values)
    for i, res in enumerate(results):
        res["adjusted_p"] = adjusted[i]
    _warn_about_undefined_tests(results)

    results = [res for res in results if abs(res["log2fc"]) >= min_log2fc]

    # Sort by adjusted p-value; genes without one go last
    results.sort(key=_adjusted_p_sort_key)

    logger.info(
        f"DE analysis complete: {len(raw_p_values)} genes tested, {len(results)} reported "
        f"(min_log2fc={min_log2fc}), "
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
    within each sample), CPM-normalizes the per-sample sums so that
    library-size differences between samples do not masquerade as
    differential expression, then performs sample-level statistical
    testing between groups. This approach accounts for the hierarchical
    structure of multi-donor/multi-sample experiments where cells within
    a sample are not independent.

    Args:
        expression_matrix: Expression matrix (cells x genes).
        cell_labels: Cell type labels for each cell. Validated for length
            but not used to filter; pre-filter the matrix to a single cell
            type before calling.
        sample_labels: Sample identifier for each cell.
        min_cells_per_sample: Minimum cells per sample required for
            inclusion in the pseudobulk analysis.

    Returns:
        Sorted list of dictionaries with DE results (same format as
        differential_expression, including ``test_status`` and the NaN convention for undefined tests: with
        a single sample in a group there is no t-test, so every gene has ``p_value = NaN``) plus a
        "normalization" key recording the sample-level normalization applied ("cpm").

    Raises:
        ValidationError: If a sample's cells are assigned to more than
            one group (group labels must be consistent within a sample).
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

    # Group labels must be consistent within each sample: a sample whose
    # cells span more than one group has no well-defined pseudobulk
    # identity, so fail closed instead of silently majority-voting.
    sample_group_values: dict[str, set[int]] = {}
    for idx, sample in enumerate(sample_labels):
        sample_group_values.setdefault(sample, set()).add(groups[idx])
    for sample in sorted(sample_group_values):
        values = sorted(sample_group_values[sample])
        if len(values) > 1:
            raise ValidationError(
                f"Sample '{sample}' has cells assigned to multiple groups {values}; "
                "group assignment must be consistent within a sample"
            )

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

        # Group is validated consistent across this sample's cells above
        sample_groups[sample] = groups[sample_indices[0]]

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

    # CPM-normalize the pseudobulk sums per sample so that library-size
    # (total count) differences between samples do not masquerade as
    # differential expression. Samples with zero total counts stay at zero.
    normalized_sums: dict[str, list[float]] = {}
    for sample, gene_sums in sample_sums.items():
        total = sum(gene_sums)
        scale = 1e6 / total if total > 0 else 0.0
        normalized_sums[sample] = [value * scale for value in gene_sums]

    # Run DE on pseudobulk expression
    results: list[dict] = []
    raw_p_values: list[float] = []

    for gene_idx in range(n_genes):
        vals_g1 = [normalized_sums[s][gene_idx] for s in g1_samples]
        vals_g2 = [normalized_sums[s][gene_idx] for s in g2_samples]

        mean_g1 = sum(vals_g1) / len(vals_g1) if vals_g1 else 0.0
        mean_g2 = sum(vals_g2) / len(vals_g2) if vals_g2 else 0.0

        log2fc = compute_log_fold_change(mean_g1, mean_g2)

        # Welch's t-test on CPM-normalized sample-level aggregates (NaN and a status when it is undefined,
        # for example with a single sample in a group)
        outcome = _welch_t_test(vals_g1, vals_g2)

        pct_g1 = (sum(1 for v in vals_g1 if v > 0) / len(vals_g1) * 100) if vals_g1 else 0.0
        pct_g2 = (sum(1 for v in vals_g2 if v > 0) / len(vals_g2) * 100) if vals_g2 else 0.0

        results.append(
            {
                "gene": gene_names[gene_idx],
                "log2fc": log2fc,
                "p_value": outcome.p_value,
                "adjusted_p": math.nan,  # filled below
                "pct_group1": pct_g1,
                "pct_group2": pct_g2,
                "mean_group1": mean_g1,
                "mean_group2": mean_g2,
                "normalization": "cpm",
                "test_status": outcome.status,
            }
        )
        raw_p_values.append(outcome.p_value)

    # Benjamini-Hochberg correction over the genes that have a p-value
    adjusted = _benjamini_hochberg(raw_p_values)
    for i, res in enumerate(results):
        res["adjusted_p"] = adjusted[i]
    _warn_about_undefined_tests(results)

    results.sort(key=_adjusted_p_sort_key)

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
    fold change and p-value thresholds. A gene without an adjusted p-value (NaN, because its test was
    undefined, or no ``adjusted_p`` key at all) cannot be called significant: it is "ns" with a NaN
    ``neg_log10_p``, never a made-up p-value of 1.

    Args:
        de_results: List of DE result dictionaries from
            differential_expression or pseudobulk_de.
        fc_threshold: Minimum absolute log2 fold change for significance.
        p_threshold: Maximum adjusted p-value for significance.

    Returns:
        Dictionary with keys:
            - genes: list[str] of gene names.
            - log2fc: list[float] of log2 fold change values (x-axis).
            - neg_log10_p: list[float] of -log10(adjusted_p) values (y-axis); NaN for genes without an
              adjusted p-value.
            - classification: list[str] of "up", "down", or "ns" per gene.
            - n_up: int count of upregulated genes.
            - n_down: int count of downregulated genes.
            - n_ns: int count of non-significant genes (including those without an adjusted p-value).
    """
    genes: list[str] = []
    log2fc_values: list[float] = []
    neg_log10_p: list[float] = []
    classification: list[str] = []

    for res in de_results:
        gene = res.get("gene", "unknown")
        fc = res.get("log2fc", 0.0)
        adj_p = res.get("adjusted_p", math.nan)

        genes.append(gene)
        log2fc_values.append(fc)

        # Compute -log10(p), clamping to avoid log(0); no p-value, no y value
        if math.isnan(adj_p):
            neg_log10_p.append(math.nan)
        else:
            neg_log10_p.append(-math.log10(max(adj_p, 1e-300)))

        # Classify (a NaN p-value compares False, so it can only be "ns")
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
              scores (one float per cell). A gene set with none of its
              genes in the expression matrix cannot be scored: its scores
              are all NaN (not 0.0, which would read as a real score).
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
        raise ValueError(f"gene_names length ({len(gene_names)}) must match columns ({n_genes})")

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
            logger.warning(f"Gene set '{gs_name}' has none of its genes in the expression matrix: scores are NaN")
            scores[gs_name] = [math.nan] * n_cells
            continue

        non_gs = set(gs_indices)
        non_gs_indices = [i for i in all_indices if i not in non_gs]
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
        converted: list[list[float]] = data.tolist()
        return converted

    if hasattr(data, "toarray"):
        dense: list[list[float]] = data.toarray().tolist()
        return dense

    raise TypeError(f"Unsupported expression matrix type: {type(data)}")


class _TestOutcome(NamedTuple):
    """Outcome of one two-sample test.

    ``p_value`` is NaN, and ``status`` names the reason, when the test cannot be computed; ``status`` is
    ``"ok"`` otherwise. A p-value is never invented to stand in for a test that did not run.
    """

    statistic: float
    p_value: float
    status: str


_STATUS_OK = "ok"


def _undefined(status: str) -> _TestOutcome:
    """The outcome of a test that could not be computed: no statistic, no p-value, and the reason."""
    return _TestOutcome(math.nan, math.nan, status)


def _all_finite(values: Sequence[float]) -> bool:
    """Whether every value is a finite number."""
    if HAS_NUMPY:
        return bool(np.isfinite(np.asarray(values, dtype=float)).all())
    return all(math.isfinite(v) for v in values)


def _mean_and_variance(values: Sequence[float]) -> tuple[float, float]:
    """Mean and unbiased sample variance of at least two values.

    Constant values have variance exactly 0.0: a variance computed as a sum of squared deviations from a
    rounded mean can come out as 1e-33 instead, and a t statistic built on that would be astronomically large.
    """
    n = len(values)
    if HAS_NUMPY:
        array = np.asarray(values, dtype=float)
        mean = float(array.mean())
        variance = 0.0 if array.min() == array.max() else float(array.var(ddof=1))
        return mean, variance
    mean = sum(values) / n
    if min(values) == max(values):
        return mean, 0.0
    return mean, sum((v - mean) ** 2 for v in values) / (n - 1)


def _wilcoxon_rank_sum(a: Sequence[float], b: Sequence[float]) -> _TestOutcome:
    """Perform Wilcoxon rank-sum test (Mann-Whitney U test).

    Tie-aware by construction: with scipy, uses
    ``scipy.stats.mannwhitneyu`` (exact for small samples without ties,
    tie-corrected asymptotic normal approximation otherwise), which is
    essential because single-cell counts are massively tied. The pure
    Python fallback applies the tie correction term
    ``sum(t**3 - t) / (12 * n * (n - 1))`` to the normal approximation.

    The test is undefined, and gives NaN with a status instead of a p-value, when a group is empty, when a
    value is not finite, and when every value of both groups is tied (the rank-sum statistic then has zero
    variance).

    Args:
        a: Values from group A.
        b: Values from group B.

    Returns:
        The U statistic of group A, the two-sided p-value and ``"ok"``; or NaN, NaN and the reason.
    """
    n_a = len(a)
    n_b = len(b)
    if n_a == 0 or n_b == 0:
        return _undefined("empty_group: a group has no observations")
    if not (_all_finite(a) and _all_finite(b)):
        return _undefined("non_finite_values: the values contain NaN or infinity")
    if min(min(a), min(b)) == max(max(a), max(b)):
        return _undefined("zero_variance: every value of both groups is tied, so the rank-sum statistic has no spread")

    if HAS_SCIPY:
        try:
            result = scipy_stats.mannwhitneyu(a, b, alternative="two-sided")
        except (ValueError, ArithmeticError) as exc:
            return _undefined(f"error: {type(exc).__name__}: {exc}")
        p_value = float(result.pvalue)
        if math.isnan(p_value):
            return _undefined("numerical_failure: scipy.stats.mannwhitneyu returned NaN")
        return _TestOutcome(float(result.statistic), p_value, _STATUS_OK)

    # Pure Python fallback: Mann-Whitney U with tie-corrected normal
    # approximation

    # Combine and rank
    combined = [(v, 0) for v in a] + [(v, 1) for v in b]
    combined.sort(key=lambda x: x[0])

    # Assign ranks (average ranks for ties) and accumulate the tie term
    ranks: list[float] = [0.0] * len(combined)
    tie_term = 0.0
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        t = j - i
        avg_rank = (i + j + 1) / 2.0  # 1-indexed average
        for k in range(i, j):
            ranks[k] = avg_rank
        if t > 1:
            tie_term += t**3 - t
        i = j

    # Sum of ranks for group A
    rank_sum_a = sum(ranks[i] for i in range(len(combined)) if combined[i][1] == 0)

    # U statistic
    u_a = rank_sum_a - n_a * (n_a + 1) / 2.0

    # Normal approximation with tie correction
    n = n_a + n_b
    mu = n_a * n_b / 2.0
    sigma_sq = (n_a * n_b / 12.0) * ((n + 1) - tie_term / (n * (n - 1)))
    if sigma_sq <= 0:
        return _undefined("zero_variance: the rank-sum statistic has no spread")

    z = (u_a - mu) / math.sqrt(sigma_sq)

    # Two-sided p-value from the standard normal
    return _TestOutcome(u_a, 2.0 * _standard_normal_cdf(-abs(z)), _STATUS_OK)


def _welch_t_test(a: Sequence[float], b: Sequence[float]) -> _TestOutcome:
    """Perform Welch's t-test for unequal variances.

    The statistic is ``(mean_a - mean_b) / sqrt(var_a / n_a + var_b / n_b)`` on the Welch-Satterthwaite
    degrees of freedom, with a two-sided p-value from the exact Student t distribution (scipy's, or the
    in-module incomplete-beta evaluation when scipy is unavailable; never a normal approximation).

    The test is undefined, and gives NaN with a status instead of a p-value, with fewer than two values in
    a group, with non-finite values, and when both groups are constant (the standard error is then zero).

    Args:
        a: Values from group A.
        b: Values from group B.

    Returns:
        The t statistic, the two-sided p-value and ``"ok"``; or NaN, NaN and the reason.
    """
    n_a = len(a)
    n_b = len(b)
    if n_a < 2 or n_b < 2:
        return _undefined("insufficient_observations: Welch's t-test needs at least 2 values in each group")
    if not (_all_finite(a) and _all_finite(b)):
        return _undefined("non_finite_values: the values contain NaN or infinity")

    mean_a, var_a = _mean_and_variance(a)
    mean_b, var_b = _mean_and_variance(b)
    squared_error_a = var_a / n_a
    squared_error_b = var_b / n_b
    squared_error = squared_error_a + squared_error_b
    if squared_error == 0.0:
        return _undefined("zero_variance: both groups are constant, so the t statistic is undefined")

    t_statistic = (mean_a - mean_b) / math.sqrt(squared_error)
    degrees_of_freedom = squared_error**2 / (squared_error_a**2 / (n_a - 1) + squared_error_b**2 / (n_b - 1))
    p_value = _student_t_two_sided_p(t_statistic, degrees_of_freedom)
    if math.isnan(p_value):
        return _undefined("numerical_failure: the Student t tail probability did not converge")
    return _TestOutcome(t_statistic, p_value, _STATUS_OK)


def _student_t_two_sided_p(t_statistic: float, degrees_of_freedom: float) -> float:
    """Two-sided tail probability ``P(|T| >= |t|)`` of Student's t distribution (NaN if it cannot be evaluated)."""
    if HAS_SCIPY:
        return min(1.0, float(2.0 * scipy_stats.t.sf(abs(t_statistic), degrees_of_freedom)))
    # P(|T| >= |t|) = I_x(df / 2, 1 / 2) with x = df / (df + t^2)
    x = degrees_of_freedom / (degrees_of_freedom + t_statistic * t_statistic)
    return _regularized_incomplete_beta(degrees_of_freedom / 2.0, 0.5, x)


def _regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta ``I_x(a, b)`` by Lentz's continued fraction (NaN if it does not converge)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_front = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x)
    front = math.exp(log_front)
    # The continued fraction converges quickly below the mean of the beta distribution; use the symmetry above it
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _beta_continued_fraction(a, b, x) / a
    return 1.0 - front * _beta_continued_fraction(b, a, 1.0 - x) / b


def _beta_continued_fraction(a: float, b: float, x: float, max_iterations: int = 100_000) -> float:
    """Continued fraction of the incomplete beta function (modified Lentz); NaN if it does not converge."""
    tiny = 1e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, max_iterations + 1):
        m2 = 2 * m
        # Even step
        step = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + step * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + step / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        # Odd step
        step = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + step * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + step / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            return h
    return math.nan


def _standard_normal_cdf(x: float) -> float:
    """CDF of the standard normal distribution, ``P(Z <= x)``, from the complementary error function.

    Args:
        x: Value at which to evaluate the CDF.

    Returns:
        ``P(Z <= x)`` for ``Z ~ N(0, 1)``, accurate to machine precision in both tails (a rational erf
        approximation is only good to about 1e-7 absolute, which is useless for p-values below that).
    """
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def _benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Apply Benjamini-Hochberg FDR correction to a list of p-values.

    NaN p-values (tests that were undefined) take no part in the correction: the family size is the
    number of p-values that exist, and a NaN input yields a NaN output at its position.

    Args:
        p_values: Raw p-values (NaN where there is none).

    Returns:
        Adjusted p-values (same length and order as input).
    """
    adjusted = [math.nan] * len(p_values)
    defined = [(index, p_value) for index, p_value in enumerate(p_values) if not math.isnan(p_value)]
    n = len(defined)
    if n == 0:
        return adjusted

    # Sort p-values and track original indices
    indexed = sorted(defined, key=lambda x: x[1])

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


def _adjusted_p_sort_key(result: dict) -> tuple[bool, float]:
    """Sort key that orders results by adjusted p-value and puts the ones without a p-value (NaN) last."""
    adjusted_p = result["adjusted_p"]
    return (math.isnan(adjusted_p), 0.0 if math.isnan(adjusted_p) else adjusted_p)


def _warn_about_undefined_tests(results: list[dict]) -> None:
    """Log one warning that counts the genes whose test was undefined, by reason."""
    reasons: dict[str, int] = {}
    for result in results:
        status = result["test_status"]
        if status != _STATUS_OK:
            reason = status.split(":", 1)[0]
            reasons[reason] = reasons.get(reason, 0) + 1
    if reasons:
        summary = ", ".join(f"{count} {reason}" for reason, count in sorted(reasons.items()))
        logger.warning(
            f"{sum(reasons.values())} of {len(results)} genes have no p-value ({summary}): "
            "p_value and adjusted_p are NaN and they are excluded from the Benjamini-Hochberg correction"
        )


__all__ = [
    "compute_log_fold_change",
    "differential_expression",
    "gene_set_scoring",
    "pseudobulk_de",
    "volcano_data",
]
