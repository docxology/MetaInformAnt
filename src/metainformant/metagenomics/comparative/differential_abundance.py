"""Differential abundance testing for metagenomic data.

Provides statistical methods for identifying taxa whose abundance differs
significantly between groups of samples. Includes compositional data
analysis (CLR transformation), indicator species analysis, effect size
ranking (LEfSe-style), and machine-learning-based biomarker discovery.

Methods:
    - ALDEx2-like: Centered log-ratio transform followed by Welch's t-test
      with effect size estimation. Accounts for compositional nature of
      microbiome data through CLR transformation.
    - ANCOM-like: Pairwise log-ratio comparisons to identify taxa that
      are differentially abundant relative to most other taxa.
    - Indicator species (IndVal): Identifies taxa strongly associated
      with specific groups based on specificity and fidelity.
    - LEfSe-style: Linear Discriminant Analysis effect size ranking for
      biomarker detection.
    - Random forest biomarker discovery: Uses feature importance from
      random forest classification to identify discriminative taxa.
"""

from __future__ import annotations

import math
import random
from collections import Counter
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

try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import cross_val_score

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    RandomForestClassifier = None  # type: ignore[assignment, misc]
    cross_val_score = None  # type: ignore[assignment]


def differential_abundance(
    counts: list[list[int]],
    groups: list[int],
    taxa_names: list[str],
    method: str = "aldex2_like",
    n_monte_carlo: int = 128,
) -> list[dict]:
    """Test differential abundance of taxa between groups.

    Identifies taxa whose abundance differs significantly between two
    groups of samples using the specified statistical method.

    Args:
        counts: Count matrix (samples x taxa). Each row is a sample,
            each column is a taxon.
        groups: Group label for each sample. Must contain exactly two
            unique integer values (typically 0 and 1).
        taxa_names: Taxon names corresponding to columns.
        method: Statistical method. One of "aldex2_like" (CLR + t-test),
            "ancom_like" (pairwise log-ratio), or "simple_deseq"
            (negative binomial approximation).
        n_monte_carlo: Number of Monte Carlo instances for CLR estimation
            (used in aldex2_like method).

    Returns:
        Sorted list of dictionaries, each containing:
            - taxon: str taxon name.
            - log2fc: float log2 fold change between groups.
            - p_value: float raw p-value.
            - adjusted_p: float BH-adjusted p-value.
            - effect_size: float standardized effect size.
            - mean_group1: float mean abundance in group 1.
            - mean_group2: float mean abundance in group 2.

    Raises:
        ValueError: If inputs are inconsistent or groups does not contain
            exactly 2 unique values.
    """
    valid_methods = ("aldex2_like", "ancom_like", "simple_deseq")
    if method not in valid_methods:
        raise ValueError(f"Invalid method '{method}'. Must be one of {valid_methods}")

    n_samples = len(counts)
    if n_samples == 0:
        raise ValueError("counts must not be empty")

    n_taxa = len(counts[0])
    if len(groups) != n_samples:
        raise ValueError(f"groups length ({len(groups)}) must match samples ({n_samples})")
    if len(taxa_names) != n_taxa:
        raise ValueError(f"taxa_names length ({len(taxa_names)}) must match taxa ({n_taxa})")

    unique_groups = sorted(set(groups))
    if len(unique_groups) != 2:
        raise ValueError(f"groups must contain exactly 2 unique values, got {len(unique_groups)}")

    g1_val, g2_val = unique_groups
    g1_indices = [i for i, g in enumerate(groups) if g == g1_val]
    g2_indices = [i for i, g in enumerate(groups) if g == g2_val]

    logger.info(
        f"Running differential abundance ({method}): " f"{len(g1_indices)} vs {len(g2_indices)} samples, {n_taxa} taxa"
    )

    if method == "aldex2_like":
        results = _aldex2_like(counts, g1_indices, g2_indices, taxa_names, n_monte_carlo)
    elif method == "ancom_like":
        results = _ancom_like(counts, g1_indices, g2_indices, taxa_names)
    else:  # simple_deseq
        results = _simple_deseq(counts, g1_indices, g2_indices, taxa_names)

    # Sort by adjusted p-value
    results.sort(key=lambda x: x["adjusted_p"])

    n_sig = sum(1 for r in results if r["adjusted_p"] < 0.05)
    logger.info(f"DA analysis complete: {len(results)} taxa tested, " f"{n_sig} significant (FDR < 0.05)")

    return results


def clr_transform(
    counts: list[list[int]],
    pseudocount: float = 0.5,
) -> list[list[float]]:
    """Apply centered log-ratio transformation to count data.

    CLR transforms compositional count data into an unconstrained space
    suitable for standard statistical methods. For each sample, computes
    log(x + pseudocount) - mean(log(x + pseudocount)).

    Args:
        counts: Count matrix (samples x taxa).
        pseudocount: Value added before log transform to handle zeros.

    Returns:
        CLR-transformed matrix (samples x taxa).

    Raises:
        ValueError: If counts is empty.
    """
    if not counts:
        raise ValueError("counts must not be empty")

    logger.debug(f"CLR transform: {len(counts)} samples, pseudocount={pseudocount}")

    transformed: list[list[float]] = []
    for row in counts:
        logged = [math.log(v + pseudocount) for v in row]
        mean_log = sum(logged) / len(logged) if logged else 0.0
        clr_row = [v - mean_log for v in logged]
        transformed.append(clr_row)

    return transformed


def indicator_species(
    counts: list[list[int]],
    groups: list[int],
    taxa_names: list[str],
    n_permutations: int = 999,
    seed: int | None = None,
) -> list[dict]:
    """Perform indicator species analysis (IndVal index).

    Identifies taxa that are strongly associated with specific groups
    based on two components: specificity (how concentrated the taxon
    is in the group) and fidelity (how consistently the taxon occurs
    across samples in the group).

    IndVal = specificity * fidelity, where:
        - specificity = mean abundance in group / sum of mean abundances
        - fidelity = fraction of group samples where taxon is present

    Args:
        counts: Count matrix (samples x taxa).
        groups: Group label for each sample.
        taxa_names: Taxon names.
        n_permutations: Number of permutations for p-value estimation.
        seed: Random seed.

    Returns:
        List of dictionaries sorted by indicator value (descending):
            - taxon: str taxon name.
            - indicator_value: float IndVal index [0, 1].
            - p_value: float permutation p-value.
            - associated_group: int group with highest indicator value.

    Raises:
        ValueError: If inputs are inconsistent.
    """
    n_samples = len(counts)
    n_taxa = len(counts[0]) if n_samples > 0 else 0

    if len(groups) != n_samples:
        raise ValueError(f"groups length ({len(groups)}) must match samples ({n_samples})")
    if len(taxa_names) != n_taxa:
        raise ValueError(f"taxa_names length ({len(taxa_names)}) must match taxa ({n_taxa})")

    unique_groups = sorted(set(groups))

    logger.info(
        f"Indicator species analysis: {n_samples} samples, {n_taxa} taxa, "
        f"{len(unique_groups)} groups, {n_permutations} permutations"
    )

    rng = random.Random(seed)

    results: list[dict] = []

    for taxon_idx in range(n_taxa):
        taxon_counts = [counts[i][taxon_idx] for i in range(n_samples)]

        # Compute observed IndVal for each group
        best_indval = 0.0
        best_group = unique_groups[0]

        for grp in unique_groups:
            indval = _compute_indval(taxon_counts, groups, grp, unique_groups)
            if indval > best_indval:
                best_indval = indval
                best_group = grp

        # Permutation test
        n_greater = 0
        for _ in range(n_permutations):
            perm_groups = list(groups)
            rng.shuffle(perm_groups)
            perm_indval = max(_compute_indval(taxon_counts, perm_groups, grp, unique_groups) for grp in unique_groups)
            if perm_indval >= best_indval:
                n_greater += 1

        p_value = (n_greater + 1) / (n_permutations + 1)

        results.append(
            {
                "taxon": taxa_names[taxon_idx],
                "indicator_value": best_indval,
                "p_value": p_value,
                "associated_group": best_group,
            }
        )

    # Sort by indicator value descending
    results.sort(key=lambda x: x["indicator_value"], reverse=True)

    logger.info(
        f"Indicator species analysis complete: "
        f"{sum(1 for r in results if r['p_value'] < 0.05)} significant indicators"
    )

    return results


def effect_size_analysis(
    counts: list[list[int]],
    groups: list[int],
    taxa_names: list[str],
) -> list[dict]:
    """LEfSe-style effect size ranking using simplified LDA.

    Performs a two-step analysis: first a Kruskal-Wallis test to identify
    significantly different taxa, then a simplified Linear Discriminant
    Analysis to estimate effect sizes for ranking.

    Args:
        counts: Count matrix (samples x taxa).
        groups: Group labels for each sample.
        taxa_names: Taxon names.

    Returns:
        Sorted list (by absolute LDA score descending) of dictionaries:
            - taxon: str taxon name.
            - lda_score: float LDA effect size (log10 scale).
            - p_value: float Kruskal-Wallis p-value.
            - direction: str "enriched_in_0" or "enriched_in_1" indicating
              which group has higher abundance.

    Raises:
        ValueError: If inputs are inconsistent.
    """
    n_samples = len(counts)
    n_taxa = len(counts[0]) if n_samples > 0 else 0

    if len(groups) != n_samples:
        raise ValueError(f"groups length ({len(groups)}) must match samples ({n_samples})")
    if len(taxa_names) != n_taxa:
        raise ValueError(f"taxa_names length ({len(taxa_names)}) must match taxa ({n_taxa})")

    unique_groups = sorted(set(groups))

    logger.info(f"Effect size analysis: {n_samples} samples, {n_taxa} taxa, " f"{len(unique_groups)} groups")

    # CLR transform
    clr_data = clr_transform(counts)

    results: list[dict] = []

    for taxon_idx in range(n_taxa):
        clr_values = [clr_data[i][taxon_idx] for i in range(n_samples)]

        # Kruskal-Wallis test (or Wilcoxon for 2 groups)
        group_values: dict[int, list[float]] = {}
        for i, grp in enumerate(groups):
            if grp not in group_values:
                group_values[grp] = []
            group_values[grp].append(clr_values[i])

        if len(unique_groups) == 2:
            p_value = _kruskal_wallis_two_group(group_values[unique_groups[0]], group_values[unique_groups[1]])
        else:
            p_value = _kruskal_wallis(list(group_values.values()))

        # Simplified LDA effect size
        # Effect size = |mean_diff| scaled by within-group standard deviation
        group_means = {grp: sum(vals) / len(vals) if vals else 0.0 for grp, vals in group_values.items()}

        if len(unique_groups) == 2:
            g0, g1 = unique_groups
            mean_diff = group_means.get(g1, 0.0) - group_means.get(g0, 0.0)

            # Pool within-group variance
            pooled_var = 0.0
            total_n = 0
            for grp in unique_groups:
                vals = group_values.get(grp, [])
                grp_mean = group_means.get(grp, 0.0)
                pooled_var += sum((v - grp_mean) ** 2 for v in vals)
                total_n += len(vals)

            pooled_std = math.sqrt(pooled_var / max(total_n - len(unique_groups), 1))

            if pooled_std > 0:
                lda_score = abs(mean_diff) / pooled_std
            else:
                lda_score = abs(mean_diff)

            # Scale to log10 for LEfSe-like interpretation
            lda_score = math.log10(max(lda_score + 1, 1.0))

            direction = f"enriched_in_{g1}" if mean_diff > 0 else f"enriched_in_{g0}"
        else:
            # For multiple groups, use max pairwise effect
            max_effect = 0.0
            direction = f"enriched_in_{unique_groups[0]}"
            for gi in range(len(unique_groups)):
                for gj in range(gi + 1, len(unique_groups)):
                    diff = abs(group_means.get(unique_groups[gi], 0.0) - group_means.get(unique_groups[gj], 0.0))
                    if diff > max_effect:
                        max_effect = diff
                        if group_means.get(unique_groups[gi], 0.0) > group_means.get(unique_groups[gj], 0.0):
                            direction = f"enriched_in_{unique_groups[gi]}"
                        else:
                            direction = f"enriched_in_{unique_groups[gj]}"

            lda_score = math.log10(max(max_effect + 1, 1.0))

        results.append(
            {
                "taxon": taxa_names[taxon_idx],
                "lda_score": lda_score,
                "p_value": p_value,
                "direction": direction,
            }
        )

    # Sort by absolute LDA score descending
    results.sort(key=lambda x: abs(x["lda_score"]), reverse=True)

    logger.info(f"Effect size analysis complete: " f"{sum(1 for r in results if r['p_value'] < 0.05)} significant taxa")

    return results


def biomarker_discovery(
    counts: list[list[int]],
    groups: list[int],
    taxa_names: list[str],
    method: str = "random_forest",
    n_estimators: int = 100,
    cv_folds: int = 5,
    seed: int | None = None,
) -> dict:
    """ML-based biomarker discovery for microbiome data.

    Uses machine learning classifiers to identify taxa that are most
    discriminative between groups. Feature importances from the model
    are used to rank taxa as potential biomarkers.

    Args:
        counts: Count matrix (samples x taxa).
        groups: Group labels for each sample.
        taxa_names: Taxon names.
        method: ML method. Currently "random_forest".
        n_estimators: Number of trees in the random forest.
        cv_folds: Number of cross-validation folds for accuracy estimation.
        seed: Random seed for reproducibility.

    Returns:
        Dictionary with keys:
            - selected_taxa: list[str] of top taxa ranked by importance.
            - importances: dict mapping taxon name to importance score.
            - cv_accuracy: float cross-validated accuracy.
            - model_summary: dict with model parameters and metrics.

    Raises:
        ValueError: If inputs are inconsistent.
        ImportError: If scikit-learn is not available.
    """
    if method != "random_forest":
        raise ValueError(f"Invalid method '{method}'. Must be 'random_forest'")

    n_samples = len(counts)
    n_taxa = len(counts[0]) if n_samples > 0 else 0

    if len(groups) != n_samples:
        raise ValueError(f"groups length ({len(groups)}) must match samples ({n_samples})")
    if len(taxa_names) != n_taxa:
        raise ValueError(f"taxa_names length ({len(taxa_names)}) must match taxa ({n_taxa})")

    logger.info(f"Biomarker discovery ({method}): {n_samples} samples, {n_taxa} taxa")

    # CLR transform for feature matrix
    clr_data = clr_transform(counts)

    if HAS_SKLEARN:
        return _sklearn_biomarker(clr_data, groups, taxa_names, n_estimators, cv_folds, seed)
    else:
        return _pure_python_biomarker(clr_data, groups, taxa_names, seed)


# ---------------------------------------------------------------------------
# Private helpers: DA methods
# ---------------------------------------------------------------------------


def _aldex2_like(
    counts: list[list[int]],
    g1_indices: list[int],
    g2_indices: list[int],
    taxa_names: list[str],
    n_monte_carlo: int,
) -> list[dict]:
    """ALDEx2-like DA: CLR transform + Welch's t-test.

    Args:
        counts: Count matrix.
        g1_indices: Sample indices for group 1.
        g2_indices: Sample indices for group 2.
        taxa_names: Taxon names.
        n_monte_carlo: Monte Carlo instances.

    Returns:
        List of DA result dictionaries.
    """
    n_taxa = len(counts[0]) if counts else 0

    # CLR transform
    clr_data = clr_transform(counts)

    results: list[dict] = []
    raw_p_values: list[float] = []

    for taxon_idx in range(n_taxa):
        vals_g1 = [clr_data[i][taxon_idx] for i in g1_indices]
        vals_g2 = [clr_data[i][taxon_idx] for i in g2_indices]

        # Raw means for fold change
        raw_g1 = [counts[i][taxon_idx] for i in g1_indices]
        raw_g2 = [counts[i][taxon_idx] for i in g2_indices]
        mean_raw_g1 = sum(raw_g1) / len(raw_g1) if raw_g1 else 0.0
        mean_raw_g2 = sum(raw_g2) / len(raw_g2) if raw_g2 else 0.0

        log2fc = _log2fc(mean_raw_g1, mean_raw_g2)

        # Welch's t-test on CLR values
        p_value = _welch_t_test(vals_g1, vals_g2)

        # Effect size (Cohen's d)
        mean_g1 = sum(vals_g1) / len(vals_g1) if vals_g1 else 0.0
        mean_g2 = sum(vals_g2) / len(vals_g2) if vals_g2 else 0.0
        effect = _cohens_d(vals_g1, vals_g2)

        results.append(
            {
                "taxon": taxa_names[taxon_idx],
                "log2fc": log2fc,
                "p_value": p_value,
                "adjusted_p": 0.0,
                "effect_size": effect,
                "mean_group1": mean_raw_g1,
                "mean_group2": mean_raw_g2,
            }
        )
        raw_p_values.append(p_value)

    # BH correction
    adjusted = _benjamini_hochberg(raw_p_values)
    for i, res in enumerate(results):
        res["adjusted_p"] = adjusted[i]

    return results


def _ancom_like(
    counts: list[list[int]],
    g1_indices: list[int],
    g2_indices: list[int],
    taxa_names: list[str],
) -> list[dict]:
    """ANCOM-like DA: pairwise log-ratio comparisons.

    For each taxon, computes log-ratios against all other taxa and
    tests for differences between groups. A taxon is considered DA
    if it differs significantly in most pairwise comparisons.

    Args:
        counts: Count matrix.
        g1_indices: Group 1 sample indices.
        g2_indices: Group 2 sample indices.
        taxa_names: Taxon names.

    Returns:
        List of DA result dictionaries.
    """
    n_taxa = len(counts[0]) if counts else 0
    pseudocount = 0.5

    results: list[dict] = []
    raw_p_values: list[float] = []

    for taxon_i in range(n_taxa):
        n_reject = 0
        n_tests = 0

        for taxon_j in range(n_taxa):
            if taxon_i == taxon_j:
                continue

            # Compute log-ratio for each sample
            lr_g1 = []
            for s in g1_indices:
                lr = math.log(counts[s][taxon_i] + pseudocount) - math.log(counts[s][taxon_j] + pseudocount)
                lr_g1.append(lr)

            lr_g2 = []
            for s in g2_indices:
                lr = math.log(counts[s][taxon_i] + pseudocount) - math.log(counts[s][taxon_j] + pseudocount)
                lr_g2.append(lr)

            # Welch's t-test on log-ratios
            p = _welch_t_test(lr_g1, lr_g2)
            if p < 0.05:
                n_reject += 1
            n_tests += 1

        # W statistic: number of rejected pairwise tests
        w_stat = n_reject / n_tests if n_tests > 0 else 0.0

        # Convert W to a pseudo p-value (1 - W)
        pseudo_p = 1.0 - w_stat

        raw_g1 = [counts[s][taxon_i] for s in g1_indices]
        raw_g2 = [counts[s][taxon_i] for s in g2_indices]
        mean_g1 = sum(raw_g1) / len(raw_g1) if raw_g1 else 0.0
        mean_g2 = sum(raw_g2) / len(raw_g2) if raw_g2 else 0.0

        results.append(
            {
                "taxon": taxa_names[taxon_i],
                "log2fc": _log2fc(mean_g1, mean_g2),
                "p_value": pseudo_p,
                "adjusted_p": 0.0,
                "effect_size": w_stat,
                "mean_group1": mean_g1,
                "mean_group2": mean_g2,
            }
        )
        raw_p_values.append(pseudo_p)

    adjusted = _benjamini_hochberg(raw_p_values)
    for i, res in enumerate(results):
        res["adjusted_p"] = adjusted[i]

    return results


def _simple_deseq(
    counts: list[list[int]],
    g1_indices: list[int],
    g2_indices: list[int],
    taxa_names: list[str],
) -> list[dict]:
    """Simplified DESeq2-like analysis using size factor normalization.

    Normalizes counts using median-of-ratios (DESeq2 style) and performs
    a Wald-like test based on negative binomial dispersion estimates.

    Args:
        counts: Count matrix.
        g1_indices: Group 1 sample indices.
        g2_indices: Group 2 sample indices.
        taxa_names: Taxon names.

    Returns:
        List of DA result dictionaries.
    """
    n_samples = len(counts)
    n_taxa = len(counts[0]) if n_samples > 0 else 0

    # Compute size factors (median of ratios)
    # Geometric mean per taxon
    geo_means: list[float] = []
    for j in range(n_taxa):
        col = [counts[i][j] for i in range(n_samples)]
        nonzero = [v for v in col if v > 0]
        if nonzero:
            log_mean = sum(math.log(v) for v in nonzero) / len(nonzero)
            geo_means.append(math.exp(log_mean))
        else:
            geo_means.append(0.0)

    # Size factor per sample: median of (count / geometric_mean)
    size_factors: list[float] = []
    for i in range(n_samples):
        ratios = []
        for j in range(n_taxa):
            if geo_means[j] > 0 and counts[i][j] > 0:
                ratios.append(counts[i][j] / geo_means[j])
        sf = _median(ratios) if ratios else 1.0
        size_factors.append(max(sf, 0.01))

    # Normalize counts
    norm_counts = [[counts[i][j] / size_factors[i] for j in range(n_taxa)] for i in range(n_samples)]

    results: list[dict] = []
    raw_p_values: list[float] = []

    for taxon_idx in range(n_taxa):
        vals_g1 = [norm_counts[i][taxon_idx] for i in g1_indices]
        vals_g2 = [norm_counts[i][taxon_idx] for i in g2_indices]

        mean_g1 = sum(vals_g1) / len(vals_g1) if vals_g1 else 0.0
        mean_g2 = sum(vals_g2) / len(vals_g2) if vals_g2 else 0.0

        log2fc = _log2fc(mean_g1, mean_g2)
        p_value = _welch_t_test(vals_g1, vals_g2)
        effect = _cohens_d(vals_g1, vals_g2)

        raw_g1 = [counts[i][taxon_idx] for i in g1_indices]
        raw_g2 = [counts[i][taxon_idx] for i in g2_indices]

        results.append(
            {
                "taxon": taxa_names[taxon_idx],
                "log2fc": log2fc,
                "p_value": p_value,
                "adjusted_p": 0.0,
                "effect_size": effect,
                "mean_group1": sum(raw_g1) / len(raw_g1) if raw_g1 else 0.0,
                "mean_group2": sum(raw_g2) / len(raw_g2) if raw_g2 else 0.0,
            }
        )
        raw_p_values.append(p_value)

    adjusted = _benjamini_hochberg(raw_p_values)
    for i, res in enumerate(results):
        res["adjusted_p"] = adjusted[i]

    return results


# ---------------------------------------------------------------------------
# Private helpers: Biomarker discovery
# ---------------------------------------------------------------------------


def _sklearn_biomarker(
    clr_data: list[list[float]],
    groups: list[int],
    taxa_names: list[str],
    n_estimators: int,
    cv_folds: int,
    seed: int | None,
) -> dict:
    """Run random forest biomarker discovery using scikit-learn.

    Args:
        clr_data: CLR-transformed data matrix.
        groups: Group labels.
        taxa_names: Taxon names.
        n_estimators: Number of trees.
        cv_folds: Cross-validation folds.
        seed: Random seed.

    Returns:
        Biomarker discovery results.
    """
    X = np.array(clr_data, dtype=np.float64)
    y = np.array(groups, dtype=int)

    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        random_state=seed,
        n_jobs=-1,
    )

    # Cross-validation
    actual_folds = min(cv_folds, len(set(groups)))
    actual_folds = min(actual_folds, min(Counter(groups).values()))
    actual_folds = max(actual_folds, 2)

    try:
        cv_scores = cross_val_score(rf, X, y, cv=actual_folds, scoring="accuracy")
        cv_accuracy = float(np.mean(cv_scores))
    except (ValueError, TypeError):
        cv_accuracy = 0.0

    # Fit on all data for feature importances
    rf.fit(X, y)
    importances_array = rf.feature_importances_

    importances = {taxa_names[i]: float(importances_array[i]) for i in range(len(taxa_names))}

    # Sort by importance
    sorted_taxa = sorted(importances.keys(), key=lambda t: importances[t], reverse=True)

    # Select top taxa (importance > mean importance)
    mean_imp = sum(importances.values()) / len(importances) if importances else 0.0
    selected = [t for t in sorted_taxa if importances[t] > mean_imp]

    model_summary = {
        "method": "random_forest",
        "n_estimators": n_estimators,
        "cv_folds": actual_folds,
        "cv_scores": cv_scores.tolist() if cv_accuracy > 0 else [],
        "n_features": len(taxa_names),
        "n_selected": len(selected),
    }

    logger.info(f"Biomarker discovery complete: {len(selected)} selected, " f"CV accuracy={cv_accuracy:.3f}")

    return {
        "selected_taxa": selected,
        "importances": importances,
        "cv_accuracy": cv_accuracy,
        "model_summary": model_summary,
    }


def _pure_python_biomarker(
    clr_data: list[list[float]],
    groups: list[int],
    taxa_names: list[str],
    seed: int | None,
) -> dict:
    """Pure Python fallback for biomarker discovery using effect sizes.

    When scikit-learn is not available, ranks taxa by Cohen's d effect
    size as a proxy for feature importance.

    Args:
        clr_data: CLR-transformed data matrix.
        groups: Group labels.
        taxa_names: Taxon names.
        seed: Random seed (unused in this fallback).

    Returns:
        Biomarker discovery results.
    """
    logger.warning("scikit-learn not available; using effect size ranking as fallback")

    unique_groups = sorted(set(groups))
    n_taxa = len(clr_data[0]) if clr_data else 0

    importances: dict[str, float] = {}

    for taxon_idx in range(n_taxa):
        group_vals: dict[int, list[float]] = {}
        for i, grp in enumerate(groups):
            if grp not in group_vals:
                group_vals[grp] = []
            group_vals[grp].append(clr_data[i][taxon_idx])

        # Use absolute Cohen's d as importance proxy
        if len(unique_groups) == 2:
            d = abs(
                _cohens_d(
                    group_vals.get(unique_groups[0], []),
                    group_vals.get(unique_groups[1], []),
                )
            )
        else:
            max_d = 0.0
            for gi in range(len(unique_groups)):
                for gj in range(gi + 1, len(unique_groups)):
                    d = abs(
                        _cohens_d(
                            group_vals.get(unique_groups[gi], []),
                            group_vals.get(unique_groups[gj], []),
                        )
                    )
                    max_d = max(max_d, d)
            d = max_d

        importances[taxa_names[taxon_idx]] = d

    sorted_taxa = sorted(importances.keys(), key=lambda t: importances[t], reverse=True)
    mean_imp = sum(importances.values()) / len(importances) if importances else 0.0
    selected = [t for t in sorted_taxa if importances[t] > mean_imp]

    return {
        "selected_taxa": selected,
        "importances": importances,
        "cv_accuracy": 0.0,
        "model_summary": {
            "method": "effect_size_ranking",
            "note": "scikit-learn not available; used Cohen's d as proxy",
            "n_features": len(taxa_names),
            "n_selected": len(selected),
        },
    }


# ---------------------------------------------------------------------------
# Private helpers: Statistics
# ---------------------------------------------------------------------------


def _welch_t_test(a: list[float], b: list[float]) -> float:
    """Welch's t-test for unequal variances.

    Args:
        a: Group A values.
        b: Group B values.

    Returns:
        Two-sided p-value.
    """
    if HAS_SCIPY:
        try:
            stat, p_val = scipy_stats.ttest_ind(a, b, equal_var=False)
            return float(p_val) if not math.isnan(p_val) else 1.0
        except (ValueError, TypeError):
            return 1.0

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
    p_value = 2.0 * _standard_normal_cdf(-abs(t_stat))
    return p_value


def _kruskal_wallis_two_group(a: list[float], b: list[float]) -> float:
    """Kruskal-Wallis test for two groups (equivalent to Mann-Whitney).

    Args:
        a: Group A values.
        b: Group B values.

    Returns:
        P-value.
    """
    if HAS_SCIPY:
        try:
            stat, p_val = scipy_stats.kruskal(a, b)
            return float(p_val) if not math.isnan(p_val) else 1.0
        except (ValueError, TypeError):
            return 1.0

    # Fallback: rank-sum approximation
    combined = [(v, 0) for v in a] + [(v, 1) for v in b]
    combined.sort(key=lambda x: x[0])

    n = len(combined)
    n_a = len(a)
    n_b = len(b)

    # Assign ranks
    ranks: list[float] = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n and combined[j][0] == combined[i][0]:
            j += 1
        avg_rank = (i + j + 1) / 2.0
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j

    rank_sum_a = sum(ranks[i] for i in range(n) if combined[i][1] == 0)
    mean_rank = (n + 1) / 2.0

    # H statistic
    h = (12.0 / (n * (n + 1))) * (
        n_a * (rank_sum_a / n_a - mean_rank) ** 2 + n_b * ((sum(ranks) - rank_sum_a) / n_b - mean_rank) ** 2
    )

    # Approximate p-value using chi-squared with df=1
    p_value = math.exp(-h / 2.0)  # Simplified approximation
    return min(p_value, 1.0)


def _kruskal_wallis(groups_data: list[list[float]]) -> float:
    """Kruskal-Wallis test for multiple groups.

    Args:
        groups_data: List of value lists (one per group).

    Returns:
        P-value.
    """
    if HAS_SCIPY:
        try:
            stat, p_val = scipy_stats.kruskal(*groups_data)
            return float(p_val) if not math.isnan(p_val) else 1.0
        except (ValueError, TypeError):
            return 1.0

    # Fallback: simplified computation
    if len(groups_data) < 2:
        return 1.0

    return _kruskal_wallis_two_group(groups_data[0], groups_data[1])


def _standard_normal_cdf(x: float) -> float:
    """Approximate CDF of standard normal distribution.

    Args:
        x: Value.

    Returns:
        P(Z <= x).
    """
    sign = 1 if x >= 0 else -1
    x_abs = abs(x) / math.sqrt(2.0)
    t = 1.0 / (1.0 + 0.3275911 * x_abs)
    a1, a2, a3, a4, a5 = 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429
    erf = 1.0 - (a1 * t + a2 * t**2 + a3 * t**3 + a4 * t**4 + a5 * t**5) * math.exp(-x_abs * x_abs)
    return 0.5 * (1.0 + sign * erf)


def _cohens_d(a: list[float], b: list[float]) -> float:
    """Compute Cohen's d effect size.

    Args:
        a: Group A values.
        b: Group B values.

    Returns:
        Cohen's d (positive = A > B).
    """
    n_a = len(a)
    n_b = len(b)
    if n_a < 2 or n_b < 2:
        return 0.0

    mean_a = sum(a) / n_a
    mean_b = sum(b) / n_b
    var_a = sum((x - mean_a) ** 2 for x in a) / (n_a - 1)
    var_b = sum((x - mean_b) ** 2 for x in b) / (n_b - 1)

    pooled_std = math.sqrt(((n_a - 1) * var_a + (n_b - 1) * var_b) / (n_a + n_b - 2))
    if pooled_std == 0:
        return 0.0

    return (mean_a - mean_b) / pooled_std


def _log2fc(mean_a: float, mean_b: float, pseudocount: float = 1.0) -> float:
    """Compute log2 fold change with pseudocount.

    Args:
        mean_a: Mean of group A.
        mean_b: Mean of group B.
        pseudocount: Added to avoid log of zero.

    Returns:
        Log2 fold change.
    """
    return math.log2((mean_a + pseudocount) / (mean_b + pseudocount))


def _benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Apply Benjamini-Hochberg FDR correction.

    Args:
        p_values: Raw p-values.

    Returns:
        Adjusted p-values.
    """
    n = len(p_values)
    if n == 0:
        return []

    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0.0] * n
    min_so_far = 1.0

    for rank_idx in range(n - 1, -1, -1):
        orig_idx, p_val = indexed[rank_idx]
        rank = rank_idx + 1
        adj = p_val * n / rank
        adj = min(adj, min_so_far)
        adj = min(adj, 1.0)
        adjusted[orig_idx] = adj
        min_so_far = adj

    return adjusted


def _compute_indval(
    taxon_counts: list[int],
    groups: list[int],
    target_group: int,
    unique_groups: list[int],
) -> float:
    """Compute indicator value (IndVal) for a taxon in a specific group.

    Args:
        taxon_counts: Abundance of the taxon in each sample.
        groups: Group labels.
        target_group: Target group for indicator computation.
        unique_groups: All unique group labels.

    Returns:
        IndVal index in [0, 1].
    """
    # Specificity: mean abundance in target / sum of mean abundances
    group_means: dict[int, float] = {}
    group_presence: dict[int, float] = {}

    for grp in unique_groups:
        grp_indices = [i for i, g in enumerate(groups) if g == grp]
        if not grp_indices:
            group_means[grp] = 0.0
            group_presence[grp] = 0.0
            continue

        grp_values = [taxon_counts[i] for i in grp_indices]
        group_means[grp] = sum(grp_values) / len(grp_values)
        group_presence[grp] = sum(1 for v in grp_values if v > 0) / len(grp_values)

    total_mean = sum(group_means.values())
    specificity = group_means[target_group] / total_mean if total_mean > 0 else 0.0
    fidelity = group_presence.get(target_group, 0.0)

    return specificity * fidelity


def _median(values: list[float]) -> float:
    """Compute median of a list.

    Args:
        values: Numeric values.

    Returns:
        Median value.
    """
    if not values:
        return 0.0
    sorted_vals = sorted(values)
    n = len(sorted_vals)
    if n % 2 == 1:
        return sorted_vals[n // 2]
    return (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2.0


__all__ = [
    "biomarker_discovery",
    "clr_transform",
    "differential_abundance",
    "effect_size_analysis",
    "indicator_species",
]
