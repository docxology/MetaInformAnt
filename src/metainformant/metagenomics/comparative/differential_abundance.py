"""Differential abundance testing for metagenomic data.

Provides statistical methods for identifying taxa whose abundance differs
significantly between groups of samples. Includes compositional data
analysis (CLR transformation), indicator species analysis, effect size
ranking (LEfSe-style), and machine-learning-based biomarker discovery.

Methods:
    - ALDEx2-like: Per-sample Dirichlet Monte Carlo instances of the
      compositional proportions (``n_monte_carlo`` seeded draws), averaged
      in CLR space, followed by Welch's t-test with effect size estimation.
      Accounts for both the compositional nature and the count-sampling
      uncertainty of microbiome data.
    - ANCOM-like: Pairwise log-ratio comparisons to identify taxa that
      are differentially abundant relative to most other taxa. The result is
      the W statistic (the number of rejected pairwise tests); ANCOM defines
      no p-value, so none is reported.
    - Simple DESeq-style: median-of-ratios size factors followed by Welch's
      t-test on the normalised counts (no negative binomial model).
    - Indicator species (IndVal): Identifies taxa strongly associated
      with specific groups based on specificity and fidelity.
    - LEfSe-style effect size ranking: a univariate standardised mean
      difference per taxon with a Kruskal-Wallis test; it is not a
      multivariate Linear Discriminant Analysis.
    - Random forest biomarker discovery: Uses feature importance from
      random forest classification to identify discriminative taxa.

Undefined tests: a test that cannot be computed (a group with fewer than two
samples, no variance in either group, identical values for Kruskal-Wallis)
reports a NaN p-value, never 1.0. NaN p-values are left out of the
Benjamini-Hochberg correction (it ranks and scales over the tests that ran)
and sort after every defined p-value.
"""

from __future__ import annotations

import math
import random
from collections import Counter
from typing import Any, cast

from metainformant.core.utils.logging import get_logger
from metainformant.metagenomics._special import chi2_sf, student_t_two_sided_p

logger = get_logger(__name__)

# Optional scientific dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = cast(Any, None)  # optional dependency; callers guard on HAS_NUMPY

try:
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = cast(Any, None)  # optional dependency; callers guard on HAS_SCIPY

try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import cross_val_score

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    RandomForestClassifier = None
    cross_val_score = None

# Share of the pairwise log-ratio tests a taxon must win to count as detected by the ANCOM-like W statistic
_ANCOM_DETECTED_FRACTION = 0.7


def _as_matrix(counts: Any) -> list[list[Any]]:
    """Return a plain nested list for lists of rows, numpy arrays or pandas frames.

    Truth-value tests such as ``if not counts`` raise for arrays with more than one
    element, so every public entry point normalises its input once and then only
    uses ``len()``.

    Args:
        counts: Matrix-like object (samples x taxa).

    Returns:
        The rows as lists.

    Raises:
        ValueError: If the input is not two-dimensional.
    """
    if hasattr(counts, "to_numpy"):
        counts = counts.to_numpy()
    if hasattr(counts, "tolist"):
        counts = counts.tolist()
    matrix: list[list[Any]] = []
    for row in counts:
        if isinstance(row, (str, bytes)) or not hasattr(row, "__iter__"):
            raise ValueError("counts must be a 2-D matrix (samples x taxa)")
        matrix.append(list(row))
    return matrix


def _as_list(values: Any) -> list[Any]:
    """Return a plain list (Python scalars) for lists, numpy arrays or pandas series."""
    if hasattr(values, "to_numpy"):
        values = values.to_numpy()
    if hasattr(values, "tolist"):
        values = values.tolist()
    return list(values)


def differential_abundance(
    counts: list[list[int]],
    groups: list[int],
    taxa_names: list[str],
    method: str = "aldex2_like",
    n_monte_carlo: int = 128,
    seed: int | None = None,
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
        method: Statistical method. One of "aldex2_like" (Dirichlet Monte
            Carlo CLR + Welch t-test), "ancom_like" (pairwise log-ratio W
            statistic; reports no p-value), or "simple_deseq" (median-of-ratios
            size factors + Welch t-test on the normalised counts).
        n_monte_carlo: Number of Monte Carlo instances for CLR estimation
            (used in aldex2_like method). Must be >= 1.
        seed: Random seed for the Dirichlet Monte Carlo draws (aldex2_like).

    Returns:
        List of dictionaries, sorted by adjusted p-value (NaN last) or, for
        ``ancom_like``, by W statistic (descending). Each contains:
            - taxon: str taxon name.
            - log2fc: float log2 fold change between groups.
            - p_value: float raw two-sided p-value; NaN when the test is
              undefined (a group with fewer than two samples, or no variance in
              either group) and for ``ancom_like``, which has none.
            - adjusted_p: float BH-adjusted p-value (NaN where p_value is NaN;
              NaN p-values are excluded from the correction).
            - effect_size: float standardized effect size (Cohen's d; for
              ``ancom_like`` the W statistic as a fraction of the other taxa).
            - mean_group1: float mean abundance in group 1.
            - mean_group2: float mean abundance in group 2.
        ``ancom_like`` rows add ``w_statistic`` (int, the number of other taxa
        whose log-ratio to this taxon differs between the groups at p < 0.05),
        ``n_ratio_tests`` (int, the pairwise tests that could be computed) and
        ``detected`` (bool, ``w_statistic`` at least 70% of ``n_ratio_tests``).

    Raises:
        ValueError: If inputs are inconsistent or groups does not contain
            exactly 2 unique values.
    """
    valid_methods = ("aldex2_like", "ancom_like", "simple_deseq")
    if method not in valid_methods:
        raise ValueError(f"Invalid method '{method}'. Must be one of {valid_methods}")

    counts = _as_matrix(counts)
    groups = _as_list(groups)
    taxa_names = _as_list(taxa_names)

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
        results = _aldex2_like(counts, g1_indices, g2_indices, taxa_names, n_monte_carlo, seed=seed)
    elif method == "ancom_like":
        results = _ancom_like(counts, g1_indices, g2_indices, taxa_names)
    else:  # simple_deseq
        results = _simple_deseq(counts, g1_indices, g2_indices, taxa_names)

    if method == "ancom_like":
        # ANCOM defines no p-value: rank by the W statistic
        results.sort(
            key=lambda x: (math.isnan(x["effect_size"]), 0.0 if math.isnan(x["effect_size"]) else -x["effect_size"])
        )
        n_sig = sum(1 for r in results if r["detected"])
        logger.info(f"DA analysis complete: {len(results)} taxa, {n_sig} detected (W at least 70% of ratio tests)")
        return results

    # Sort by adjusted p-value; undefined tests (NaN) go last
    results.sort(key=lambda x: (math.isnan(x["adjusted_p"]), 0.0 if math.isnan(x["adjusted_p"]) else x["adjusted_p"]))

    n_undefined = sum(1 for r in results if math.isnan(r["p_value"]))
    n_sig = sum(1 for r in results if r["adjusted_p"] < 0.05)
    logger.info(f"DA analysis complete: {len(results)} taxa tested, {n_sig} significant (FDR < 0.05)")
    if n_undefined:
        logger.warning(f"{n_undefined} of {len(results)} taxa have an undefined test (NaN p-value)")

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
    counts = _as_matrix(counts)
    if len(counts) == 0:
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
            - associated_group: int group with highest indicator value; None when the
              taxon is absent from every sample (all its indicator values are zero).

    Raises:
        ValueError: If inputs are inconsistent.
    """
    counts = _as_matrix(counts)
    groups = _as_list(groups)
    taxa_names = _as_list(taxa_names)

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
        best_group: int | None = None

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
    """Rank taxa by a univariate effect size with a Kruskal-Wallis test (LEfSe-inspired).

    Each taxon is analysed independently on CLR-transformed counts:

    - ``p_value`` is the Kruskal-Wallis test across the groups (NaN when the
      test is undefined, for example when every value is identical).
    - ``lda_score`` is ``log10(1 + effect)``. With two groups ``effect`` is the
      absolute difference of the group means divided by the pooled standard
      deviation (the unstandardised difference when that is zero); with more
      groups it is the largest absolute difference between two group means,
      in CLR units, without standardisation.

    This is not the multivariate linear discriminant analysis of LEfSe: no
    discriminant direction is fitted across taxa, so the scores are not
    comparable with LEfSe's LDA scores.

    Args:
        counts: Count matrix (samples x taxa).
        groups: Group labels for each sample.
        taxa_names: Taxon names.

    Returns:
        Sorted list (by absolute ``lda_score`` descending) of dictionaries:
            - taxon: str taxon name.
            - lda_score: float univariate effect size (log10 scale).
            - p_value: float Kruskal-Wallis p-value (NaN when undefined).
            - direction: str "enriched_in_<group>" naming the group with the
              higher mean CLR abundance; None when the group means are equal.

    Raises:
        ValueError: If inputs are inconsistent.
    """
    counts = _as_matrix(counts)
    groups = _as_list(groups)
    taxa_names = _as_list(taxa_names)

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

            direction = None if mean_diff == 0 else (f"enriched_in_{g1}" if mean_diff > 0 else f"enriched_in_{g0}")
        else:
            # For multiple groups, use max pairwise effect
            max_effect = 0.0
            direction = None
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
        cv_folds: Requested number of stratified cross-validation folds (at
            least 2). It is capped at the size of the smallest group, and the
            number actually used is reported in ``model_summary["cv_folds"]``.
        seed: Random seed for reproducibility.

    Returns:
        Dictionary with keys:
            - selected_taxa: list[str] of the taxa whose importance exceeds the
              mean importance, ranked by importance.
            - importances: dict mapping taxon name to importance score.
            - cv_accuracy: float mean cross-validated accuracy, or None when
              cross-validation could not run. It is never reported as 0.0 for a
              missing estimate.
            - cv_unavailable_reason: str saying why ``cv_accuracy`` is None
              (smallest group has fewer than 2 samples, or scikit-learn is
              not installed); None when cross-validation ran.
            - model_summary: dict with model parameters and metrics
              (``cv_folds`` is None and ``cv_scores`` empty when cross-validation
              did not run).

    Raises:
        ValueError: If inputs are inconsistent, ``groups`` has fewer than two
            distinct labels, or ``cv_folds`` is below 2.
    """
    if method != "random_forest":
        raise ValueError(f"Invalid method '{method}'. Must be 'random_forest'")
    if cv_folds < 2:
        raise ValueError(f"cv_folds must be at least 2, got {cv_folds}")

    counts = _as_matrix(counts)
    groups = _as_list(groups)
    taxa_names = _as_list(taxa_names)

    n_samples = len(counts)
    n_taxa = len(counts[0]) if n_samples > 0 else 0

    if len(groups) != n_samples:
        raise ValueError(f"groups length ({len(groups)}) must match samples ({n_samples})")
    if len(taxa_names) != n_taxa:
        raise ValueError(f"taxa_names length ({len(taxa_names)}) must match taxa ({n_taxa})")
    if len(set(groups)) < 2:
        raise ValueError("groups must contain at least two distinct labels to find discriminative taxa")

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
    seed: int | None = None,
) -> list[dict]:
    """ALDEx2-like DA: Dirichlet Monte Carlo CLR + Welch's t-test.

    For every sample, ``n_monte_carlo`` Dirichlet draws (0.5 pseudocount
    prior per taxon, as in ALDEx2) simulate compositional sampling
    uncertainty. Each draw is CLR-transformed and averaged within the
    sample; Welch's t-test and Cohen's d then operate on these
    Monte-Carlo-averaged CLR values, so the reported p-values reflect both
    compositional structure and MC-estimated sampling noise.

    Args:
        counts: Count matrix.
        g1_indices: Sample indices for group 1.
        g2_indices: Sample indices for group 2.
        taxa_names: Taxon names.
        n_monte_carlo: Monte Carlo instances (must be >= 1).
        seed: Random seed for reproducibility (None = nondeterministic).

    Returns:
        List of DA result dictionaries.
    """
    if n_monte_carlo < 1:
        raise ValueError(f"n_monte_carlo must be >= 1, got {n_monte_carlo}")

    n_taxa = len(counts[0]) if len(counts) > 0 else 0

    if HAS_NUMPY:
        rng = np.random.default_rng(seed)
        pseudocount = 0.5
        test_values: list[list[float]] = []
        for sample_row in counts:
            alpha = np.asarray(sample_row, dtype=float) + pseudocount
            draws = rng.dirichlet(alpha, size=n_monte_carlo)
            log_draws = np.log(draws)
            clr_draws = log_draws - log_draws.mean(axis=1, keepdims=True)
            test_values.append(clr_draws.mean(axis=0).tolist())
    else:
        # numpy unavailable: draw the same Dirichlet instances with the standard library
        test_values = _dirichlet_clr_means(counts, n_monte_carlo, seed)

    results: list[dict] = []
    raw_p_values: list[float] = []

    for taxon_idx in range(n_taxa):
        vals_g1 = [test_values[i][taxon_idx] for i in g1_indices]
        vals_g2 = [test_values[i][taxon_idx] for i in g2_indices]

        # Raw means for fold change
        raw_g1 = [counts[i][taxon_idx] for i in g1_indices]
        raw_g2 = [counts[i][taxon_idx] for i in g2_indices]
        mean_raw_g1 = sum(raw_g1) / len(raw_g1) if raw_g1 else 0.0
        mean_raw_g2 = sum(raw_g2) / len(raw_g2) if raw_g2 else 0.0

        log2fc = _log2fc(mean_raw_g1, mean_raw_g2)

        # Welch's t-test on Monte-Carlo-averaged CLR values
        p_value = _welch_t_test(vals_g1, vals_g2)

        # Effect size (Cohen's d)
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


def _dirichlet_clr_means(counts: list[list[int]], n_monte_carlo: int, seed: int | None) -> list[list[float]]:
    """Monte-Carlo-averaged CLR of Dirichlet(count + 0.5) draws, using only the standard library.

    A Dirichlet draw is a vector of independent Gamma(alpha_j) variates divided by their sum. The
    sum cancels in the centred log-ratio (``log g_j - mean_k log g_k``), so the gamma variates are
    used directly. This is the same Monte Carlo as the NumPy path, with a different random stream.

    Args:
        counts: Count matrix (samples x taxa).
        n_monte_carlo: Number of Dirichlet instances per sample (>= 1).
        seed: Seed of the random generator (None draws fresh entropy).

    Returns:
        One row per sample: the mean CLR over the draws.
    """
    rng = random.Random(seed)
    means: list[list[float]] = []
    for row in counts:
        alphas = [float(v) + 0.5 for v in row]
        total = [0.0] * len(alphas)
        for _ in range(n_monte_carlo):
            # the floor only guards log(0) for a gamma variate that underflowed
            logs = [math.log(max(rng.gammavariate(alpha, 1.0), 5e-324)) for alpha in alphas]
            centre = sum(logs) / len(logs)
            for j, value in enumerate(logs):
                total[j] += value - centre
        means.append([value / n_monte_carlo for value in total])
    return means


def _ancom_like(
    counts: list[list[int]],
    g1_indices: list[int],
    g2_indices: list[int],
    taxa_names: list[str],
) -> list[dict]:
    """ANCOM-like DA: pairwise log-ratio comparisons and the W statistic.

    For each taxon, computes log-ratios against all other taxa and tests
    each for a difference between the groups (Welch's t-test, p < 0.05).
    ``w_statistic`` counts the rejections. ANCOM defines no p-value, so
    ``p_value`` and ``adjusted_p`` are NaN; the evidence is the W statistic
    (``effect_size`` is W as a fraction of the tests that could be computed)
    and ``detected`` marks taxa with W of at least 70% of those tests. A
    pairwise test that is undefined (no variance in either group, or a group
    with a single sample) is neither a rejection nor counted as a test.

    Args:
        counts: Count matrix.
        g1_indices: Group 1 sample indices.
        g2_indices: Group 2 sample indices.
        taxa_names: Taxon names.

    Returns:
        List of DA result dictionaries (see ``differential_abundance``).
    """
    n_taxa = len(counts[0]) if len(counts) > 0 else 0
    pseudocount = 0.5

    results: list[dict] = []

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
            if math.isnan(p):
                continue
            if p < 0.05:
                n_reject += 1
            n_tests += 1

        # W as a fraction of the tests that ran; undefined when none could run
        w_fraction = n_reject / n_tests if n_tests > 0 else math.nan

        raw_g1 = [counts[s][taxon_i] for s in g1_indices]
        raw_g2 = [counts[s][taxon_i] for s in g2_indices]
        mean_g1 = sum(raw_g1) / len(raw_g1) if raw_g1 else 0.0
        mean_g2 = sum(raw_g2) / len(raw_g2) if raw_g2 else 0.0

        results.append(
            {
                "taxon": taxa_names[taxon_i],
                "log2fc": _log2fc(mean_g1, mean_g2),
                "p_value": math.nan,
                "adjusted_p": math.nan,
                "effect_size": w_fraction,
                "mean_group1": mean_g1,
                "mean_group2": mean_g2,
                "w_statistic": n_reject,
                "n_ratio_tests": n_tests,
                "detected": n_tests > 0 and n_reject >= _ANCOM_DETECTED_FRACTION * n_tests,
            }
        )

    return results


def _simple_deseq(
    counts: list[list[int]],
    g1_indices: list[int],
    g2_indices: list[int],
    taxa_names: list[str],
) -> list[dict]:
    """Simplified DESeq2-style analysis using size factor normalization.

    Normalizes counts using median-of-ratios size factors (the geometric mean of a
    taxon is taken over its non-zero samples) and compares the normalized counts
    of the two groups with Welch's t-test. No negative binomial model or
    dispersion estimate is fitted, so this is not DESeq2's Wald test.

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

    cv_scores: list[float] = []
    cv_accuracy: float | None = None
    cv_unavailable_reason: str | None = None
    actual_folds: int | None = None
    # Stratified k-fold needs every class to have at least ``n_splits`` members, so the
    # smallest class bounds the folds; the number of classes is irrelevant. A class with a
    # single sample cannot be split between a training and a test fold at all.
    smallest_class = min(Counter(groups).values())
    if smallest_class < 2:
        cv_unavailable_reason = (
            f"cross-validation needs at least 2 samples in every group; the smallest group has {smallest_class}"
        )
        logger.warning(f"Cross-validation skipped: {cv_unavailable_reason}")
    else:
        actual_folds = min(cv_folds, smallest_class)
        cv_scores = [float(s) for s in cross_val_score(rf, X, y, cv=actual_folds, scoring="accuracy")]
        cv_accuracy = float(np.mean(cv_scores))
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
        "cv_scores": cv_scores,
        "n_features": len(taxa_names),
        "n_selected": len(selected),
    }

    accuracy_text = "unavailable" if cv_accuracy is None else f"{cv_accuracy:.3f}"
    logger.info(f"Biomarker discovery complete: {len(selected)} selected, CV accuracy={accuracy_text}")

    return {
        "selected_taxa": selected,
        "importances": importances,
        "cv_accuracy": cv_accuracy,
        "cv_unavailable_reason": cv_unavailable_reason,
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
    size as a proxy for feature importance. No classifier is fitted, so
    ``cv_accuracy`` is None and ``cv_unavailable_reason`` says so.

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
    n_taxa = len(clr_data[0]) if len(clr_data) > 0 else 0

    importances: dict[str, float] = {}

    for taxon_idx in range(n_taxa):
        group_vals: dict[int, list[float]] = {}
        for i, grp in enumerate(groups):
            if grp not in group_vals:
                group_vals[grp] = []
            group_vals[grp].append(clr_data[i][taxon_idx])

        # Use absolute Cohen's d as importance proxy. A taxon whose effect size is undefined (no
        # variance, or too few samples) carries no evidence of discriminating the groups, so it
        # ranks with importance 0 rather than disturbing the ordering with NaN.
        pairwise_effects: list[float] = []
        for gi in range(len(unique_groups)):
            for gj in range(gi + 1, len(unique_groups)):
                effect = _cohens_d(group_vals.get(unique_groups[gi], []), group_vals.get(unique_groups[gj], []))
                if not math.isnan(effect):
                    pairwise_effects.append(abs(effect))

        importances[taxa_names[taxon_idx]] = max(pairwise_effects, default=0.0)

    sorted_taxa = sorted(importances.keys(), key=lambda t: importances[t], reverse=True)
    mean_imp = sum(importances.values()) / len(importances) if importances else 0.0
    selected = [t for t in sorted_taxa if importances[t] > mean_imp]

    return {
        "selected_taxa": selected,
        "importances": importances,
        "cv_accuracy": None,
        "cv_unavailable_reason": (
            "scikit-learn is not installed: the effect-size fallback fits no classifier, "
            "so there is nothing to cross-validate"
        ),
        "model_summary": {
            "method": "effect_size_ranking",
            "note": "scikit-learn not available; used Cohen's d as proxy",
            "cv_folds": None,
            "cv_scores": [],
            "n_features": len(taxa_names),
            "n_selected": len(selected),
        },
    }


# ---------------------------------------------------------------------------
# Private helpers: Statistics
# ---------------------------------------------------------------------------


def _welch_t_test(a: list[float], b: list[float]) -> float:
    """Welch's t-test for unequal variances (two-sided).

    The statistic and the Welch-Satterthwaite degrees of freedom are computed here and the
    p-value is the exact Student's t tail (a regularized incomplete beta function), so the
    result is identical with and without SciPy and matches ``scipy.stats.ttest_ind(...,
    equal_var=False)``.

    Args:
        a: Group A values.
        b: Group B values.

    Returns:
        Two-sided p-value, or NaN when the test is undefined: a group has fewer than two
        values, neither group has any variance (the statistic is then 0/0 or x/0), or a
        value is not finite.
    """
    n_a = len(a)
    n_b = len(b)
    if n_a < 2 or n_b < 2:
        return math.nan

    # The statistic does not depend on the unit of the data; working with unit-scale values keeps
    # the squared variances below from overflowing or underflowing.
    scale = max(abs(x) for x in (*a, *b))
    if not (math.isfinite(scale) and scale > 0.0):
        return math.nan
    a = [x / scale for x in a]
    b = [x / scale for x in b]

    mean_a = sum(a) / n_a
    mean_b = sum(b) / n_b
    var_of_mean_a = sum((x - mean_a) ** 2 for x in a) / (n_a - 1) / n_a
    var_of_mean_b = sum((x - mean_b) ** 2 for x in b) / (n_b - 1) / n_b

    variance_of_difference = var_of_mean_a + var_of_mean_b
    if not variance_of_difference > 0.0:
        return math.nan

    welch_denominator = var_of_mean_a**2 / (n_a - 1) + var_of_mean_b**2 / (n_b - 1)
    if not welch_denominator > 0.0:
        return math.nan  # variances so small that their squares underflow: the degrees of freedom are 0 / 0

    t_stat = (mean_a - mean_b) / math.sqrt(variance_of_difference)
    return student_t_two_sided_p(t_stat, variance_of_difference**2 / welch_denominator)


def _kruskal_wallis_undefined(groups_data: list[list[float]]) -> bool:
    """True when the Kruskal-Wallis statistic does not exist.

    That is the case with fewer than two groups, an empty group, or every value identical.
    """
    if len(groups_data) < 2 or any(len(g) == 0 for g in groups_data):
        return True
    return len({value for g in groups_data for value in g}) < 2


def _kruskal_wallis_pure(groups_data: list[list[float]]) -> float:
    """Kruskal-Wallis H test for any number of groups, without SciPy.

    Pools all values, assigns mid-ranks to ties, and computes
    ``H = 12 / (N (N + 1)) * sum(R_i ** 2 / n_i) - 3 (N + 1)`` divided by the tie
    correction ``1 - sum(t ** 3 - t) / (N ** 3 - N)``. The p-value is the chi-square
    survival function with ``k - 1`` degrees of freedom, matching
    ``scipy.stats.kruskal``.

    Args:
        groups_data: List of value lists (one per group).

    Returns:
        P-value; NaN when the test is undefined (fewer than two non-empty groups, or
        all values identical).
    """
    groups = [list(g) for g in groups_data]
    if _kruskal_wallis_undefined(groups):
        return math.nan

    n_total = sum(len(g) for g in groups)
    pooled = sorted((value, idx) for idx, g in enumerate(groups) for value in g)

    rank_sums = [0.0] * len(groups)
    tie_term = 0.0
    i = 0
    while i < n_total:
        j = i
        while j < n_total and pooled[j][0] == pooled[i][0]:
            j += 1
        mid_rank = (i + j + 1) / 2.0  # ranks are 1-based: positions i..j-1 hold ranks i+1..j
        for k in range(i, j):
            rank_sums[pooled[k][1]] += mid_rank
        tie_size = j - i
        tie_term += tie_size**3 - tie_size
        i = j

    tie_correction = 1.0 - tie_term / (n_total**3 - n_total) if n_total > 1 else 0.0
    if tie_correction <= 0.0:
        return math.nan  # all values identical: H is undefined

    h = 12.0 / (n_total * (n_total + 1)) * sum(r**2 / len(g) for r, g in zip(rank_sums, groups)) - 3.0 * (n_total + 1)
    h = max(h / tie_correction, 0.0)
    return min(chi2_sf(h, len(groups) - 1), 1.0)


def _kruskal_wallis_two_group(a: list[float], b: list[float]) -> float:
    """Kruskal-Wallis test for two groups (equivalent to Mann-Whitney).

    Args:
        a: Group A values.
        b: Group B values.

    Returns:
        P-value; NaN when the test is undefined.
    """
    return _kruskal_wallis([a, b])


def _kruskal_wallis(groups_data: list[list[float]]) -> float:
    """Kruskal-Wallis test for any number of groups.

    Uses ``scipy.stats.kruskal`` when SciPy is installed and the exact SciPy-free
    implementation otherwise; the two agree.

    Args:
        groups_data: List of value lists (one per group).

    Returns:
        P-value; NaN when the test is undefined.
    """
    if _kruskal_wallis_undefined(groups_data):
        return math.nan
    if HAS_SCIPY:
        return float(scipy_stats.kruskal(*groups_data).pvalue)
    return _kruskal_wallis_pure(groups_data)


def _cohens_d(a: list[float], b: list[float]) -> float:
    """Compute Cohen's d effect size.

    Args:
        a: Group A values.
        b: Group B values.

    Returns:
        Cohen's d (positive = A > B); NaN when it is undefined (a group has fewer than
        two values, or the pooled standard deviation is zero).
    """
    n_a = len(a)
    n_b = len(b)
    if n_a < 2 or n_b < 2:
        return math.nan

    mean_a = sum(a) / n_a
    mean_b = sum(b) / n_b
    var_a = sum((x - mean_a) ** 2 for x in a) / (n_a - 1)
    var_b = sum((x - mean_b) ** 2 for x in b) / (n_b - 1)

    pooled_std = math.sqrt(((n_a - 1) * var_a + (n_b - 1) * var_b) / (n_a + n_b - 2))
    if pooled_std == 0:
        return math.nan

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

    NaN p-values (tests that could not be computed) are left out of the correction:
    they stay NaN, and the ranks and the number of tests ``m`` count only the defined
    p-values, so a failed test cannot inflate or dilute the family.

    Args:
        p_values: Raw p-values (NaN allowed).

    Returns:
        Adjusted p-values (NaN where the input was NaN).
    """
    adjusted = [math.nan] * len(p_values)
    defined = sorted((i for i, p in enumerate(p_values) if not math.isnan(p)), key=lambda i: p_values[i])
    n = len(defined)
    min_so_far = 1.0

    for rank_idx in range(n - 1, -1, -1):
        orig_idx = defined[rank_idx]
        adj = min(p_values[orig_idx] * n / (rank_idx + 1), min_so_far, 1.0)
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
