"""Community diversity metrics for metagenomic data.

Provides alpha diversity metrics (within-sample diversity), beta diversity
distance matrices (between-sample dissimilarity), rarefaction analysis,
statistical testing of community differences (PERMANOVA), and ordination
methods for visualization (PCoA, NMDS).

Alpha diversity metrics:
    - Shannon: Entropy-based diversity accounting for richness and evenness.
    - Simpson: Probability that two randomly chosen individuals are different.
    - Inverse Simpson: Reciprocal of Simpson; effective number of species.
    - Chao1: Richness estimator accounting for unseen species.
    - ACE: Abundance-based Coverage Estimator for richness.
    - Observed: Simple count of species with non-zero abundance.
    - Fisher alpha: Fisher's log-series alpha diversity parameter.
    - Pielou evenness: Shannon diversity normalized by log of richness.

Beta diversity metrics:
    - Bray-Curtis: Quantitative dissimilarity based on abundance differences.
    - Jaccard: Presence/absence dissimilarity (shared vs. total species).
    - Aitchison: Euclidean distance on CLR-transformed compositions.
"""

from __future__ import annotations

import math
import random
from collections import Counter
from typing import Any, cast

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional scientific dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = cast(Any, None)  # optional dependency; callers guard on HAS_NUMPY


# Estimators defined on integer counts (singletons, doubletons, log-series)
_COUNT_METRICS = ("chao1", "ace", "fisher_alpha")


def _as_list(values: Any) -> list[Any]:
    """Return a plain list (Python scalars) for lists, numpy arrays or pandas series.

    Truth-value tests such as ``if not values`` raise for arrays with more than one
    element, so entry points normalise once and then only use ``len()``.
    """
    if hasattr(values, "to_numpy"):
        values = values.to_numpy()
    if hasattr(values, "tolist"):
        values = values.tolist()
    return list(values)


def _require_integer_counts(values: list[float], metric: str) -> None:
    """Raise ValueError unless every value is an integer count.

    Rounding proportions to integers (0.7 -> 1, 0.2 -> 0, ...) produces meaningless
    singleton/doubleton statistics, so count-based estimators refuse them.

    Args:
        values: Positive abundances.
        metric: Metric name for the error message.

    Raises:
        ValueError: If any value is not within 1e-9 of an integer.
    """
    for value in values:
        if abs(value - round(value)) > 1e-9:
            raise ValueError(
                f"Metric '{metric}' requires integer counts but got the non-integer abundance {value!r}. "
                "Relative abundances (proportions) cannot be used; supply raw or rarefied counts."
            )


def alpha_diversity(
    abundances: list[float] | list[int],
    metric: str = "shannon",
) -> dict:
    """Compute alpha diversity for a single sample.

    Calculates within-sample diversity using the specified metric. Input
    abundances can be raw counts or relative abundances.

    Args:
        abundances: Abundance values for each taxon/OTU in the sample (a list or a
            numpy array). The count-based estimators (chao1, ace, fisher_alpha)
            require integer-valued counts (floats such as 3.0 are fine); proportions
            are rejected. For entropy-based metrics (shannon, simpson,
            invsimpson, pielou_evenness) both counts and proportions are accepted.
        metric: Diversity metric. One of "shannon", "simpson",
            "invsimpson", "chao1", "ace", "observed", "fisher_alpha",
            "pielou_evenness".

    Returns:
        Dictionary with keys:
            - value: float computed diversity value (``inf`` for ``fisher_alpha`` and
              ``ace`` when the estimator has no finite value: every species a singleton,
              or every rare species a singleton).
            - metric: str metric name used.
            - n_species: int number of species with non-zero abundance.
            - total_count: float sum of all abundances.

    Raises:
        ValueError: If metric is not recognized, abundances is empty, or a count-based
            estimator (chao1, ace, fisher_alpha) is given non-integer abundances such
            as relative abundances.
    """
    valid_metrics = (
        "shannon",
        "simpson",
        "invsimpson",
        "chao1",
        "ace",
        "observed",
        "fisher_alpha",
        "pielou_evenness",
    )
    if metric not in valid_metrics:
        raise ValueError(f"Invalid metric '{metric}'. Must be one of {valid_metrics}")

    abundances = _as_list(abundances)
    if len(abundances) == 0:
        raise ValueError("abundances must not be empty")

    # Filter to non-zero values
    nonzero = [a for a in abundances if a > 0]
    n_species = len(nonzero)
    total = sum(nonzero)

    if metric in _COUNT_METRICS:
        _require_integer_counts(nonzero, metric)

    if total == 0:
        return {"value": 0.0, "metric": metric, "n_species": 0, "total_count": 0.0}

    # Convert to proportions for entropy-based metrics
    proportions = [a / total for a in nonzero]

    if metric == "shannon":
        value = _shannon_entropy(proportions)
    elif metric == "simpson":
        value = _simpson_index(proportions)
    elif metric == "invsimpson":
        si = _simpson_index(proportions)
        value = 1.0 / si if si > 0 else float("inf")
    elif metric == "chao1":
        counts = [int(round(a)) for a in nonzero]
        value = _chao1(counts)
    elif metric == "ace":
        counts = [int(round(a)) for a in nonzero]
        value = _ace(counts)
    elif metric == "observed":
        value = float(n_species)
    elif metric == "fisher_alpha":
        counts = [int(round(a)) for a in nonzero]
        value = _fisher_alpha(counts)
    elif metric == "pielou_evenness":
        if n_species <= 1:
            value = 1.0
        else:
            h = _shannon_entropy(proportions)
            value = h / math.log(n_species)
    else:
        value = 0.0

    logger.debug(f"Alpha diversity ({metric}): {value:.4f} " f"({n_species} species, total={total:.1f})")

    return {
        "value": value,
        "metric": metric,
        "n_species": n_species,
        "total_count": float(total),
    }


def beta_diversity(
    samples: list[list[float]],
    metric: str = "bray_curtis",
) -> dict:
    """Compute beta diversity distance matrix between samples.

    Calculates pairwise dissimilarity between all pairs of samples
    using the specified distance metric.

    Args:
        samples: List of abundance vectors (one per sample). All vectors
            must have the same length (number of taxa).
        metric: Distance metric. One of "bray_curtis", "jaccard",
            "aitchison".

    Returns:
        Dictionary with keys:
            - distance_matrix: list[list[float]] symmetric distance matrix.
            - metric: str metric name.
            - n_samples: int number of samples.

    Raises:
        ValueError: If fewer than 2 samples, inconsistent lengths, or
            invalid metric.
    """
    valid_metrics = ("bray_curtis", "jaccard", "aitchison")
    if metric not in valid_metrics:
        raise ValueError(f"Invalid metric '{metric}'. Must be one of {valid_metrics}")

    if len(samples) < 2:
        raise ValueError("At least 2 samples required for beta diversity")

    n_taxa = len(samples[0])
    for i, s in enumerate(samples):
        if len(s) != n_taxa:
            raise ValueError(
                f"Sample {i} has {len(s)} taxa, expected {n_taxa}. " f"All samples must have the same number of taxa."
            )

    n_samples = len(samples)
    logger.info(f"Computing beta diversity ({metric}): {n_samples} samples, {n_taxa} taxa")

    # Pre-transform for Aitchison
    transformed = samples
    if metric == "aitchison":
        transformed = _clr_transform_samples(samples)

    # Compute pairwise distances
    dist_matrix: list[list[float]] = [[0.0] * n_samples for _ in range(n_samples)]

    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            if metric == "bray_curtis":
                d = _bray_curtis(samples[i], samples[j])
            elif metric == "jaccard":
                d = _jaccard_distance(samples[i], samples[j])
            else:  # aitchison
                d = _euclidean(transformed[i], transformed[j])

            dist_matrix[i][j] = d
            dist_matrix[j][i] = d

    logger.info(f"Beta diversity matrix computed: {n_samples}x{n_samples}")

    return {
        "distance_matrix": dist_matrix,
        "metric": metric,
        "n_samples": n_samples,
    }


def rarefaction_curve(
    abundances: list[int],
    depths: list[int] | None = None,
    n_iterations: int = 10,
    seed: int | None = None,
) -> dict:
    """Compute rarefaction curve for a single sample.

    Subsamples the community at increasing sequencing depths and computes
    the mean number of observed species at each depth. Used to assess
    whether sequencing depth is sufficient to capture community diversity.

    Args:
        abundances: Integer counts per taxon in the sample.
        depths: Sampling depths to evaluate. If None, automatically
            generates 20 evenly spaced depths from 1 to total count.
        n_iterations: Number of random subsampling iterations per depth
            for computing mean and standard deviation.
        seed: Random seed for reproducibility.

    Returns:
        Dictionary with keys:
            - depths: list[int] sampling depths evaluated.
            - mean_species: list[float] mean observed species at each depth.
            - std_species: list[float] standard deviation at each depth.
            - is_saturated: bool whether the curve has plateaued (last three
              points within 5% of each other).

    Raises:
        ValueError: If abundances is empty or contains negative values.
    """
    abundances = _as_list(abundances)
    if len(abundances) == 0:
        raise ValueError("abundances must not be empty")

    if any(a < 0 for a in abundances):
        raise ValueError("abundances must not contain negative values")

    total = sum(abundances)
    if total == 0:
        raise ValueError("Total abundance must be > 0")

    # Build pool of individual organisms
    pool: list[int] = []
    for taxon_idx, count in enumerate(abundances):
        pool.extend([taxon_idx] * count)

    if depths is None:
        max_depth = total
        n_points = 20
        depths = [max(1, int(max_depth * i / n_points)) for i in range(1, n_points + 1)]
        # Ensure no duplicates and sorted
        depths = sorted(set(depths))

    # Filter depths that exceed total
    depths = [d for d in depths if d <= total]

    logger.info(
        f"Computing rarefaction curve: total={total}, " f"{len(depths)} depth points, {n_iterations} iterations"
    )

    rng = random.Random(seed)
    mean_species: list[float] = []
    std_species: list[float] = []

    for depth in depths:
        observed_counts: list[int] = []
        for _ in range(n_iterations):
            subsample = rng.sample(pool, depth)
            n_observed = len(set(subsample))
            observed_counts.append(n_observed)

        mean_obs = sum(observed_counts) / len(observed_counts)
        mean_species.append(mean_obs)

        if len(observed_counts) > 1:
            variance = sum((x - mean_obs) ** 2 for x in observed_counts) / (len(observed_counts) - 1)
            std_species.append(math.sqrt(variance))
        else:
            std_species.append(0.0)

    # Check saturation: last 3 points within 5% of each other
    is_saturated = False
    if len(mean_species) >= 3:
        last_three = mean_species[-3:]
        max_val = max(last_three)
        min_val = min(last_three)
        if max_val > 0 and (max_val - min_val) / max_val < 0.05:
            is_saturated = True

    logger.info(
        f"Rarefaction complete: saturated={is_saturated}, " f"max species observed={mean_species[-1]:.1f}"
        if mean_species
        else "Rarefaction complete"
    )

    return {
        "depths": depths,
        "mean_species": mean_species,
        "std_species": std_species,
        "is_saturated": is_saturated,
    }


def rarefy(
    abundances: list[int],
    depth: int,
    seed: int | None = None,
) -> list[int]:
    """Subsample a community to a given sequencing depth.

    Randomly samples individuals without replacement from the community
    and returns the rarefied count vector.

    Args:
        abundances: Integer counts per taxon.
        depth: Target sequencing depth (total reads after rarefaction).
        seed: Random seed for reproducibility.

    Returns:
        List of rarefied counts (same length as input, summing to depth).

    Raises:
        ValueError: If depth exceeds total abundance or is negative.
    """
    total = sum(abundances)
    if depth < 0:
        raise ValueError(f"depth must be non-negative, got {depth}")
    if depth > total:
        raise ValueError(
            f"depth ({depth}) exceeds total abundance ({total}). " f"Cannot rarefy to a depth greater than the total."
        )

    if depth == 0:
        return [0] * len(abundances)

    # Build pool and subsample
    pool: list[int] = []
    for taxon_idx, count in enumerate(abundances):
        pool.extend([taxon_idx] * count)

    rng = random.Random(seed)
    subsample = rng.sample(pool, depth)

    # Count occurrences
    counts = Counter(subsample)
    result = [counts.get(i, 0) for i in range(len(abundances))]
    return result


def permanova(
    distance_matrix: list[list[float]],
    groups: list[str],
    n_permutations: int = 999,
    seed: int | None = None,
) -> dict:
    """PERMANOVA test for group differences in community composition.

    Performs a Permutational Multivariate Analysis of Variance to test
    whether community composition differs significantly between groups.
    The test compares the observed pseudo-F statistic to a null distribution
    generated by permuting group labels.

    Args:
        distance_matrix: Square symmetric distance matrix between samples.
        groups: Group label for each sample.
        n_permutations: Number of permutations for p-value computation.
        seed: Random seed for reproducibility.

    Returns:
        Dictionary with keys:
            - pseudo_f: float observed pseudo-F statistic (``inf`` when the groups are
              perfectly separated, i.e. every within-group distance is zero; NaN when all
              distances are zero).
            - p_value: float permutation p-value (NaN when all distances are zero).
            - r_squared: float fraction of variance explained by groups (NaN when all
              distances are zero).
            - n_permutations: int number of permutations used.

    Raises:
        ValueError: If dimensions are inconsistent, the matrix is not a valid distance
            matrix, there are fewer than 2 groups, or there are no more samples than groups
            (no residual degrees of freedom).
    """
    n = len(distance_matrix)
    if len(groups) != n:
        raise ValueError(f"groups length ({len(groups)}) must match distance matrix " f"dimension ({n})")

    unique_groups = sorted(set(groups))
    if len(unique_groups) < 2:
        raise ValueError("At least 2 groups required for PERMANOVA")
    if n <= len(unique_groups):
        raise ValueError(
            f"PERMANOVA needs more samples ({n}) than groups ({len(unique_groups)}): "
            "with none left over the residual degrees of freedom are zero"
        )
    distance_matrix = _validate_distance_matrix(distance_matrix)

    logger.info(f"Running PERMANOVA: {n} samples, {len(unique_groups)} groups, " f"{n_permutations} permutations")

    # Compute observed pseudo-F
    observed_f = _compute_pseudo_f(distance_matrix, groups, unique_groups)
    if math.isnan(observed_f):
        logger.warning("PERMANOVA is undefined: every distance is zero")
        return {
            "pseudo_f": math.nan,
            "p_value": math.nan,
            "r_squared": math.nan,
            "n_permutations": n_permutations,
        }

    # Permutation test. The comparison allows for rounding: a permutation that reproduces the
    # observed partition can differ from it in the last bits.
    rng = random.Random(seed)
    n_greater = 0
    threshold = observed_f - 1e-12 * max(1.0, abs(observed_f)) if math.isfinite(observed_f) else observed_f

    for _ in range(n_permutations):
        perm_groups = list(groups)
        rng.shuffle(perm_groups)
        perm_f = _compute_pseudo_f(distance_matrix, perm_groups, unique_groups)
        if perm_f >= threshold:
            n_greater += 1

    p_value = (n_greater + 1) / (n_permutations + 1)

    # R-squared: SS_between / SS_total
    ss_total = _compute_ss_total(distance_matrix)
    ss_within = _compute_ss_within(distance_matrix, groups, unique_groups)
    ss_between = ss_total - ss_within
    r_squared = ss_between / ss_total

    logger.info(f"PERMANOVA: pseudo-F={observed_f:.4f}, p={p_value:.4f}, " f"R^2={r_squared:.4f}")

    return {
        "pseudo_f": observed_f,
        "p_value": p_value,
        "r_squared": r_squared,
        "n_permutations": n_permutations,
    }


def ordination(
    distance_matrix: list[list[float]],
    method: str = "pcoa",
    n_components: int = 2,
    seed: int | None = 42,
    n_init: int = 5,
) -> dict:
    """Perform ordination on a distance matrix.

    Reduces the distance matrix to a low-dimensional representation
    for visualization.

    PCoA (classical scaling) eigendecomposes the double-centred squared distance matrix. A
    distance matrix that is not Euclidean (Bray-Curtis, Jaccard and most ecological
    dissimilarities) gives that matrix negative eigenvalues; they have no real coordinates
    and are never hidden. ``all_eigenvalues`` lists the whole spectrum,
    ``n_negative_eigenvalues`` and ``negative_eigenvalue_fraction`` quantify the
    non-Euclidean part, ``variance_explained`` is relative to the sum of the *positive*
    eigenvalues (so it never exceeds 1), and axes whose eigenvalue is not positive get
    coordinates of zero.

    NMDS is Kruskal's non-metric multidimensional scaling: the configuration is fitted so
    that its inter-point distances are a monotone (isotonic) function of the
    dissimilarities, minimising Kruskal's stress-1 by SMACOF majorisation. Ties in the
    dissimilarities are not forced onto equal distances (Kruskal's primary approach).
    Because stress has local minima, the best of ``n_init`` seeded random starts and one
    classical-scaling start is returned. The configuration is centred, rotated to its
    principal axes (sign fixed) and scaled so that its root-mean-square distance equals the
    root-mean-square dissimilarity. NMDS axes have no eigenvalues and no share of variance,
    so none is reported; judge the fit by ``stress`` (below 0.05 excellent, below 0.1 good,
    above 0.2 poor, after Clarke 1993).

    Args:
        distance_matrix: Square symmetric distance matrix with a zero diagonal.
        method: Ordination method. "pcoa" (principal coordinates
            analysis) or "nmds" (non-metric multidimensional scaling).
        n_components: Number of output dimensions.
        seed: Seed for the random starting configurations of NMDS, so repeated calls
            give identical results (default 42; None draws fresh entropy). PCoA is
            deterministic and ignores it.
        n_init: Number of random starting configurations for NMDS, tried in addition to
            the classical-scaling start (default 5; 0 uses the classical start only).
            PCoA ignores it.

    Returns:
        PCoA: dictionary with keys
            - coordinates: list[list[float]] of sample coordinates
              (n_samples x n_components; fewer columns if n_components > n_samples - 1).
            - eigenvalues: list[float] of the leading ``n_components`` eigenvalues.
            - variance_explained: list[float] eigenvalue / sum of positive eigenvalues
              for each retained axis (0.0 for a non-positive eigenvalue).
            - all_eigenvalues: list[float] the complete spectrum, descending, negative
              eigenvalues included.
            - n_negative_eigenvalues: int eigenvalues below zero beyond rounding noise.
            - negative_eigenvalue_fraction: float share of the total absolute eigenvalue
              mass that is negative (0.0 for a Euclidean distance matrix).

        NMDS: dictionary with keys
            - coordinates: list[list[float]] of sample coordinates.
            - stress: float Kruskal stress-1 of the returned configuration (0 is a perfect
              rank-order fit).
            - stress_by_start: list[float] final stress of every start, the classical
              start first.
            - n_iter: int iterations used by the returned start.
            - converged: bool whether the returned start met the stress tolerance before
              the iteration limit.

    Raises:
        ValueError: If method is not recognized, the matrix is not square, symmetric,
            finite, non-negative and zero on its diagonal, ``n_components`` or ``n_init``
            is out of range, or NMDS is given fewer than 3 samples or all-zero
            dissimilarities.
    """
    valid_methods = ("pcoa", "nmds")
    if method not in valid_methods:
        raise ValueError(f"Invalid method '{method}'. Must be one of {valid_methods}")
    if n_components < 1:
        raise ValueError(f"n_components must be at least 1, got {n_components}")
    if n_init < 0:
        raise ValueError(f"n_init must be non-negative, got {n_init}")

    matrix = _validate_distance_matrix(distance_matrix)
    n = len(matrix)

    logger.info(f"Running ordination ({method}): {n} samples, {n_components} components")

    if method == "pcoa":
        return _pcoa(matrix, n_components)
    return _nmds(matrix, n_components, seed=seed, n_init=n_init)


# ---------------------------------------------------------------------------
# Private helpers: Alpha diversity
# ---------------------------------------------------------------------------


def _shannon_entropy(proportions: list[float]) -> float:
    """Compute Shannon entropy H = -sum(p * ln(p)).

    Args:
        proportions: Relative abundance proportions (must sum to ~1.0).

    Returns:
        Shannon entropy (natural log).
    """
    h = 0.0
    for p in proportions:
        if p > 0:
            h -= p * math.log(p)
    return h


def _simpson_index(proportions: list[float]) -> float:
    """Compute Simpson's diversity index D = sum(p^2).

    Note: This is the dominance form. Simpson's diversity = 1 - D.

    Args:
        proportions: Relative abundance proportions.

    Returns:
        Simpson's index (dominance).
    """
    return sum(p * p for p in proportions)


def _chao1(counts: list[int]) -> float:
    """Compute Chao1 richness estimator.

    Chao1 = S_obs + f1^2 / (2 * f2), where f1 = singletons, f2 = doubletons.

    Args:
        counts: Non-zero integer counts per species.

    Returns:
        Chao1 estimated richness.
    """
    s_obs = len(counts)
    f1 = sum(1 for c in counts if c == 1)  # singletons
    f2 = sum(1 for c in counts if c == 2)  # doubletons

    if f2 == 0:
        # Bias-corrected form when f2 = 0
        return s_obs + f1 * (f1 - 1) / 2.0 if f1 > 1 else float(s_obs)

    return s_obs + f1 * f1 / (2.0 * f2)


def _ace(counts: list[int]) -> float:
    """Compute Abundance-based Coverage Estimator (ACE).

    Args:
        counts: Non-zero integer counts per species.

    Returns:
        ACE estimated richness. When every rare species (abundance of at most 10) is a
        singleton the sample coverage is zero, the estimate diverges and ``math.inf`` is
        returned (as for Fisher's alpha with no finite solution), not the observed richness.
    """
    threshold = 10  # Abundance threshold for rare/abundant
    rare = [c for c in counts if c <= threshold]
    abundant = [c for c in counts if c > threshold]

    s_rare = len(rare)
    s_abundant = len(abundant)
    n_rare = sum(rare)

    if n_rare == 0 or s_rare == 0:
        return float(len(counts))

    f1 = sum(1 for c in rare if c == 1)
    c_ace = 1.0 - f1 / n_rare

    if c_ace <= 0:
        return math.inf

    # Coefficient of variation squared
    sum_fi_i = sum(c * (c - 1) for c in rare)
    gamma_sq = (
        max(
            (s_rare * sum_fi_i) / (c_ace * n_rare * (n_rare - 1)) - 1.0,
            0.0,
        )
        if n_rare > 1
        else 0.0
    )

    ace = s_abundant + s_rare / c_ace + f1 * gamma_sq / c_ace
    return ace


def _fisher_alpha(counts: list[int]) -> float:
    """Compute Fisher's log-series alpha.

    Solves ``S = alpha * ln(1 + N / alpha)`` for ``alpha``, where ``S`` is the number
    of species and ``N`` the number of individuals. The left side minus ``S`` is
    strictly increasing in ``alpha`` (from ``-S`` at 0 to ``N - S`` at infinity), so
    there is a unique root exactly when ``N > S``; it is found by bisection on a
    bracket, which cannot diverge. When every species is a singleton (``N == S``)
    no finite alpha exists and ``math.inf`` is returned.

    Args:
        counts: Integer counts per species (zeros are ignored).

    Returns:
        Fisher's alpha parameter (0.0 for an empty community, ``inf`` if all
        species are singletons).
    """
    positive = [c for c in counts if c > 0]
    s = len(positive)
    n = sum(positive)

    if s == 0 or n == 0:
        return 0.0
    if n <= s:
        return math.inf

    def excess(alpha: float) -> float:
        return alpha * math.log1p(n / alpha) - s

    lo = 1e-12  # excess(lo) is about -s < 0
    hi = float(s)
    for _ in range(2000):
        if excess(hi) >= 0.0:
            break
        lo = hi
        hi *= 2.0
    else:
        return math.inf

    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if excess(mid) < 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 1e-14 * hi:
            break

    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Private helpers: Beta diversity
# ---------------------------------------------------------------------------


def _bray_curtis(a: list[float], b: list[float]) -> float:
    """Compute Bray-Curtis dissimilarity.

    BC = sum(|a_i - b_i|) / sum(a_i + b_i).

    Args:
        a: Abundances for sample A.
        b: Abundances for sample B.

    Returns:
        Bray-Curtis dissimilarity in [0, 1].
    """
    num = sum(abs(a[i] - b[i]) for i in range(len(a)))
    denom = sum(a[i] + b[i] for i in range(len(a)))
    return num / denom if denom > 0 else 0.0


def _jaccard_distance(a: list[float], b: list[float]) -> float:
    """Compute Jaccard distance (presence/absence).

    J = 1 - |A intersection B| / |A union B|.

    Args:
        a: Abundances for sample A.
        b: Abundances for sample B.

    Returns:
        Jaccard distance in [0, 1].
    """
    present_a = set(i for i, v in enumerate(a) if v > 0)
    present_b = set(i for i, v in enumerate(b) if v > 0)

    intersection = len(present_a & present_b)
    union = len(present_a | present_b)

    return 1.0 - intersection / union if union > 0 else 0.0


def _clr_transform_samples(samples: list[list[float]]) -> list[list[float]]:
    """Apply centered log-ratio transform to all samples.

    Args:
        samples: List of abundance vectors.

    Returns:
        CLR-transformed samples.
    """
    pseudocount = 0.5
    transformed: list[list[float]] = []

    for sample in samples:
        # Add pseudocount and take log
        logged = [math.log(max(v, 0) + pseudocount) for v in sample]
        mean_log = sum(logged) / len(logged) if logged else 0.0
        clr = [v - mean_log for v in logged]
        transformed.append(clr)

    return transformed


def _euclidean(a: list[float], b: list[float]) -> float:
    """Compute Euclidean distance.

    Args:
        a: First vector.
        b: Second vector.

    Returns:
        Euclidean distance.
    """
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(len(a))))


# ---------------------------------------------------------------------------
# Private helpers: PERMANOVA
# ---------------------------------------------------------------------------


def _compute_pseudo_f(
    dist_matrix: list[list[float]],
    groups: list[str],
    unique_groups: list[str],
) -> float:
    """Compute pseudo-F statistic for PERMANOVA.

    F = (SS_between / (a - 1)) / (SS_within / (N - a))
    where a = number of groups, N = total samples.

    Args:
        dist_matrix: Distance matrix.
        groups: Group labels.
        unique_groups: Sorted unique group labels.

    Returns:
        Pseudo-F statistic. With no within-group variation (``SS_within`` of zero) it is
        ``inf`` when the groups differ and NaN when every distance is zero; it is never
        replaced by 0, which would give a perfectly separated design a p-value of 1.
    """
    n = len(groups)
    a = len(unique_groups)

    ss_total = _compute_ss_total(dist_matrix)
    ss_within = _compute_ss_within(dist_matrix, groups, unique_groups)
    ss_between = max(ss_total - ss_within, 0.0)

    if ss_within > 0:
        return (ss_between / (a - 1)) / (ss_within / (n - a))
    return math.inf if ss_between > 0 else math.nan


def _compute_ss_total(dist_matrix: list[list[float]]) -> float:
    """Compute total sum of squares from distance matrix.

    SS_total = sum(d_ij^2) / N.

    Args:
        dist_matrix: Distance matrix.

    Returns:
        Total sum of squared distances.
    """
    n = len(dist_matrix)
    if n == 0:
        return 0.0

    total = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            total += dist_matrix[i][j] ** 2

    return total / n


def _compute_ss_within(
    dist_matrix: list[list[float]],
    groups: list[str],
    unique_groups: list[str],
) -> float:
    """Compute within-group sum of squares.

    Args:
        dist_matrix: Distance matrix.
        groups: Group labels.
        unique_groups: Sorted unique group labels.

    Returns:
        Within-group sum of squared distances.
    """
    ss_within = 0.0
    for grp in unique_groups:
        grp_indices = [i for i, g in enumerate(groups) if g == grp]
        n_grp = len(grp_indices)
        if n_grp <= 1:
            continue

        grp_ss = 0.0
        for ii in range(len(grp_indices)):
            for jj in range(ii + 1, len(grp_indices)):
                i = grp_indices[ii]
                j = grp_indices[jj]
                grp_ss += dist_matrix[i][j] ** 2

        ss_within += grp_ss / n_grp

    return ss_within


# ---------------------------------------------------------------------------
# Private helpers: Ordination
# ---------------------------------------------------------------------------

# Relative size below which an eigenvalue, or the asymmetry of a distance matrix, is rounding noise
_EIGEN_NOISE = 1e-10

# SMACOF settings for NMDS: iteration limit and the relative stress improvement that counts as converged
_NMDS_MAX_ITER = 500
_NMDS_TOLERANCE = 1e-6


def _validate_distance_matrix(distance_matrix: Any) -> list[list[float]]:
    """Return ``distance_matrix`` as float rows after checking that it is a distance matrix.

    Args:
        distance_matrix: Square matrix as nested lists, a numpy array or a pandas frame.

    Returns:
        The matrix as a list of lists of floats.

    Raises:
        ValueError: If the matrix is empty, not square, has a non-finite or negative entry,
            has a non-zero diagonal, or is not symmetric (beyond rounding noise).
    """
    rows = _as_list(distance_matrix)
    n = len(rows)
    if n == 0:
        raise ValueError("Distance matrix must not be empty")

    matrix: list[list[float]] = []
    for row in rows:
        if isinstance(row, (str, bytes)) or not hasattr(row, "__iter__"):
            raise ValueError("Distance matrix must be square (a 2-D matrix)")
        cells = _as_list(row)
        if len(cells) != n:
            raise ValueError("Distance matrix must be square")
        matrix.append([float(value) for value in cells])

    scale = 1.0
    for i in range(n):
        for j in range(n):
            value = matrix[i][j]
            if not math.isfinite(value):
                raise ValueError(f"Distance matrix entry ({i}, {j}) is not finite: {value!r}")
            if value < 0:
                raise ValueError(f"Distance matrix entry ({i}, {j}) is negative: {value!r}")
            scale = max(scale, value)

    tolerance = _EIGEN_NOISE * scale
    for i in range(n):
        if matrix[i][i] > tolerance:
            raise ValueError(f"Distance matrix must have a zero diagonal, but entry ({i}, {i}) is {matrix[i][i]!r}")
        for j in range(i + 1, n):
            if abs(matrix[i][j] - matrix[j][i]) > tolerance:
                raise ValueError(
                    f"Distance matrix must be symmetric, but entries ({i}, {j}) and ({j}, {i}) are "
                    f"{matrix[i][j]!r} and {matrix[j][i]!r}"
                )
    return matrix


def _jacobi_eigh(matrix: list[list[float]], max_sweeps: int = 100) -> tuple[list[float], list[list[float]]]:
    """Eigendecomposition of a real symmetric matrix by cyclic Jacobi rotations (pure Python).

    Unlike power iteration, which finds the eigenvalue of largest *magnitude* one at a time
    (a large negative eigenvalue of a non-Euclidean PCoA matrix would be returned first),
    Jacobi diagonalises the whole matrix, so every eigenvalue, negative ones included, is
    exact to rounding.

    Args:
        matrix: Symmetric matrix (not modified).
        max_sweeps: Upper bound on full sweeps over the off-diagonal (a handful suffice).

    Returns:
        ``(eigenvalues, eigenvectors)`` with eigenvalues ascending and ``eigenvectors[i][k]`` the
        i-th component of the eigenvector of eigenvalue ``k``.
    """
    n = len(matrix)
    a = [list(row) for row in matrix]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    for _ in range(max_sweeps):
        off_diagonal = sum(a[i][j] ** 2 for i in range(n) for j in range(i + 1, n))
        diagonal = sum(a[i][i] ** 2 for i in range(n))
        if off_diagonal <= 1e-28 * (diagonal + off_diagonal):
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                a_pq = a[p][q]
                if a_pq == 0.0:
                    continue
                theta = (a[q][q] - a[p][p]) / (2.0 * a_pq)
                t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):  # columns p, q of A
                    a_kp, a_kq = a[k][p], a[k][q]
                    a[k][p] = c * a_kp - s * a_kq
                    a[k][q] = s * a_kp + c * a_kq
                for k in range(n):  # rows p, q of A
                    a_pk, a_qk = a[p][k], a[q][k]
                    a[p][k] = c * a_pk - s * a_qk
                    a[q][k] = s * a_pk + c * a_qk
                for k in range(n):  # accumulate the eigenvectors
                    v_kp, v_kq = v[k][p], v[k][q]
                    v[k][p] = c * v_kp - s * v_kq
                    v[k][q] = s * v_kp + c * v_kq

    order = sorted(range(n), key=lambda k: a[k][k])
    return [a[k][k] for k in order], [[v[i][k] for k in order] for i in range(n)]


def _symmetric_eigh(matrix: list[list[float]]) -> tuple[list[float], list[list[float]]]:
    """Eigenvalues (descending) and eigenvectors of a real symmetric matrix.

    Uses ``numpy.linalg.eigh`` when NumPy is installed and ``_jacobi_eigh`` otherwise. The sign
    of each eigenvector is fixed (its largest-magnitude component is positive) because the
    sign is arbitrary and differs between LAPACK builds.

    Args:
        matrix: Symmetric matrix.

    Returns:
        ``(eigenvalues, eigenvectors)``, eigenvalues descending and ``eigenvectors[i][k]`` the
        i-th component of the eigenvector of eigenvalue ``k``.
    """
    n = len(matrix)
    if HAS_NUMPY:
        values_np, vectors_np = np.linalg.eigh(np.array(matrix, dtype=np.float64))
        values = [float(x) for x in values_np[::-1]]
        vectors = vectors_np[:, ::-1].tolist()
    else:
        ascending, ascending_vectors = _jacobi_eigh(matrix)
        values = ascending[::-1]
        vectors = [row[::-1] for row in ascending_vectors]

    for k in range(n):
        pivot = max(range(n), key=lambda i: abs(vectors[i][k]))
        if vectors[pivot][k] < 0:
            for i in range(n):
                vectors[i][k] = -vectors[i][k]
    return values, vectors


def _double_centre(matrix: list[list[float]]) -> list[list[float]]:
    """Gower's double-centred matrix ``B = -1/2 J D**2 J`` of a symmetric distance matrix."""
    n = len(matrix)
    squared = [[value * value for value in row] for row in matrix]
    row_means = [sum(row) / n for row in squared]
    grand_mean = sum(row_means) / n
    return [[-0.5 * (squared[i][j] - row_means[i] - row_means[j] + grand_mean) for j in range(n)] for i in range(n)]


def _pcoa(
    distance_matrix: list[list[float]],
    n_components: int,
    warn_non_euclidean: bool = True,
) -> dict:
    """Principal Coordinates Analysis (classical multidimensional scaling).

    Eigendecomposes the double-centred squared distance matrix (``numpy.linalg.eigh``, or exact
    Jacobi rotations without NumPy). The whole spectrum is kept so that negative eigenvalues,
    which arise whenever the distances are not Euclidean, are reported instead of being
    dropped.

    Args:
        distance_matrix: Validated square distance matrix.
        n_components: Number of output dimensions.
        warn_non_euclidean: Log a warning when there are negative eigenvalues (NMDS uses PCoA only
            as a starting point, where they do not matter).

    Returns:
        Dictionary with coordinates, eigenvalues, variance_explained, all_eigenvalues,
        n_negative_eigenvalues and negative_eigenvalue_fraction (see ``ordination``).
    """
    n = len(distance_matrix)
    actual_components = min(n_components, n - 1) if n > 1 else 1

    eigenvalues, eigenvectors = _symmetric_eigh(_double_centre(distance_matrix))

    noise = _EIGEN_NOISE * max((abs(e) for e in eigenvalues), default=0.0)
    positive_mass = sum(e for e in eigenvalues if e > noise)
    negative_mass = -sum(e for e in eigenvalues if e < -noise)
    n_negative = sum(1 for e in eigenvalues if e < -noise)
    total_mass = positive_mass + negative_mass
    negative_fraction = negative_mass / total_mass if total_mass > 0 else 0.0

    # An axis needs a positive eigenvalue to have real coordinates; the others stay at zero
    coords = [
        [
            eigenvectors[i][k] * math.sqrt(eigenvalues[k]) if eigenvalues[k] > noise else 0.0
            for k in range(actual_components)
        ]
        for i in range(n)
    ]
    retained = eigenvalues[:actual_components]
    variance_explained = [e / positive_mass if e > noise and positive_mass > 0 else 0.0 for e in retained]

    if n_negative and warn_non_euclidean:
        logger.warning(
            f"PCoA: the distance matrix is not Euclidean; {n_negative} negative eigenvalues hold "
            f"{negative_fraction:.1%} of the total absolute eigenvalue mass, and the axes represent only "
            "the Euclidean part"
        )
    logger.info(
        f"PCoA complete: {actual_components} components, "
        f"variance explained = {[f'{v:.3f}' for v in variance_explained]}"
    )

    return {
        "coordinates": coords,
        "eigenvalues": retained,
        "variance_explained": variance_explained,
        "all_eigenvalues": eigenvalues,
        "n_negative_eigenvalues": n_negative,
        "negative_eigenvalue_fraction": negative_fraction,
    }


def _isotonic_fit(values: list[float]) -> list[float]:
    """Least-squares non-decreasing fit of ``values`` (pool adjacent violators).

    Args:
        values: Sequence to fit, in the order its fit must be non-decreasing.

    Returns:
        The isotonic regression of ``values``, same length.
    """
    sums: list[float] = []
    counts: list[int] = []
    for value in values:
        sums.append(value)
        counts.append(1)
        # pool the last two blocks while their means decrease (cross-multiplied to avoid divisions)
        while len(sums) > 1 and sums[-2] * counts[-1] > sums[-1] * counts[-2]:
            pooled_sum = sums.pop()
            pooled_count = counts.pop()
            sums[-1] += pooled_sum
            counts[-1] += pooled_count
    fitted: list[float] = []
    for total, count in zip(sums, counts):
        fitted.extend([total / count] * count)
    return fitted


class _PairIndex:
    """The ``n (n - 1) / 2`` unordered pairs of an NMDS problem and the order of their dissimilarities.

    Pairs are sorted by dissimilarity once. Pairs whose dissimilarities tie form blocks inside
    which the order is free (Kruskal's primary approach: tied dissimilarities need not map to
    equal distances), so each iteration only reorders those blocks by the current distances.
    """

    def __init__(self, dissimilarities: list[list[float]]) -> None:
        n = len(dissimilarities)
        self.n = n
        self.i = [i for i in range(n - 1) for _ in range(i + 1, n)]
        self.j = [j for i in range(n - 1) for j in range(i + 1, n)]
        self.delta = [dissimilarities[i][j] for i, j in zip(self.i, self.j)]
        self.size = len(self.delta)
        self.base_order = sorted(range(self.size), key=self.delta.__getitem__)
        self.tie_blocks: list[tuple[int, int]] = []
        start = 0
        for position in range(1, self.size + 1):
            if position == self.size or self.delta[self.base_order[position]] != self.delta[self.base_order[start]]:
                if position - start > 1:
                    self.tie_blocks.append((start, position))
                start = position

    def distances(self, config: list[list[float]]) -> list[float]:
        """Euclidean distances of the configuration for every pair."""
        return [math.dist(config[i], config[j]) for i, j in zip(self.i, self.j)]

    def disparities(self, distances: list[float]) -> list[float]:
        """Isotonic regression of ``distances`` on the dissimilarities (the disparities)."""
        order = list(self.base_order)
        for lo, hi in self.tie_blocks:
            order[lo:hi] = sorted(order[lo:hi], key=distances.__getitem__)
        fitted = _isotonic_fit([distances[k] for k in order])
        disparities = [0.0] * self.size
        for position, k in enumerate(order):
            disparities[k] = fitted[position]
        return disparities


def _kruskal_stress1(distances: list[float], disparities: list[float]) -> float:
    """Kruskal's stress-1, ``sqrt(sum((d - d_hat)**2) / sum(d**2))``."""
    total = sum(d * d for d in distances)
    if total <= 0.0:
        return 0.0
    return math.sqrt(sum((d - h) ** 2 for d, h in zip(distances, disparities)) / total)


def _nmds_stress(pairs: _PairIndex, config: list[list[float]]) -> float:
    """Kruskal's stress-1 of a configuration, with its optimal monotone disparities."""
    distances = pairs.distances(config)
    return _kruskal_stress1(distances, pairs.disparities(distances))


def _smacof_nmds(
    pairs: _PairIndex,
    config: list[list[float]],
    max_iter: int,
    tolerance: float,
) -> tuple[list[list[float]], float, int, bool]:
    """Minimise Kruskal's stress-1 from one starting configuration (SMACOF with isotonic disparities).

    Each iteration fits the disparities (isotonic regression of the current distances on the
    dissimilarities), rescales them to a fixed sum of squares so the configuration cannot
    shrink to nothing, and applies the Guttman transform, which never increases the raw
    stress for fixed disparities.

    Args:
        pairs: Pair index of the dissimilarities.
        config: Starting configuration (n x p), not modified.
        max_iter: Iteration limit.
        tolerance: Stop when the stress improves by less than this fraction of itself.

    Returns:
        ``(configuration, stress, iterations, converged)`` for the best configuration visited.
    """
    n = pairs.n
    dim = len(config[0])
    current = [list(row) for row in config]
    best = current
    best_stress = math.inf
    previous = math.inf
    converged = False
    iterations = 0

    for iterations in range(1, max_iter + 1):
        distances = pairs.distances(current)
        disparities = pairs.disparities(distances)
        stress = _kruskal_stress1(distances, disparities)
        if stress < best_stress:
            best, best_stress = current, stress
        if stress <= 1e-12 or (previous < math.inf and previous - stress <= tolerance * previous):
            converged = True
            break
        previous = stress

        total_sq = sum(h * h for h in disparities)
        if total_sq <= 0.0:
            break
        scale = math.sqrt(pairs.size / total_sq)
        updated = [[0.0] * dim for _ in range(n)]
        for i, j, d, h in zip(pairs.i, pairs.j, distances, disparities):
            if d > 0.0:
                factor = scale * h / d
                xi, xj, ui, uj = current[i], current[j], updated[i], updated[j]
                for c in range(dim):
                    step = factor * (xi[c] - xj[c])
                    ui[c] += step
                    uj[c] -= step
        current = [[value / n for value in row] for row in updated]

    return best, best_stress, iterations, converged


def _principal_axes(config: list[list[float]]) -> list[list[float]]:
    """Centre a configuration and rotate it to its principal axes (largest variance first, sign fixed)."""
    n = len(config)
    dim = len(config[0])
    means = [sum(row[c] for row in config) / n for c in range(dim)]
    centred = [[row[c] - means[c] for c in range(dim)] for row in config]
    covariance = [[sum(row[a] * row[b] for row in centred) for b in range(dim)] for a in range(dim)]
    _, axes = _symmetric_eigh(covariance)
    return [[sum(row[a] * axes[a][k] for a in range(dim)) for k in range(dim)] for row in centred]


def _nmds(
    distance_matrix: list[list[float]],
    n_components: int,
    seed: int | None = 42,
    n_init: int = 5,
) -> dict:
    """Kruskal's non-metric multidimensional scaling with isotonic regression.

    Runs SMACOF (see ``_smacof_nmds``) from a classical-scaling start and from ``n_init``
    random starts drawn from a generator seeded with ``seed``, and returns the configuration
    with the lowest stress-1.

    Args:
        distance_matrix: Validated square distance matrix.
        n_components: Number of output dimensions.
        seed: Seed of the random starts (None draws fresh entropy).
        n_init: Number of random starts besides the classical one.

    Returns:
        Dictionary with coordinates, stress, stress_by_start, n_iter and converged
        (see ``ordination``).

    Raises:
        ValueError: If there are fewer than 3 samples or every dissimilarity is zero.
    """
    n = len(distance_matrix)
    if n < 3:
        raise ValueError(f"NMDS needs at least 3 samples, got {n}")

    pairs = _PairIndex(distance_matrix)
    if not any(value > 0.0 for value in pairs.delta):
        raise ValueError("All dissimilarities are zero; NMDS is undefined")

    rng = random.Random(seed)
    mean_delta_sq = sum(value * value for value in pairs.delta) / pairs.size

    # Start 0: classical scaling. Axes without a positive eigenvalue are filled with small noise
    # so that the starting points do not lie in a lower-dimensional subspace.
    classical = _pcoa(distance_matrix, n_components, warn_non_euclidean=False)["coordinates"]
    jitter = 1e-3 * math.sqrt(mean_delta_sq)
    starts = [
        [
            [value + rng.gauss(0.0, jitter) for value in row]
            + [rng.gauss(0.0, jitter) for _ in range(n_components - len(row))]
            for row in classical
        ]
    ]
    for _ in range(n_init):
        starts.append([[rng.gauss(0.0, 1.0) for _ in range(n_components)] for _ in range(n)])

    outcomes = [_smacof_nmds(pairs, start, _NMDS_MAX_ITER, _NMDS_TOLERANCE) for start in starts]
    stress_by_start = [outcome[1] for outcome in outcomes]
    best = min(outcomes, key=lambda outcome: outcome[1])  # the earliest start wins a tie

    coordinates = _principal_axes(best[0])
    # Scale so the configuration's root-mean-square distance equals that of the dissimilarities
    mean_d_sq = sum(d * d for d in pairs.distances(coordinates)) / pairs.size
    if mean_d_sq > 0.0:
        factor = math.sqrt(mean_delta_sq / mean_d_sq)
        coordinates = [[value * factor for value in row] for row in coordinates]

    stress = _nmds_stress(pairs, coordinates)
    logger.info(
        f"NMDS complete: {len(starts)} starts, best stress-1={stress:.4f} after {best[2]} iterations "
        f"({'converged' if best[3] else 'iteration limit'})"
    )

    return {
        "coordinates": coordinates,
        "stress": stress,
        "stress_by_start": stress_by_start,
        "n_iter": best[2],
        "converged": best[3],
    }


__all__ = [
    "alpha_diversity",
    "beta_diversity",
    "ordination",
    "permanova",
    "rarefaction_curve",
    "rarefy",
]
