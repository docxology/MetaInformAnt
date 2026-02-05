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


def alpha_diversity(
    abundances: list[float] | list[int],
    metric: str = "shannon",
) -> dict:
    """Compute alpha diversity for a single sample.

    Calculates within-sample diversity using the specified metric. Input
    abundances can be raw counts or relative abundances.

    Args:
        abundances: Abundance values for each taxon/OTU in the sample.
            For count-based metrics (chao1, ace, observed), integer counts
            are expected. For entropy-based metrics (shannon, simpson),
            both counts and proportions are accepted.
        metric: Diversity metric. One of "shannon", "simpson",
            "invsimpson", "chao1", "ace", "observed", "fisher_alpha",
            "pielou_evenness".

    Returns:
        Dictionary with keys:
            - value: float computed diversity value.
            - metric: str metric name used.
            - n_species: int number of species with non-zero abundance.
            - total_count: float sum of all abundances.

    Raises:
        ValueError: If metric is not recognized or abundances is empty.
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

    if not abundances:
        raise ValueError("abundances must not be empty")

    # Filter to non-zero values
    nonzero = [a for a in abundances if a > 0]
    n_species = len(nonzero)
    total = sum(nonzero)

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
    if not abundances:
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
            - pseudo_f: float observed pseudo-F statistic.
            - p_value: float permutation p-value.
            - r_squared: float fraction of variance explained by groups.
            - n_permutations: int number of permutations used.

    Raises:
        ValueError: If dimensions are inconsistent or fewer than 2 groups.
    """
    n = len(distance_matrix)
    if len(groups) != n:
        raise ValueError(f"groups length ({len(groups)}) must match distance matrix " f"dimension ({n})")

    unique_groups = sorted(set(groups))
    if len(unique_groups) < 2:
        raise ValueError("At least 2 groups required for PERMANOVA")

    logger.info(f"Running PERMANOVA: {n} samples, {len(unique_groups)} groups, " f"{n_permutations} permutations")

    # Compute observed pseudo-F
    observed_f = _compute_pseudo_f(distance_matrix, groups, unique_groups)

    # Permutation test
    rng = random.Random(seed)
    n_greater = 0

    for _ in range(n_permutations):
        perm_groups = list(groups)
        rng.shuffle(perm_groups)
        perm_f = _compute_pseudo_f(distance_matrix, perm_groups, unique_groups)
        if perm_f >= observed_f:
            n_greater += 1

    p_value = (n_greater + 1) / (n_permutations + 1)

    # R-squared: SS_between / SS_total
    ss_total = _compute_ss_total(distance_matrix)
    ss_within = _compute_ss_within(distance_matrix, groups, unique_groups)
    ss_between = ss_total - ss_within
    r_squared = ss_between / ss_total if ss_total > 0 else 0.0

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
) -> dict:
    """Perform ordination on a distance matrix.

    Reduces the distance matrix to a low-dimensional representation
    for visualization.

    Args:
        distance_matrix: Square symmetric distance matrix.
        method: Ordination method. "pcoa" (principal coordinates
            analysis) or "nmds" (non-metric multidimensional scaling).
        n_components: Number of output dimensions.

    Returns:
        Dictionary with keys:
            - coordinates: list[list[float]] of sample coordinates
              (n_samples x n_components).
            - eigenvalues: list[float] eigenvalues (pcoa) or stress values
              (nmds) per component.
            - variance_explained: list[float] fraction of variance
              explained by each component (pcoa only).

    Raises:
        ValueError: If method is not recognized or matrix is not square.
    """
    valid_methods = ("pcoa", "nmds")
    if method not in valid_methods:
        raise ValueError(f"Invalid method '{method}'. Must be one of {valid_methods}")

    n = len(distance_matrix)
    if any(len(row) != n for row in distance_matrix):
        raise ValueError("Distance matrix must be square")

    logger.info(f"Running ordination ({method}): {n} samples, {n_components} components")

    if method == "pcoa":
        return _pcoa(distance_matrix, n_components)
    else:
        return _nmds(distance_matrix, n_components)


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
        ACE estimated richness.
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
        return float(len(counts))

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

    Uses Newton's method to solve: S = alpha * ln(1 + N/alpha).

    Args:
        counts: Non-zero integer counts per species.

    Returns:
        Fisher's alpha parameter.
    """
    s = len(counts)
    n = sum(counts)

    if s == 0 or n == 0:
        return 0.0

    # Initial estimate
    alpha = s / math.log(n) if n > 1 else float(s)

    # Newton's method
    for _ in range(100):
        if alpha <= 0:
            alpha = 0.01
        f_val = alpha * math.log(1.0 + n / alpha) - s
        # Derivative: ln(1 + N/alpha) - N/(alpha + N)
        f_prime = math.log(1.0 + n / alpha) - n / (alpha + n)

        if abs(f_prime) < 1e-15:
            break

        alpha_new = alpha - f_val / f_prime
        if abs(alpha_new - alpha) < 1e-8:
            break
        alpha = alpha_new

    return max(alpha, 0.0)


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
        Pseudo-F statistic.
    """
    n = len(groups)
    a = len(unique_groups)

    if a <= 1 or n <= a:
        return 0.0

    ss_total = _compute_ss_total(dist_matrix)
    ss_within = _compute_ss_within(dist_matrix, groups, unique_groups)
    ss_between = ss_total - ss_within

    f_stat = (ss_between / (a - 1)) / (ss_within / (n - a)) if ss_within > 0 else 0.0
    return f_stat


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

        for ii in range(len(grp_indices)):
            for jj in range(ii + 1, len(grp_indices)):
                i = grp_indices[ii]
                j = grp_indices[jj]
                ss_within += dist_matrix[i][j] ** 2

        ss_within /= 1  # Already summed; normalize per group below

    # Normalize: sum(d^2) / n_i for each group
    # Re-compute with proper normalization
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


def _pcoa(
    distance_matrix: list[list[float]],
    n_components: int,
) -> dict:
    """Principal Coordinates Analysis (classical multidimensional scaling).

    Performs eigendecomposition of the double-centered distance matrix
    to find coordinates that preserve the original distances.

    Args:
        distance_matrix: Square distance matrix.
        n_components: Number of output dimensions.

    Returns:
        Dictionary with coordinates, eigenvalues, and variance_explained.
    """
    n = len(distance_matrix)
    actual_components = min(n_components, n - 1) if n > 1 else 1

    # Step 1: Square the distances
    D_sq = [[distance_matrix[i][j] ** 2 for j in range(n)] for i in range(n)]

    # Step 2: Double centering
    # B = -0.5 * J * D^2 * J, where J = I - 1/n * ones
    row_means = [sum(D_sq[i]) / n for i in range(n)]
    col_means = [sum(D_sq[i][j] for i in range(n)) / n for j in range(n)]
    grand_mean = sum(sum(row) for row in D_sq) / (n * n)

    B = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            B[i][j] = -0.5 * (D_sq[i][j] - row_means[i] - col_means[j] + grand_mean)

    # Step 3: Eigendecomposition
    if HAS_NUMPY:
        B_np = np.array(B, dtype=np.float64)
        eigenvalues, eigenvectors = np.linalg.eigh(B_np)
        # Sort by decreasing eigenvalue
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        # Take top components
        evals = eigenvalues[:actual_components].tolist()
        coords = []
        for i in range(n):
            row = []
            for j in range(actual_components):
                ev = max(eigenvalues[j], 0.0)
                row.append(float(eigenvectors[i, j] * math.sqrt(ev)))
            coords.append(row)
    else:
        # Pure Python fallback: power iteration for top eigenvectors
        coords, evals = _power_iteration_eigen(B, actual_components)

    # Compute variance explained
    total_var = sum(max(e, 0.0) for e in evals)
    variance_explained = [max(e, 0.0) / total_var if total_var > 0 else 0.0 for e in evals]

    logger.info(
        f"PCoA complete: {actual_components} components, "
        f"variance explained = {[f'{v:.3f}' for v in variance_explained]}"
    )

    return {
        "coordinates": coords,
        "eigenvalues": evals,
        "variance_explained": variance_explained,
    }


def _nmds(
    distance_matrix: list[list[float]],
    n_components: int,
    max_iter: int = 300,
    seed: int | None = None,
) -> dict:
    """Non-metric multidimensional scaling (NMDS).

    Iteratively optimizes coordinates to preserve the rank order of
    pairwise distances. Uses stress as the goodness-of-fit metric.

    Args:
        distance_matrix: Square distance matrix.
        n_components: Number of output dimensions.
        max_iter: Maximum iterations.
        seed: Random seed.

    Returns:
        Dictionary with coordinates, eigenvalues (stress per iteration),
        and variance_explained.
    """
    n = len(distance_matrix)
    rng = random.Random(seed)

    # Initialize with random coordinates
    coords = [[rng.gauss(0, 1) for _ in range(n_components)] for _ in range(n)]

    stress_history: list[float] = []

    for iteration in range(max_iter):
        # Compute current distances
        current_dists = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                d = _euclidean(coords[i], coords[j])
                current_dists[i][j] = d
                current_dists[j][i] = d

        # Compute stress
        stress_num = 0.0
        stress_denom = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                delta = distance_matrix[i][j]
                d = current_dists[i][j]
                stress_num += (d - delta) ** 2
                stress_denom += delta**2

        stress = math.sqrt(stress_num / stress_denom) if stress_denom > 0 else 0.0
        stress_history.append(stress)

        # Check convergence
        if iteration > 0 and abs(stress_history[-1] - stress_history[-2]) < 1e-6:
            break

        # Gradient step: move points to reduce stress
        learning_rate = 0.1 / (1 + iteration * 0.01)
        new_coords = [list(row) for row in coords]

        for i in range(n):
            for dim in range(n_components):
                gradient = 0.0
                for j in range(n):
                    if i == j:
                        continue
                    d = current_dists[i][j]
                    delta = distance_matrix[i][j]
                    if d > 0:
                        gradient += (d - delta) / d * (coords[i][dim] - coords[j][dim])

                new_coords[i][dim] -= learning_rate * gradient

        coords = new_coords

    final_stress = stress_history[-1] if stress_history else 0.0
    logger.info(f"NMDS complete: {len(stress_history)} iterations, " f"final stress={final_stress:.4f}")

    return {
        "coordinates": coords,
        "eigenvalues": stress_history[-min(n_components, len(stress_history)) :],
        "variance_explained": [1.0 - final_stress] * n_components,
    }


def _power_iteration_eigen(
    matrix: list[list[float]],
    n_components: int,
    max_iter: int = 100,
) -> tuple[list[list[float]], list[float]]:
    """Extract top eigenvectors via power iteration (pure Python fallback).

    Args:
        matrix: Symmetric matrix.
        n_components: Number of eigenvectors to extract.
        max_iter: Maximum iterations per eigenvector.

    Returns:
        Tuple of (coordinates, eigenvalues).
    """
    n = len(matrix)
    eigenvalues: list[float] = []
    eigenvectors: list[list[float]] = []
    work_matrix = [list(row) for row in matrix]

    for comp in range(n_components):
        # Random initialization
        rng = random.Random(42 + comp)
        vec = [rng.gauss(0, 1) for _ in range(n)]
        norm = math.sqrt(sum(v * v for v in vec))
        vec = [v / norm for v in vec]

        eigenvalue = 0.0
        for _ in range(max_iter):
            # Matrix-vector multiplication
            new_vec = [sum(work_matrix[i][j] * vec[j] for j in range(n)) for i in range(n)]

            # Compute eigenvalue (Rayleigh quotient)
            eigenvalue = sum(new_vec[i] * vec[i] for i in range(n))

            # Normalize
            norm = math.sqrt(sum(v * v for v in new_vec))
            if norm == 0:
                break
            vec = [v / norm for v in new_vec]

        eigenvalues.append(eigenvalue)
        eigenvectors.append(vec)

        # Deflate matrix
        for i in range(n):
            for j in range(n):
                work_matrix[i][j] -= eigenvalue * vec[i] * vec[j]

    # Convert eigenvectors to coordinates
    coords: list[list[float]] = []
    for i in range(n):
        row = []
        for comp in range(n_components):
            ev = max(eigenvalues[comp], 0.0)
            row.append(eigenvectors[comp][i] * math.sqrt(ev))
        coords.append(row)

    return coords, eigenvalues


__all__ = [
    "alpha_diversity",
    "beta_diversity",
    "ordination",
    "permanova",
    "rarefaction_curve",
    "rarefy",
]
