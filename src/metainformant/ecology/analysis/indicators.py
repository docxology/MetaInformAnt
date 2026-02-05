"""Ecological indicator and multivariate community analysis.

This module implements classical ecological statistical methods for analyzing
community composition across groups, including indicator species analysis,
permutational hypothesis testing, hierarchical clustering, and similarity
decomposition. All algorithms are pure-Python implementations following the
original published methods.

References:
    - Dufrene, M. & Legendre, P. (1997). Species assemblages and indicator species.
    - Clarke, K.R. (1993). Non-parametric multivariate analyses of changes in community structure.
    - Anderson, M.J. (2001). A new non-parametric multivariate analysis of variance.
    - Clarke, K.R. (1993). SIMPER: Similarity Percentages.
    - Anderson, M.J. (2006). Distance-based tests for homogeneity of multivariate dispersions.
"""

from __future__ import annotations

import math
import random
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core import errors, logging, validation

logger = logging.get_logger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _validate_square_matrix(matrix: List[List[float]], name: str = "distance_matrix") -> int:
    """Validate that a matrix is square and symmetric-ish.

    Args:
        matrix: The matrix to validate.
        name: Human-readable name for error messages.

    Returns:
        The dimension (number of rows/columns).

    Raises:
        ValueError: If the matrix is empty, non-square, or contains negative values.
    """
    validation.validate_not_empty(matrix, name)
    n = len(matrix)
    for i, row in enumerate(matrix):
        if len(row) != n:
            raise ValueError(f"{name} must be square: row {i} has {len(row)} elements, expected {n}")
        for j, val in enumerate(row):
            if val < 0:
                raise ValueError(f"{name} contains negative value at ({i}, {j}): {val}")
    return n


def _validate_labels_match(labels: List[str], n: int, name: str = "group_labels") -> None:
    """Validate that label count matches the expected dimension.

    Args:
        labels: Group labels to validate.
        n: Expected number of labels.
        name: Human-readable name for error messages.

    Raises:
        ValueError: If the number of labels does not match *n*.
    """
    validation.validate_not_empty(labels, name)
    if len(labels) != n:
        raise ValueError(f"{name} length ({len(labels)}) must match matrix dimension ({n})")


def _unique_groups(labels: List[str]) -> List[str]:
    """Return unique group labels preserving first-occurrence order."""
    seen: set[str] = set()
    groups: List[str] = []
    for label in labels:
        if label not in seen:
            seen.add(label)
            groups.append(label)
    return groups


def _bray_curtis_single(a: List[float], b: List[float]) -> float:
    """Compute Bray-Curtis dissimilarity between two abundance vectors.

    Args:
        a: First abundance vector.
        b: Second abundance vector.

    Returns:
        Bray-Curtis dissimilarity in [0, 1].
    """
    numerator = sum(abs(x - y) for x, y in zip(a, b))
    denominator = sum(x + y for x, y in zip(a, b))
    return numerator / denominator if denominator > 0 else 0.0


def _permute_labels(labels: List[str], rng: random.Random) -> List[str]:
    """Return a randomly shuffled copy of *labels* using *rng*."""
    shuffled = list(labels)
    rng.shuffle(shuffled)
    return shuffled


# ---------------------------------------------------------------------------
# 1. Indicator Value (IndVal) Analysis
# ---------------------------------------------------------------------------


def indval(
    abundance_matrix: List[List[float]],
    group_labels: List[str],
    *,
    n_permutations: int = 999,
    seed: int = 42,
) -> List[Dict[str, Any]]:
    """Indicator Value (IndVal) analysis (Dufrene & Legendre, 1997).

    IndVal combines a species' *specificity* (relative abundance concentration
    in a target group) and *fidelity* (frequency of occurrence in the target
    group) to identify indicator species for each habitat or treatment.

    IndVal = specificity * fidelity * 100

    Args:
        abundance_matrix: Sites-by-species matrix where ``abundance_matrix[i][j]``
            is the abundance of species *j* at site *i*.
        group_labels: Group membership for each site (length must equal the
            number of rows in *abundance_matrix*).
        n_permutations: Number of permutations for the significance test.
        seed: Random seed for reproducibility.

    Returns:
        A list of dicts (one per species) each containing:
            - ``species_idx``: Column index of the species.
            - ``indval``: Maximum IndVal score (0--100).
            - ``specificity``: Specificity component for the best group.
            - ``fidelity``: Fidelity component for the best group.
            - ``best_group``: Group label with the highest IndVal.
            - ``p_value``: Permutation-based *p*-value.

    Raises:
        ValueError: If dimensions are inconsistent or inputs are empty.

    Example:
        >>> mat = [[10, 0, 3], [8, 1, 2], [0, 7, 1], [1, 9, 0]]
        >>> labels = ["forest", "forest", "grassland", "grassland"]
        >>> results = indval(mat, labels, n_permutations=99, seed=0)
        >>> results[0]["best_group"]
        'forest'
    """
    validation.validate_not_empty(abundance_matrix, "abundance_matrix")
    n_sites = len(abundance_matrix)
    _validate_labels_match(group_labels, n_sites, "group_labels")
    n_species = len(abundance_matrix[0])
    for i, row in enumerate(abundance_matrix):
        if len(row) != n_species:
            raise ValueError(
                f"All rows in abundance_matrix must have the same length; row {i} has {len(row)}, expected {n_species}"
            )

    groups = _unique_groups(group_labels)
    if len(groups) < 2:
        raise ValueError("At least two distinct groups are required for IndVal analysis")

    rng = random.Random(seed)
    logger.info("Running IndVal analysis: %d sites, %d species, %d groups", n_sites, n_species, len(groups))

    def _compute_indval_for_species(sp: int, labels: List[str]) -> Tuple[float, float, float, str]:
        """Return (indval, specificity, fidelity, best_group) for species *sp*."""
        # Group-level statistics
        group_abundances: Dict[str, List[float]] = defaultdict(list)
        for site_idx, grp in enumerate(labels):
            group_abundances[grp].append(abundance_matrix[site_idx][sp])

        best_iv = -1.0
        best_spec = 0.0
        best_fid = 0.0
        best_grp = groups[0]

        # Mean abundance across ALL groups (sum of group means)
        group_means = {g: (sum(vals) / len(vals)) if vals else 0.0 for g, vals in group_abundances.items()}
        total_mean = sum(group_means.values())

        for grp in groups:
            vals = group_abundances[grp]
            if not vals:
                continue

            # Specificity: mean abundance in target group / sum of mean abundances across all groups
            mean_in_group = group_means[grp]
            specificity = mean_in_group / total_mean if total_mean > 0 else 0.0

            # Fidelity: proportion of sites in target group where the species is present
            n_present = sum(1 for v in vals if v > 0)
            fidelity = n_present / len(vals)

            iv = specificity * fidelity * 100.0

            if iv > best_iv:
                best_iv = iv
                best_spec = specificity
                best_fid = fidelity
                best_grp = grp

        return best_iv, best_spec, best_fid, best_grp

    results: List[Dict[str, Any]] = []

    for sp in range(n_species):
        observed_iv, obs_spec, obs_fid, obs_grp = _compute_indval_for_species(sp, group_labels)

        # Permutation test
        n_geq = 0
        for _ in range(n_permutations):
            perm_labels = _permute_labels(group_labels, rng)
            perm_iv, _, _, _ = _compute_indval_for_species(sp, perm_labels)
            if perm_iv >= observed_iv:
                n_geq += 1

        p_value = (n_geq + 1) / (n_permutations + 1)

        results.append(
            {
                "species_idx": sp,
                "indval": round(observed_iv, 4),
                "specificity": round(obs_spec, 4),
                "fidelity": round(obs_fid, 4),
                "best_group": obs_grp,
                "p_value": round(p_value, 4),
            }
        )

    logger.info("IndVal analysis complete for %d species", n_species)
    return results


# ---------------------------------------------------------------------------
# 2. ANOSIM
# ---------------------------------------------------------------------------


def anosim(
    distance_matrix: List[List[float]],
    group_labels: List[str],
    *,
    n_permutations: int = 999,
    seed: int = 42,
) -> Dict[str, Any]:
    """Analysis of Similarities (Clarke, 1993).

    ANOSIM tests whether between-group dissimilarities are significantly
    larger than within-group dissimilarities using a rank-based R statistic
    and a permutation test.

    R = (mean_rank_between - mean_rank_within) / (N*(N-1)/4)

    Args:
        distance_matrix: Symmetric pairwise distance matrix (N x N).
        group_labels: Group membership for each observation (length N).
        n_permutations: Number of permutations for the significance test.
        seed: Random seed for reproducibility.

    Returns:
        Dictionary with keys:
            - ``r_statistic``: The ANOSIM R value (-1 to 1; values near 1
              indicate strong group separation).
            - ``p_value``: Proportion of permuted R values >= observed R.
            - ``n_permutations``: Number of permutations performed.

    Raises:
        ValueError: If matrix is not square or labels do not match.

    Example:
        >>> dm = [[0, 1, 5, 5], [1, 0, 5, 5], [5, 5, 0, 1], [5, 5, 1, 0]]
        >>> labels = ["A", "A", "B", "B"]
        >>> result = anosim(dm, labels, n_permutations=99, seed=0)
        >>> result["r_statistic"] > 0
        True
    """
    n = _validate_square_matrix(distance_matrix)
    _validate_labels_match(group_labels, n)

    rng = random.Random(seed)
    logger.info("Running ANOSIM: %d observations, %d permutations", n, n_permutations)

    # Collect all pairwise distances and rank them
    pairs: List[Tuple[int, int, float]] = []
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append((i, j, distance_matrix[i][j]))

    # Rank distances (average ranks for ties)
    sorted_indices = sorted(range(len(pairs)), key=lambda k: pairs[k][2])
    ranks = [0.0] * len(pairs)
    i = 0
    while i < len(sorted_indices):
        j = i
        while j < len(sorted_indices) and pairs[sorted_indices[j]][2] == pairs[sorted_indices[i]][2]:
            j += 1
        avg_rank = (i + j - 1) / 2.0 + 1.0  # 1-based average rank
        for k in range(i, j):
            ranks[sorted_indices[k]] = avg_rank
        i = j

    # Build a lookup: (i, j) -> rank  (i < j)
    rank_lookup: Dict[Tuple[int, int], float] = {}
    for idx, (ii, jj, _) in enumerate(pairs):
        rank_lookup[(ii, jj)] = ranks[idx]

    def _compute_r(labels: List[str]) -> float:
        within_ranks: List[float] = []
        between_ranks: List[float] = []
        for ii in range(n):
            for jj in range(ii + 1, n):
                r = rank_lookup[(ii, jj)]
                if labels[ii] == labels[jj]:
                    within_ranks.append(r)
                else:
                    between_ranks.append(r)

        if not within_ranks or not between_ranks:
            return 0.0

        mean_within = sum(within_ranks) / len(within_ranks)
        mean_between = sum(between_ranks) / len(between_ranks)
        total_pairs = n * (n - 1) / 2.0
        divisor = total_pairs / 2.0  # N*(N-1)/4

        return (mean_between - mean_within) / divisor if divisor > 0 else 0.0

    observed_r = _compute_r(group_labels)

    n_geq = 0
    for _ in range(n_permutations):
        perm_labels = _permute_labels(group_labels, rng)
        if _compute_r(perm_labels) >= observed_r:
            n_geq += 1

    p_value = (n_geq + 1) / (n_permutations + 1)

    logger.info("ANOSIM complete: R=%.4f, p=%.4f", observed_r, p_value)
    return {
        "r_statistic": round(observed_r, 6),
        "p_value": round(p_value, 4),
        "n_permutations": n_permutations,
    }


# ---------------------------------------------------------------------------
# 3. PERMANOVA
# ---------------------------------------------------------------------------


def permanova(
    distance_matrix: List[List[float]],
    group_labels: List[str],
    *,
    n_permutations: int = 999,
    seed: int = 42,
) -> Dict[str, Any]:
    """Permutational Multivariate Analysis of Variance (Anderson, 2001).

    PERMANOVA partitions the total sum of squared distances into
    within-group and between-group components and tests for significant
    differences among groups using a pseudo-F statistic.

    Args:
        distance_matrix: Symmetric pairwise distance matrix (N x N).
        group_labels: Group membership for each observation (length N).
        n_permutations: Number of permutations for the significance test.
        seed: Random seed for reproducibility.

    Returns:
        Dictionary with keys:
            - ``f_statistic``: Pseudo-F value.
            - ``p_value``: Proportion of permuted F values >= observed F.
            - ``r_squared``: Proportion of variance explained by grouping.
            - ``n_permutations``: Number of permutations performed.

    Raises:
        ValueError: If matrix is not square or labels do not match.

    Example:
        >>> dm = [[0, 1, 5, 5], [1, 0, 5, 5], [5, 5, 0, 1], [5, 5, 1, 0]]
        >>> labels = ["A", "A", "B", "B"]
        >>> result = permanova(dm, labels, n_permutations=99, seed=0)
        >>> result["f_statistic"] > 0
        True
    """
    n = _validate_square_matrix(distance_matrix)
    _validate_labels_match(group_labels, n)

    groups = _unique_groups(group_labels)
    if len(groups) < 2:
        raise ValueError("At least two distinct groups are required for PERMANOVA")

    rng = random.Random(seed)
    logger.info("Running PERMANOVA: %d observations, %d groups, %d permutations", n, len(groups), n_permutations)

    # Pre-compute squared distances
    sq_dist: List[List[float]] = [[d * d for d in row] for row in distance_matrix]

    def _compute_f(labels: List[str]) -> Tuple[float, float]:
        """Return (pseudo_F, r_squared) for given labels."""
        # Total sum of squared distances / N
        ss_total = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                ss_total += sq_dist[i][j]
        ss_total /= n  # Gower's centring

        # Within-group sum of squares
        ss_within = 0.0
        group_sizes: Dict[str, int] = defaultdict(int)
        for lb in labels:
            group_sizes[lb] += 1

        for grp in groups:
            members = [idx for idx, lb in enumerate(labels) if lb == grp]
            ng = len(members)
            if ng < 1:
                continue
            grp_ss = 0.0
            for a_idx in range(len(members)):
                for b_idx in range(a_idx + 1, len(members)):
                    grp_ss += sq_dist[members[a_idx]][members[b_idx]]
            ss_within += grp_ss / ng

        ss_between = ss_total - ss_within

        a = len(groups)  # number of groups
        df_between = a - 1
        df_within = n - a

        if df_within <= 0 or ss_within == 0:
            return 0.0, 0.0

        f_stat = (ss_between / df_between) / (ss_within / df_within)
        r_sq = ss_between / ss_total if ss_total > 0 else 0.0
        return f_stat, r_sq

    observed_f, observed_r2 = _compute_f(group_labels)

    n_geq = 0
    for _ in range(n_permutations):
        perm_labels = _permute_labels(group_labels, rng)
        perm_f, _ = _compute_f(perm_labels)
        if perm_f >= observed_f:
            n_geq += 1

    p_value = (n_geq + 1) / (n_permutations + 1)

    logger.info("PERMANOVA complete: F=%.4f, R2=%.4f, p=%.4f", observed_f, observed_r2, p_value)
    return {
        "f_statistic": round(observed_f, 6),
        "p_value": round(p_value, 4),
        "r_squared": round(observed_r2, 6),
        "n_permutations": n_permutations,
    }


# ---------------------------------------------------------------------------
# 4. Hierarchical clustering
# ---------------------------------------------------------------------------


def cluster_communities(
    distance_matrix: List[List[float]],
    method: str = "upgma",
    n_clusters: Optional[int] = None,
) -> Dict[str, Any]:
    """Hierarchical agglomerative clustering of communities.

    Pure-Python implementation supporting UPGMA (average), single-linkage,
    and complete-linkage methods. Optionally cuts the dendrogram to produce
    a specified number of clusters and computes the cophenetic correlation
    coefficient.

    Args:
        distance_matrix: Symmetric pairwise distance matrix (N x N).
        method: Linkage method -- ``"upgma"`` (average), ``"single"``, or
            ``"complete"``.
        n_clusters: If provided, cut the dendrogram to yield this many
            clusters. Must satisfy ``2 <= n_clusters <= N``.

    Returns:
        Dictionary with keys:
            - ``dendrogram``: List of merge steps, each a dict with
              ``cluster_a``, ``cluster_b``, ``distance``, and ``size``.
            - ``cluster_labels``: List of integer cluster assignments
              (only meaningful when *n_clusters* is given; otherwise every
              observation is its own cluster label equal to its index).
            - ``cophenetic_correlation``: Pearson correlation between
              original distances and cophenetic distances.

    Raises:
        ValueError: If *method* is unsupported or *n_clusters* is out of range.

    Example:
        >>> dm = [[0, 2, 6], [2, 0, 5], [6, 5, 0]]
        >>> result = cluster_communities(dm, method="upgma", n_clusters=2)
        >>> result["cluster_labels"]
        [0, 0, 1]
    """
    supported = ("upgma", "single", "complete")
    if method not in supported:
        raise ValueError(f"Unsupported linkage method '{method}'. Choose from {supported}")

    n = _validate_square_matrix(distance_matrix)

    if n_clusters is not None:
        if n_clusters < 1 or n_clusters > n:
            raise ValueError(f"n_clusters must be between 1 and {n}, got {n_clusters}")

    logger.info("Clustering %d observations with %s linkage", n, method)

    # Working copy of distance matrix (dict-based for dynamic cluster merging)
    # Each active cluster is keyed by an integer id. Initially ids 0..n-1.
    active: Dict[int, List[int]] = {i: [i] for i in range(n)}  # cluster_id -> member indices
    dist: Dict[Tuple[int, int], float] = {}
    for i in range(n):
        for j in range(i + 1, n):
            dist[(i, j)] = distance_matrix[i][j]

    def _get_dist(a: int, b: int) -> float:
        key = (min(a, b), max(a, b))
        return dist.get(key, 0.0)

    dendrogram: List[Dict[str, Any]] = []
    next_id = n  # id counter for new clusters

    for _step in range(n - 1):
        # Find closest pair among active clusters
        best_pair: Optional[Tuple[int, int]] = None
        best_dist = math.inf
        ids = sorted(active.keys())
        for idx_a in range(len(ids)):
            for idx_b in range(idx_a + 1, len(ids)):
                d = _get_dist(ids[idx_a], ids[idx_b])
                if d < best_dist:
                    best_dist = d
                    best_pair = (ids[idx_a], ids[idx_b])

        if best_pair is None:
            break  # pragma: no cover

        ca, cb = best_pair
        members_a = active[ca]
        members_b = active[cb]
        merged_members = members_a + members_b

        dendrogram.append(
            {
                "cluster_a": ca,
                "cluster_b": cb,
                "distance": round(best_dist, 8),
                "size": len(merged_members),
            }
        )

        # Compute distances from new cluster to all remaining active clusters
        new_id = next_id
        next_id += 1

        remaining = [cid for cid in active if cid != ca and cid != cb]
        for cid in remaining:
            d_a = _get_dist(ca, cid)
            d_b = _get_dist(cb, cid)

            if method == "single":
                new_d = min(d_a, d_b)
            elif method == "complete":
                new_d = max(d_a, d_b)
            else:  # upgma
                na = len(members_a)
                nb = len(members_b)
                new_d = (d_a * na + d_b * nb) / (na + nb)

            key = (min(new_id, cid), max(new_id, cid))
            dist[key] = new_d

        # Update active clusters
        del active[ca]
        del active[cb]
        active[new_id] = merged_members

    # --- Cluster labels via dendrogram cutting ---
    if n_clusters is None or n_clusters >= n:
        cluster_labels = list(range(n))
    elif n_clusters == 1:
        cluster_labels = [0] * n
    else:
        # Replay merges and stop when we reach the desired number of clusters
        label_map: Dict[int, int] = {i: i for i in range(n)}
        current_n_clusters = n
        merge_idx = 0
        while current_n_clusters > n_clusters and merge_idx < len(dendrogram):
            step = dendrogram[merge_idx]
            ca_id = step["cluster_a"]
            cb_id = step["cluster_b"]
            # Merge cb into ca (relabel all cb members to ca's label)
            target_label = label_map.get(ca_id, ca_id)
            source_label = label_map.get(cb_id, cb_id)
            # Relabel all observations that have source_label
            for obs in range(n):
                if label_map[obs] == source_label:
                    label_map[obs] = target_label
            # Also update the new merged cluster id to use the target label
            new_merged_id = n + merge_idx
            label_map[new_merged_id] = target_label
            current_n_clusters -= 1
            merge_idx += 1

        # Remap to contiguous 0-based labels
        obs_labels = [label_map[i] for i in range(n)]
        unique_labels = sorted(set(obs_labels))
        remap = {old: new for new, old in enumerate(unique_labels)}
        cluster_labels = [remap[lb] for lb in obs_labels]

    # --- Cophenetic correlation ---
    cophenetic_dist = _cophenetic_distances(n, dendrogram)
    coph_corr = _pearson_correlation(distance_matrix, cophenetic_dist, n)

    logger.info("Clustering complete. Cophenetic correlation: %.4f", coph_corr)
    return {
        "dendrogram": dendrogram,
        "cluster_labels": cluster_labels,
        "cophenetic_correlation": round(coph_corr, 6),
    }


def _cophenetic_distances(n: int, dendrogram: List[Dict[str, Any]]) -> List[List[float]]:
    """Build cophenetic distance matrix from dendrogram merge history.

    The cophenetic distance between two original observations is the
    distance at which they first become members of the same cluster.

    Args:
        n: Number of original observations.
        dendrogram: Merge history from :func:`cluster_communities`.

    Returns:
        N x N cophenetic distance matrix.
    """
    coph = [[0.0] * n for _ in range(n)]

    # Track cluster membership: cluster_id -> set of original observations
    members: Dict[int, set] = {i: {i} for i in range(n)}

    for step_idx, step in enumerate(dendrogram):
        ca = step["cluster_a"]
        cb = step["cluster_b"]
        d = step["distance"]
        new_id = n + step_idx

        ma = members.get(ca, set())
        mb = members.get(cb, set())

        # Set cophenetic distance for all cross-cluster pairs
        for obs_a in ma:
            for obs_b in mb:
                coph[obs_a][obs_b] = d
                coph[obs_b][obs_a] = d

        members[new_id] = ma | mb
        # Clean up (optional; keeps memory lower)
        members.pop(ca, None)
        members.pop(cb, None)

    return coph


def _pearson_correlation(mat_a: List[List[float]], mat_b: List[List[float]], n: int) -> float:
    """Pearson correlation of the upper-triangle elements of two N x N matrices.

    Args:
        mat_a: First symmetric matrix.
        mat_b: Second symmetric matrix.
        n: Matrix dimension.

    Returns:
        Pearson r in [-1, 1], or 0.0 if undefined.
    """
    vals_a: List[float] = []
    vals_b: List[float] = []
    for i in range(n):
        for j in range(i + 1, n):
            vals_a.append(mat_a[i][j])
            vals_b.append(mat_b[i][j])

    if len(vals_a) < 2:
        return 0.0

    mean_a = sum(vals_a) / len(vals_a)
    mean_b = sum(vals_b) / len(vals_b)

    cov = sum((a - mean_a) * (b - mean_b) for a, b in zip(vals_a, vals_b))
    var_a = sum((a - mean_a) ** 2 for a in vals_a)
    var_b = sum((b - mean_b) ** 2 for b in vals_b)

    denom = math.sqrt(var_a * var_b)
    return cov / denom if denom > 0 else 0.0


# ---------------------------------------------------------------------------
# 5. SIMPER
# ---------------------------------------------------------------------------


def simper(
    abundance_matrix: List[List[float]],
    group_labels: List[str],
) -> List[Dict[str, Any]]:
    """Similarity Percentages (SIMPER) analysis (Clarke, 1993).

    Decomposes Bray-Curtis dissimilarity between groups to identify which
    species contribute most to observed differences. For every pair of
    groups, species contributions are calculated and the results across all
    group pairs are aggregated.

    Args:
        abundance_matrix: Sites-by-species matrix.
        group_labels: Group membership for each site.

    Returns:
        A list of dicts sorted by descending ``contribution_pct``, each with:
            - ``species_idx``: Column index of the species.
            - ``avg_dissimilarity``: Mean per-species Bray-Curtis contribution.
            - ``contribution_pct``: Percentage of total average dissimilarity.
            - ``cumulative_pct``: Running cumulative percentage.
            - ``sd``: Standard deviation of per-pair contributions.
            - ``avg_abundance_by_group``: Dict mapping group label to mean
              abundance of the species in that group.

    Raises:
        ValueError: If inputs are empty or inconsistent.

    Example:
        >>> mat = [[10, 0], [8, 2], [0, 9], [1, 7]]
        >>> labels = ["A", "A", "B", "B"]
        >>> res = simper(mat, labels)
        >>> res[0]["species_idx"] in (0, 1)
        True
    """
    validation.validate_not_empty(abundance_matrix, "abundance_matrix")
    n_sites = len(abundance_matrix)
    _validate_labels_match(group_labels, n_sites)
    n_species = len(abundance_matrix[0])

    groups = _unique_groups(group_labels)
    if len(groups) < 2:
        raise ValueError("At least two distinct groups are required for SIMPER")

    logger.info("Running SIMPER: %d sites, %d species, %d groups", n_sites, n_species, len(groups))

    # Group indices
    group_idx: Dict[str, List[int]] = defaultdict(list)
    for idx, grp in enumerate(group_labels):
        group_idx[grp].append(idx)

    # Accumulate species-level contributions across all group pairs
    species_contribs: Dict[int, List[float]] = defaultdict(list)

    for g_idx_a in range(len(groups)):
        for g_idx_b in range(g_idx_a + 1, len(groups)):
            ga, gb = groups[g_idx_a], groups[g_idx_b]
            sites_a = group_idx[ga]
            sites_b = group_idx[gb]

            # For each between-group site pair, decompose Bray-Curtis per species
            for sa in sites_a:
                for sb in sites_b:
                    row_a = abundance_matrix[sa]
                    row_b = abundance_matrix[sb]
                    total_sum = sum(row_a[k] + row_b[k] for k in range(n_species))
                    if total_sum == 0:
                        continue
                    for sp in range(n_species):
                        contrib = abs(row_a[sp] - row_b[sp]) / total_sum
                        species_contribs[sp].append(contrib)

    # Summarise
    results_raw: List[Dict[str, Any]] = []
    for sp in range(n_species):
        contribs = species_contribs.get(sp, [0.0])
        avg_d = sum(contribs) / len(contribs) if contribs else 0.0

        # Standard deviation
        if len(contribs) > 1:
            mean_c = avg_d
            var_c = sum((c - mean_c) ** 2 for c in contribs) / (len(contribs) - 1)
            sd_c = math.sqrt(var_c)
        else:
            sd_c = 0.0

        # Mean abundance per group
        avg_by_group: Dict[str, float] = {}
        for grp in groups:
            vals = [abundance_matrix[s][sp] for s in group_idx[grp]]
            avg_by_group[grp] = sum(vals) / len(vals) if vals else 0.0

        results_raw.append(
            {
                "species_idx": sp,
                "avg_dissimilarity": avg_d,
                "sd": sd_c,
                "avg_abundance_by_group": avg_by_group,
            }
        )

    # Sort by contribution (descending) and compute percentages
    total_dissimilarity = sum(r["avg_dissimilarity"] for r in results_raw)
    results_raw.sort(key=lambda r: r["avg_dissimilarity"], reverse=True)

    cumulative = 0.0
    for r in results_raw:
        pct = (r["avg_dissimilarity"] / total_dissimilarity * 100.0) if total_dissimilarity > 0 else 0.0
        cumulative += pct
        r["contribution_pct"] = round(pct, 4)
        r["cumulative_pct"] = round(cumulative, 4)
        r["avg_dissimilarity"] = round(r["avg_dissimilarity"], 6)
        r["sd"] = round(r["sd"], 6)

    logger.info("SIMPER complete: total avg dissimilarity = %.4f", total_dissimilarity)
    return results_raw


# ---------------------------------------------------------------------------
# 6. PERMDISP (Multivariate Dispersion)
# ---------------------------------------------------------------------------


def multivariate_dispersion(
    distance_matrix: List[List[float]],
    group_labels: List[str],
    *,
    n_permutations: int = 999,
    seed: int = 42,
) -> Dict[str, Any]:
    """Test for homogeneity of multivariate dispersions (PERMDISP / betadisper).

    Computes the distance from each observation to its group centroid (using
    the distance matrix directly via the principal-coordinate-free centroid
    approach of Anderson 2006), then tests whether group dispersions differ
    using a permutation-based F-test on the centroid distances.

    Args:
        distance_matrix: Symmetric pairwise distance matrix (N x N).
        group_labels: Group membership for each observation (length N).
        n_permutations: Number of permutations for the significance test.
        seed: Random seed for reproducibility.

    Returns:
        Dictionary with keys:
            - ``f_statistic``: F-value comparing group dispersions.
            - ``p_value``: Permutation *p*-value.
            - ``group_dispersions``: Dict mapping each group label to its
              mean distance-to-centroid.
            - ``n_permutations``: Number of permutations performed.

    Raises:
        ValueError: If matrix is not square or labels do not match.

    Example:
        >>> dm = [[0, 1, 4, 5], [1, 0, 5, 4], [4, 5, 0, 1], [5, 4, 1, 0]]
        >>> labels = ["X", "X", "Y", "Y"]
        >>> result = multivariate_dispersion(dm, labels, n_permutations=99, seed=0)
        >>> "f_statistic" in result
        True
    """
    n = _validate_square_matrix(distance_matrix)
    _validate_labels_match(group_labels, n)

    groups = _unique_groups(group_labels)
    if len(groups) < 2:
        raise ValueError("At least two distinct groups are required for PERMDISP")

    rng = random.Random(seed)
    logger.info("Running PERMDISP: %d observations, %d groups", n, len(groups))

    # Pre-compute squared distances
    sq_dist: List[List[float]] = [[d * d for d in row] for row in distance_matrix]

    def _distances_to_centroids(labels: List[str]) -> List[float]:
        """Compute distance from each observation to its group centroid.

        Uses the identity: d(i, centroid_g)^2 = (1/n_g) * sum_j d(i,j)^2
        - (1/n_g^2) * sum_{j<k} d(j,k)^2  where j,k in group g.

        This avoids explicit coordinate embedding.
        """
        g_members: Dict[str, List[int]] = defaultdict(list)
        for idx, lb in enumerate(labels):
            g_members[lb].append(idx)

        # Pre-compute within-group sum of squared distances for each group
        g_within_ss: Dict[str, float] = {}
        for grp, members in g_members.items():
            ss = 0.0
            for a_idx in range(len(members)):
                for b_idx in range(a_idx + 1, len(members)):
                    ss += sq_dist[members[a_idx]][members[b_idx]]
            g_within_ss[grp] = ss

        dists = [0.0] * n
        for obs in range(n):
            grp = labels[obs]
            members = g_members[grp]
            ng = len(members)
            if ng <= 1:
                dists[obs] = 0.0
                continue
            # sum of squared distances from obs to all other members
            sum_sq = sum(sq_dist[obs][m] for m in members)
            d_sq = sum_sq / ng - g_within_ss[grp] / (ng * ng)
            dists[obs] = math.sqrt(max(d_sq, 0.0))

        return dists

    def _compute_f(centroid_dists: List[float], labels: List[str]) -> Tuple[float, Dict[str, float]]:
        """ANOVA F-statistic on distances-to-centroid across groups."""
        g_members: Dict[str, List[int]] = defaultdict(list)
        for idx, lb in enumerate(labels):
            g_members[lb].append(idx)

        # Group means
        group_means: Dict[str, float] = {}
        for grp, members in g_members.items():
            vals = [centroid_dists[m] for m in members]
            group_means[grp] = sum(vals) / len(vals) if vals else 0.0

        grand_mean = sum(centroid_dists) / n if n > 0 else 0.0

        # SS between
        ss_between = sum(
            len(g_members[grp]) * (group_means[grp] - grand_mean) ** 2 for grp in groups if grp in g_members
        )
        # SS within
        ss_within = 0.0
        for grp, members in g_members.items():
            gm = group_means[grp]
            ss_within += sum((centroid_dists[m] - gm) ** 2 for m in members)

        a = len([g for g in groups if g in g_members and len(g_members[g]) > 0])
        df_between = a - 1
        df_within = n - a

        if df_within <= 0 or ss_within == 0:
            return 0.0, group_means

        f_val = (ss_between / df_between) / (ss_within / df_within) if df_between > 0 else 0.0
        return f_val, group_means

    observed_dists = _distances_to_centroids(group_labels)
    observed_f, group_dispersions = _compute_f(observed_dists, group_labels)

    # Round group dispersions for output
    group_dispersions_out = {g: round(v, 6) for g, v in group_dispersions.items()}

    # Permutation test: permute group labels, recompute centroid distances and F
    n_geq = 0
    for _ in range(n_permutations):
        perm_labels = _permute_labels(group_labels, rng)
        perm_dists = _distances_to_centroids(perm_labels)
        perm_f, _ = _compute_f(perm_dists, perm_labels)
        if perm_f >= observed_f:
            n_geq += 1

    p_value = (n_geq + 1) / (n_permutations + 1)

    logger.info("PERMDISP complete: F=%.4f, p=%.4f", observed_f, p_value)
    return {
        "f_statistic": round(observed_f, 6),
        "p_value": round(p_value, 4),
        "group_dispersions": group_dispersions_out,
        "n_permutations": n_permutations,
    }
