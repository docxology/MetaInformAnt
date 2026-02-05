"""Functional ecology analysis: trait-based diversity and community metrics.

This module provides comprehensive functional ecology methods for analyzing
community structure through species traits. Implements Functional Richness (FRic),
Functional Evenness (FEve), Functional Divergence (FDiv), Functional Dispersion
(FDis), Rao's Quadratic Entropy, Community Weighted Means, and related metrics.

All algorithms are pure Python implementations with no external dependencies
beyond the standard library.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

from metainformant.core import logging, errors, validation

logger = logging.get_logger(__name__)


# ---------------------------------------------------------------------------
# Helper: distance functions
# ---------------------------------------------------------------------------


def _euclidean_distance(a: List[float], b: List[float]) -> float:
    """Compute Euclidean distance between two trait vectors.

    Args:
        a: First trait vector.
        b: Second trait vector (same length as *a*).

    Returns:
        Euclidean distance.

    Raises:
        ValueError: If vectors differ in length.
    """
    if len(a) != len(b):
        raise ValueError(f"Vectors must have equal length, got {len(a)} and {len(b)}")
    return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))


def _gower_distance(a: List[float], b: List[float], ranges: List[float]) -> float:
    """Compute Gower distance between two trait vectors.

    Gower distance normalises each trait by its range so that mixed-scale traits
    contribute equally.  ``d = (1/T) * sum(|a_k - b_k| / range_k)`` where *T* is
    the number of traits with non-zero range.

    Args:
        a: First trait vector.
        b: Second trait vector.
        ranges: Pre-computed range (max - min) for each trait across the full
            trait matrix.

    Returns:
        Gower distance in [0, 1].

    Raises:
        ValueError: If vector lengths differ.
    """
    if len(a) != len(b) or len(a) != len(ranges):
        raise ValueError("Vectors and ranges must have equal length")
    total = 0.0
    valid_traits = 0
    for ai, bi, r in zip(a, b, ranges):
        if r > 0:
            total += abs(ai - bi) / r
            valid_traits += 1
    return total / valid_traits if valid_traits > 0 else 0.0


# ---------------------------------------------------------------------------
# Helper: trait distance matrix
# ---------------------------------------------------------------------------


def trait_distance_matrix(
    trait_matrix: List[List[float]], method: str = "euclidean"
) -> List[List[float]]:
    """Compute pairwise trait distance matrix for all species.

    Args:
        trait_matrix: Species-by-traits matrix where ``trait_matrix[i]`` is the
            trait vector for species *i*.
        method: Distance method -- ``"euclidean"`` or ``"gower"``.

    Returns:
        Symmetric ``n x n`` distance matrix.

    Raises:
        ValueError: If *method* is unsupported or the matrix is empty.

    Example:
        >>> dm = trait_distance_matrix([[0, 0], [3, 4]], method="euclidean")
        >>> round(dm[0][1], 1)
        5.0
    """
    supported = ("euclidean", "gower")
    if method not in supported:
        raise ValueError(f"Unsupported distance method: {method}. Choose from {list(supported)}")

    validation.validate_not_empty(trait_matrix, "trait_matrix")
    n = len(trait_matrix)

    # Pre-compute ranges for Gower
    ranges: List[float] = []
    if method == "gower":
        n_traits = len(trait_matrix[0])
        for t in range(n_traits):
            vals = [trait_matrix[i][t] for i in range(n)]
            ranges.append(max(vals) - min(vals))

    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if method == "euclidean":
                d = _euclidean_distance(trait_matrix[i], trait_matrix[j])
            else:
                d = _gower_distance(trait_matrix[i], trait_matrix[j], ranges)
            dist[i][j] = d
            dist[j][i] = d
    return dist


# ---------------------------------------------------------------------------
# Helper: 2-D convex hull (Graham scan)
# ---------------------------------------------------------------------------


def _cross(o: Tuple[float, float], a: Tuple[float, float], b: Tuple[float, float]) -> float:
    """Return the cross product of vectors OA and OB (positive if counter-clockwise)."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def convex_hull_2d(points: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """Compute the 2-D convex hull using Andrew's monotone-chain algorithm.

    Args:
        points: List of (x, y) tuples.

    Returns:
        Vertices of the convex hull in counter-clockwise order (closed polygon;
        first vertex is *not* repeated at the end).

    Example:
        >>> hull = convex_hull_2d([(0, 0), (1, 0), (0.5, 0.5), (0, 1), (1, 1)])
        >>> len(hull)
        4
    """
    pts = sorted(set(points))
    if len(pts) <= 1:
        return list(pts)

    # Build lower hull
    lower: List[Tuple[float, float]] = []
    for p in pts:
        while len(lower) >= 2 and _cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper: List[Tuple[float, float]] = []
    for p in reversed(pts):
        while len(upper) >= 2 and _cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Remove last point of each half because it repeats at the junction
    return lower[:-1] + upper[:-1]


def _polygon_area(vertices: List[Tuple[float, float]]) -> float:
    """Compute polygon area using the Shoelace formula.

    Args:
        vertices: Ordered vertices of the polygon.

    Returns:
        Absolute area of the polygon.
    """
    n = len(vertices)
    if n < 3:
        return 0.0
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += vertices[i][0] * vertices[j][1]
        area -= vertices[j][0] * vertices[i][1]
    return abs(area) / 2.0


# ---------------------------------------------------------------------------
# Helper: minimum spanning tree (Prim's algorithm)
# ---------------------------------------------------------------------------


def minimum_spanning_tree(dist_matrix: List[List[float]]) -> List[Tuple[int, int, float]]:
    """Build a minimum spanning tree using Prim's algorithm.

    Args:
        dist_matrix: Symmetric ``n x n`` distance matrix.

    Returns:
        List of ``(i, j, weight)`` edges in the MST.

    Example:
        >>> edges = minimum_spanning_tree([[0, 1, 3], [1, 0, 2], [3, 2, 0]])
        >>> sum(w for _, _, w in edges)
        3
    """
    n = len(dist_matrix)
    if n <= 1:
        return []

    in_tree = [False] * n
    min_cost = [float("inf")] * n
    min_edge = [-1] * n
    edges: List[Tuple[int, int, float]] = []

    # Start from vertex 0
    min_cost[0] = 0.0
    for _ in range(n):
        # Find the closest vertex not yet in the tree
        u = -1
        for v in range(n):
            if not in_tree[v] and (u == -1 or min_cost[v] < min_cost[u]):
                u = v
        in_tree[u] = True

        if min_edge[u] != -1:
            edges.append((min_edge[u], u, dist_matrix[min_edge[u]][u]))

        # Update candidate edges
        for v in range(n):
            if not in_tree[v] and dist_matrix[u][v] < min_cost[v]:
                min_cost[v] = dist_matrix[u][v]
                min_edge[v] = u

    return edges


# ---------------------------------------------------------------------------
# Helper: relative abundances
# ---------------------------------------------------------------------------


def _to_relative(abundances: List[float]) -> List[float]:
    """Convert raw abundances to relative proportions that sum to 1.

    Args:
        abundances: Raw abundance values (non-negative).

    Returns:
        List of proportions.
    """
    total = sum(abundances)
    if total == 0:
        return [0.0] * len(abundances)
    return [a / total for a in abundances]


# ---------------------------------------------------------------------------
# Public API: Functional diversity metrics
# ---------------------------------------------------------------------------


def functional_richness(
    trait_matrix: List[List[float]], abundances: Optional[List[float]] = None
) -> float:
    """Calculate Functional Richness (FRic).

    FRic measures the volume of functional trait space occupied by the community.
    For 1-D it is the range; for 2-D it is the convex hull area; for 3-D and
    above it is approximated from the maximum pairwise trait distance envelope.

    Args:
        trait_matrix: Species-by-traits matrix (each row is one species).
        abundances: Optional abundance vector.  If provided, species with zero
            abundance are excluded before computing FRic.

    Returns:
        Functional richness value (>= 0).

    Raises:
        ValueError: If the trait matrix is empty.

    Example:
        >>> functional_richness([[0], [5], [10]])
        10.0
        >>> round(functional_richness([[0, 0], [1, 0], [0, 1]]), 1)
        0.5
    """
    validation.validate_not_empty(trait_matrix, "trait_matrix")
    logger.info("Calculating Functional Richness (FRic)")

    # Filter by abundance if provided
    if abundances is not None:
        if len(abundances) != len(trait_matrix):
            raise ValueError("abundances length must match number of species in trait_matrix")
        filtered = [trait_matrix[i] for i in range(len(trait_matrix)) if abundances[i] > 0]
    else:
        filtered = list(trait_matrix)

    if not filtered:
        return 0.0

    n_species = len(filtered)
    n_traits = len(filtered[0])

    # Single species -- zero volume
    if n_species <= 1:
        return 0.0

    # 1-D: range of trait values
    if n_traits == 1:
        vals = [row[0] for row in filtered]
        return max(vals) - min(vals)

    # 2-D: convex hull area
    if n_traits == 2:
        points = [(row[0], row[1]) for row in filtered]
        hull = convex_hull_2d(points)
        return _polygon_area(hull)

    # 3-D+: approximate as product of per-trait ranges (hyper-volume approximation)
    # This gives a bounding-box volume which is a standard FRic approximation
    # when convex hull algorithms for > 2-D are unavailable.
    volume = 1.0
    any_nonzero = False
    for t in range(n_traits):
        vals = [filtered[i][t] for i in range(n_species)]
        r = max(vals) - min(vals)
        if r > 0:
            volume *= r
            any_nonzero = True
    return volume if any_nonzero else 0.0


def functional_evenness(trait_matrix: List[List[float]], abundances: List[float]) -> float:
    """Calculate Functional Evenness (FEve).

    FEve quantifies the regularity of the abundance distribution along the
    minimum spanning tree in trait space.  Values range from 0 (uneven) to 1
    (perfectly even).

    The algorithm (Villeger et al. 2008):
        1. Build the MST over species in trait space.
        2. For each MST branch (i, j) compute the *partial weighted evenness*:
           ``EW_b = dist_ij / (abundance_i + abundance_j)``.
        3. Normalise to ``PEW_b = EW_b / sum(EW)``.
        4. ``FEve = (sum(min(PEW_b, 1/(S-1))) - 1/(S-1)) / (1 - 1/(S-1))``.

    Args:
        trait_matrix: Species-by-traits matrix.
        abundances: Species abundance vector (must have same length as
            trait_matrix rows).

    Returns:
        FEve value in [0, 1].

    Raises:
        ValueError: If inputs are invalid.

    Example:
        >>> feve = functional_evenness([[0], [5], [10]], [10, 10, 10])
        >>> round(feve, 2)
        1.0
    """
    validation.validate_not_empty(trait_matrix, "trait_matrix")
    validation.validate_not_empty(abundances, "abundances")
    if len(abundances) != len(trait_matrix):
        raise ValueError("abundances length must match number of species in trait_matrix")
    logger.info("Calculating Functional Evenness (FEve)")

    # Filter to species with positive abundance
    idx = [i for i in range(len(abundances)) if abundances[i] > 0]
    if len(idx) < 2:
        return 0.0

    tm = [trait_matrix[i] for i in idx]
    ab = [abundances[i] for i in idx]
    s = len(tm)

    # Distance matrix among present species
    dm = trait_distance_matrix(tm, method="euclidean")

    # MST
    mst_edges = minimum_spanning_tree(dm)
    if not mst_edges:
        return 0.0

    # Weighted evenness per branch
    ew_values: List[float] = []
    for i, j, d in mst_edges:
        denom = ab[i] + ab[j]
        ew = d / denom if denom > 0 else 0.0
        ew_values.append(ew)

    total_ew = sum(ew_values)
    if total_ew == 0:
        return 0.0

    # Partial weighted evenness
    pew = [ew / total_ew for ew in ew_values]

    threshold = 1.0 / (s - 1)
    numerator = sum(min(p, threshold) for p in pew) - threshold
    denominator = 1.0 - threshold

    if denominator <= 0:
        return 0.0

    feve = numerator / denominator
    # Clamp to [0, 1]
    return max(0.0, min(1.0, feve))


def functional_divergence(trait_matrix: List[List[float]], abundances: List[float]) -> float:
    """Calculate Functional Divergence (FDiv).

    FDiv measures how species abundances are distributed relative to the centre
    of the functional trait space.  High FDiv indicates that the most abundant
    species have extreme trait values.

    The algorithm (Villeger et al. 2008):
        1. Compute the abundance-weighted centroid *G*.
        2. Compute each species' Euclidean distance to the centroid: ``dG_i``.
        3. Compute the mean distance: ``mean_dG = sum(dG_i) / S``.
        4. ``FDiv = (deltaD + mean_dG) / (deltaAbs + mean_dG)``
           where ``deltaD  = sum(w_i * (dG_i - mean_dG))``
           and   ``deltaAbs = sum(w_i * |dG_i - mean_dG|)``.

    Args:
        trait_matrix: Species-by-traits matrix.
        abundances: Species abundance vector.

    Returns:
        FDiv value in [0, 1].

    Raises:
        ValueError: If inputs are invalid.

    Example:
        >>> round(functional_divergence([[0], [5], [10]], [1, 1, 1]), 1)
        0.6
    """
    validation.validate_not_empty(trait_matrix, "trait_matrix")
    validation.validate_not_empty(abundances, "abundances")
    if len(abundances) != len(trait_matrix):
        raise ValueError("abundances length must match number of species in trait_matrix")
    logger.info("Calculating Functional Divergence (FDiv)")

    # Filter to positive abundances
    idx = [i for i in range(len(abundances)) if abundances[i] > 0]
    if len(idx) < 2:
        return 0.0

    tm = [trait_matrix[i] for i in idx]
    ab = [abundances[i] for i in idx]
    rel = _to_relative(ab)
    s = len(tm)
    n_traits = len(tm[0])

    # Abundance-weighted centroid
    centroid = [0.0] * n_traits
    for t in range(n_traits):
        centroid[t] = sum(rel[i] * tm[i][t] for i in range(s))

    # Distances to centroid
    dg = [_euclidean_distance(tm[i], centroid) for i in range(s)]
    mean_dg = sum(dg) / s if s > 0 else 0.0

    # Weighted deviations
    delta_d = sum(rel[i] * (dg[i] - mean_dg) for i in range(s))
    delta_abs = sum(rel[i] * abs(dg[i] - mean_dg) for i in range(s))

    denominator = delta_abs + mean_dg
    if denominator == 0:
        return 0.0

    fdiv = (delta_d + mean_dg) / denominator
    return max(0.0, min(1.0, fdiv))


def functional_dispersion(trait_matrix: List[List[float]], abundances: List[float]) -> float:
    """Calculate Functional Dispersion (FDis).

    FDis is the abundance-weighted mean distance of species to the weighted
    centroid of the community in trait space (Laliberte & Legendre 2010).

    Args:
        trait_matrix: Species-by-traits matrix.
        abundances: Species abundance vector.

    Returns:
        FDis value (>= 0).

    Raises:
        ValueError: If inputs are invalid.

    Example:
        >>> round(functional_dispersion([[0], [10]], [1, 1]), 1)
        5.0
    """
    validation.validate_not_empty(trait_matrix, "trait_matrix")
    validation.validate_not_empty(abundances, "abundances")
    if len(abundances) != len(trait_matrix):
        raise ValueError("abundances length must match number of species in trait_matrix")
    logger.info("Calculating Functional Dispersion (FDis)")

    # Filter to positive abundances
    idx = [i for i in range(len(abundances)) if abundances[i] > 0]
    if len(idx) < 2:
        return 0.0

    tm = [trait_matrix[i] for i in idx]
    ab = [abundances[i] for i in idx]
    rel = _to_relative(ab)
    s = len(tm)
    n_traits = len(tm[0])

    # Weighted centroid
    centroid = [0.0] * n_traits
    for t in range(n_traits):
        centroid[t] = sum(rel[i] * tm[i][t] for i in range(s))

    # Weighted mean distance to centroid
    fdis = sum(rel[i] * _euclidean_distance(tm[i], centroid) for i in range(s))
    return fdis


def raos_quadratic_entropy(trait_matrix: List[List[float]], abundances: List[float]) -> float:
    """Calculate Rao's Quadratic Entropy (Q).

    ``Q = sum_i sum_j (d_ij * p_i * p_j)`` where ``d_ij`` is the Euclidean
    trait dissimilarity and ``p_i`` is the relative abundance of species *i*.

    Args:
        trait_matrix: Species-by-traits matrix.
        abundances: Species abundance vector.

    Returns:
        Rao's Q value (>= 0).

    Raises:
        ValueError: If inputs are invalid.

    Example:
        >>> round(raos_quadratic_entropy([[0], [10]], [1, 1]), 1)
        5.0
    """
    validation.validate_not_empty(trait_matrix, "trait_matrix")
    validation.validate_not_empty(abundances, "abundances")
    if len(abundances) != len(trait_matrix):
        raise ValueError("abundances length must match number of species in trait_matrix")
    logger.info("Calculating Rao's Quadratic Entropy")

    # Filter to positive abundances
    idx = [i for i in range(len(abundances)) if abundances[i] > 0]
    if len(idx) < 2:
        return 0.0

    tm = [trait_matrix[i] for i in idx]
    ab = [abundances[i] for i in idx]
    rel = _to_relative(ab)
    s = len(tm)

    # Build distance matrix and compute Q
    q = 0.0
    for i in range(s):
        for j in range(i + 1, s):
            d = _euclidean_distance(tm[i], tm[j])
            q += d * rel[i] * rel[j] * 2  # symmetric: count (i,j) and (j,i)
    return q


def community_weighted_mean(trait_matrix: List[List[float]], abundances: List[float]) -> List[float]:
    """Calculate Community Weighted Mean (CWM) trait values.

    CWM_t = sum(p_i * trait_it) for each trait *t*, where ``p_i`` is the
    relative abundance of species *i*.

    Args:
        trait_matrix: Species-by-traits matrix.
        abundances: Species abundance vector.

    Returns:
        List of CWM values (one per trait).

    Raises:
        ValueError: If inputs are invalid.

    Example:
        >>> community_weighted_mean([[2, 4], [8, 6]], [1, 1])
        [5.0, 5.0]
    """
    validation.validate_not_empty(trait_matrix, "trait_matrix")
    validation.validate_not_empty(abundances, "abundances")
    if len(abundances) != len(trait_matrix):
        raise ValueError("abundances length must match number of species in trait_matrix")
    logger.info("Calculating Community Weighted Mean (CWM)")

    # Filter to positive abundances
    idx = [i for i in range(len(abundances)) if abundances[i] > 0]
    if not idx:
        return [0.0] * len(trait_matrix[0])

    tm = [trait_matrix[i] for i in idx]
    ab = [abundances[i] for i in idx]
    rel = _to_relative(ab)
    n_traits = len(tm[0])

    cwm = [0.0] * n_traits
    for t in range(n_traits):
        cwm[t] = sum(rel[i] * tm[i][t] for i in range(len(tm)))
    return cwm


def functional_redundancy(trait_matrix: List[List[float]], abundances: List[float]) -> float:
    """Calculate Functional Redundancy.

    Functional redundancy is defined as the difference between Simpson diversity
    and Rao's Quadratic Entropy:

        ``FR = Simpson - Rao's Q``

    High functional redundancy indicates that taxonomic diversity exceeds
    functional diversity, meaning many species share similar trait values.

    Args:
        trait_matrix: Species-by-traits matrix.
        abundances: Species abundance vector.

    Returns:
        Functional redundancy value.  Typically >= 0 when trait distances are
        normalised to [0, 1], but can be negative for un-normalised distances.

    Raises:
        ValueError: If inputs are invalid.

    Example:
        >>> fr = functional_redundancy([[0], [0]], [5, 5])
        >>> round(fr, 2)
        0.5
    """
    validation.validate_not_empty(trait_matrix, "trait_matrix")
    validation.validate_not_empty(abundances, "abundances")
    if len(abundances) != len(trait_matrix):
        raise ValueError("abundances length must match number of species in trait_matrix")
    logger.info("Calculating Functional Redundancy")

    # Simpson diversity: 1 - sum(p_i^2)
    pos_ab = [a for a in abundances if a > 0]
    if len(pos_ab) < 2:
        return 0.0

    rel = _to_relative(pos_ab)
    simpson = 1.0 - sum(p ** 2 for p in rel)

    rao_q = raos_quadratic_entropy(trait_matrix, abundances)

    return simpson - rao_q


def functional_beta_diversity(
    community1_traits: List[List[float]],
    community1_abundances: List[float],
    community2_traits: List[List[float]],
    community2_abundances: List[float],
) -> Dict[str, float]:
    """Calculate Functional Beta Diversity between two communities.

    Decomposes total functional dissimilarity into turnover and nestedness
    components following an approach inspired by Villeger et al. (2013):

    * **total**: Rao's Q of the pooled community minus the mean of the
      individual Rao's Q values.
    * **turnover**: Fraction attributable to trait replacement between
      communities.
    * **nestedness**: Fraction attributable to trait loss (one community being a
      functional subset of the other).

    Args:
        community1_traits: Trait matrix for community 1.
        community1_abundances: Abundance vector for community 1.
        community2_traits: Trait matrix for community 2.
        community2_abundances: Abundance vector for community 2.

    Returns:
        Dictionary with keys ``total``, ``turnover``, and ``nestedness``.

    Raises:
        ValueError: If inputs are invalid.

    Example:
        >>> result = functional_beta_diversity(
        ...     [[0], [5]], [1, 1],
        ...     [[10], [15]], [1, 1],
        ... )
        >>> result["total"] > 0
        True
    """
    validation.validate_not_empty(community1_traits, "community1_traits")
    validation.validate_not_empty(community1_abundances, "community1_abundances")
    validation.validate_not_empty(community2_traits, "community2_traits")
    validation.validate_not_empty(community2_abundances, "community2_abundances")
    logger.info("Calculating Functional Beta Diversity")

    # Individual Rao's Q
    q1 = raos_quadratic_entropy(community1_traits, community1_abundances)
    q2 = raos_quadratic_entropy(community2_traits, community2_abundances)

    # Pooled community -- merge traits and abundances, halve abundances to keep
    # total relative proportions balanced between the two communities.
    total_ab1 = sum(community1_abundances)
    total_ab2 = sum(community2_abundances)

    if total_ab1 == 0 and total_ab2 == 0:
        return {"total": 0.0, "turnover": 0.0, "nestedness": 0.0}

    # Scale so both communities contribute equally to the pool
    scale1 = 1.0 / (2.0 * total_ab1) if total_ab1 > 0 else 0.0
    scale2 = 1.0 / (2.0 * total_ab2) if total_ab2 > 0 else 0.0

    pooled_traits = list(community1_traits) + list(community2_traits)
    pooled_abundances = [a * scale1 for a in community1_abundances] + [
        a * scale2 for a in community2_abundances
    ]

    # Renormalise pooled abundances so they sum to a positive total for Rao
    pool_total = sum(pooled_abundances)
    if pool_total > 0:
        pooled_abundances = [a / pool_total for a in pooled_abundances]
    # Make them look like "counts" by multiplying -- Rao already normalises
    pooled_abundances = [a * 100 for a in pooled_abundances]

    q_pool = raos_quadratic_entropy(pooled_traits, pooled_abundances)

    # Total beta = pooled Q minus mean of individual Q values
    mean_q = (q1 + q2) / 2.0
    total_beta = max(0.0, q_pool - mean_q)

    # Decompose: nestedness is driven by the difference between the two Q values
    q_min = min(q1, q2)
    q_max = max(q1, q2)
    nestedness_component = q_max - q_min
    if total_beta > 0:
        # Cap nestedness at total_beta
        nestedness_component = min(nestedness_component, total_beta)
        turnover_component = total_beta - nestedness_component
    else:
        nestedness_component = 0.0
        turnover_component = 0.0

    return {
        "total": total_beta,
        "turnover": turnover_component,
        "nestedness": nestedness_component,
    }


# ---------------------------------------------------------------------------
# Convenience / composite functions
# ---------------------------------------------------------------------------


def functional_diversity_suite(
    trait_matrix: List[List[float]], abundances: List[float]
) -> Dict[str, Any]:
    """Compute all major functional diversity indices in a single pass.

    This is a convenience wrapper that returns FRic, FEve, FDiv, FDis, Rao's Q,
    CWM, and Functional Redundancy in one dictionary.

    Args:
        trait_matrix: Species-by-traits matrix.
        abundances: Species abundance vector.

    Returns:
        Dictionary with keys: ``fric``, ``feve``, ``fdiv``, ``fdis``,
        ``raos_q``, ``cwm``, ``redundancy``.

    Example:
        >>> suite = functional_diversity_suite([[0], [5], [10]], [10, 10, 10])
        >>> sorted(suite.keys())
        ['cwm', 'fdis', 'fdiv', 'feve', 'fric', 'raos_q', 'redundancy']
    """
    validation.validate_not_empty(trait_matrix, "trait_matrix")
    validation.validate_not_empty(abundances, "abundances")
    logger.info("Computing full functional diversity suite")

    return {
        "fric": functional_richness(trait_matrix, abundances),
        "feve": functional_evenness(trait_matrix, abundances),
        "fdiv": functional_divergence(trait_matrix, abundances),
        "fdis": functional_dispersion(trait_matrix, abundances),
        "raos_q": raos_quadratic_entropy(trait_matrix, abundances),
        "cwm": community_weighted_mean(trait_matrix, abundances),
        "redundancy": functional_redundancy(trait_matrix, abundances),
    }
