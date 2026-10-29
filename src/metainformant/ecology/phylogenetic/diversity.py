"""Phylogenetic diversity metrics and community structure analysis.

Provides Faith's Phylogenetic Diversity, UniFrac distances (weighted and
unweighted), Net Relatedness Index (NRI), Nearest Taxon Index (NTI),
phylogenetic signal tests (Blomberg's K, Pagel's lambda), and simple
tree construction (UPGMA/NJ) from distance matrices.

The diversity metrics are pure Python implementations using nested-dict tree
representations; the phylogenetic-signal tests also use NumPy and SciPy for the
linear algebra of the tree's covariance matrix::

    tree = {
        "name": "root",
        "branch_length": 0.0,
        "children": [
            {"name": "A", "branch_length": 0.5, "children": []},
            {"name": "internal", "branch_length": 0.3, "children": [
                {"name": "B", "branch_length": 0.2, "children": []},
                {"name": "C", "branch_length": 0.4, "children": []},
            ]},
        ],
    }
"""

from __future__ import annotations

import math
import random

import numpy as np
from scipy import linalg, optimize, stats

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Tree helpers
# ---------------------------------------------------------------------------


def _get_all_tips(node: dict) -> list[str]:
    """Recursively collect all tip (leaf) names from a tree node."""
    if not node.get("children"):
        return [node["name"]]
    tips: list[str] = []
    for child in node["children"]:
        tips.extend(_get_all_tips(child))
    return tips


def _get_all_edges(node: dict, parent_name: str | None = None) -> list[tuple[str, str, float]]:
    """Collect all edges as (parent, child, branch_length) tuples."""
    edges: list[tuple[str, str, float]] = []
    node_id = node.get("name", "unnamed")
    bl = node.get("branch_length", 0.0)
    if parent_name is not None:
        edges.append((parent_name, node_id, bl))
    for child in node.get("children", []):
        edges.extend(_get_all_edges(child, node_id))
    return edges


def _pairwise_distances(tree: dict) -> dict[tuple[str, str], float]:
    """Compute pairwise tip distances through the tree.

    Uses a simple DFS approach: for each tip pair, finds the path
    through the tree and sums branch lengths.
    """
    tips = _get_all_tips(tree)
    # Build adjacency list from edges
    edges = _get_all_edges(tree)
    adj: dict[str, list[tuple[str, float]]] = {}
    for parent, child, bl in edges:
        adj.setdefault(parent, []).append((child, bl))
        adj.setdefault(child, []).append((parent, bl))

    distances: dict[tuple[str, str], float] = {}

    for i, tip_a in enumerate(tips):
        # BFS from tip_a
        dist_from_a: dict[str, float] = {tip_a: 0.0}
        queue = [tip_a]
        visited = {tip_a}
        while queue:
            current = queue.pop(0)
            for neighbour, bl in adj.get(current, []):
                if neighbour not in visited:
                    visited.add(neighbour)
                    dist_from_a[neighbour] = dist_from_a[current] + bl
                    queue.append(neighbour)

        for j in range(i + 1, len(tips)):
            tip_b = tips[j]
            d = dist_from_a.get(tip_b, 0.0)
            distances[(tip_a, tip_b)] = d
            distances[(tip_b, tip_a)] = d

    return distances


def _total_branch_length(node: dict) -> float:
    """Compute total branch length of the tree."""
    total: float = node.get("branch_length", 0.0)
    for child in node.get("children", []):
        total += _total_branch_length(child)
    return total


def _branch_lengths_connecting(node: dict, taxa_set: set[str]) -> float:
    """Compute sum of branch lengths connecting taxa in taxa_set to the root.

    A branch is included if at least one taxon in taxa_set descends from it.
    """
    tips_below = set(_get_all_tips(node))
    # If no intersection, this branch doesn't contribute
    if not tips_below.intersection(taxa_set):
        return 0.0

    total = 0.0
    for child in node.get("children", []):
        child_tips = set(_get_all_tips(child))
        if child_tips.intersection(taxa_set):
            total += child.get("branch_length", 0.0)
            total += _branch_lengths_connecting(child, taxa_set)

    return total


# ---------------------------------------------------------------------------
# Faith's Phylogenetic Diversity
# ---------------------------------------------------------------------------


def faiths_pd(tree: dict, taxa_present: list[str]) -> float:
    """Compute Faith's Phylogenetic Diversity.

    Faith's PD is the sum of all branch lengths on the minimum spanning
    path connecting the taxa present in the community to the root of the
    phylogeny.

    Args:
        tree: Phylogenetic tree as nested dict with keys ``name``,
            ``branch_length``, ``children``.
        taxa_present: List of tip names present in the community.

    Returns:
        Faith's PD value (sum of branch lengths).

    Raises:
        ValueError: If no taxa from the list are found in the tree.
    """
    all_tips = set(_get_all_tips(tree))
    present = set(taxa_present)
    valid = present.intersection(all_tips)

    if not valid:
        raise ValueError(f"None of the taxa {taxa_present[:5]}... found in tree tips: " f"{sorted(all_tips)[:10]}...")

    pd_value = _branch_lengths_connecting(tree, valid)
    logger.info(
        "Faith's PD for %d/%d taxa: %.4f",
        len(valid),
        len(all_tips),
        pd_value,
    )
    return pd_value


# ---------------------------------------------------------------------------
# UniFrac
# ---------------------------------------------------------------------------


def _edges_with_tips(node: dict, is_root: bool = True) -> list[tuple[float, frozenset[str]]]:
    """Every branch below the root as ``(branch_length, tips descending from it)``.

    Walking the nested dicts directly, instead of looking each branch up by name, keeps trees whose internal
    nodes share a name (or have none) correct.
    """
    branches: list[tuple[float, frozenset[str]]] = []
    for child in node.get("children", []):
        branches.extend(_edges_with_tips(child, is_root=False))
    if not is_root:
        branches.append((node.get("branch_length", 0.0), frozenset(_get_all_tips(node))))
    return branches


def compute_unifrac(
    tree: dict,
    abundances_a: dict,
    abundances_b: dict,
    weighted: bool = False,
) -> float:
    """Compute UniFrac distance between two communities.

    Unweighted UniFrac is the fraction of the branch length spanned by either community that leads to
    descendants of only one of them: ``sum(l_i * |A_i - B_i|) / sum(l_i * max(A_i, B_i))`` for branches ``i`` with
    presence indicators ``A_i`` and ``B_i`` (Lozupone and Knight 2005). Branches with no descendant in either
    community are outside the denominator; counting them (the whole tree's length) shrinks the distance towards 0
    for every taxon the tree holds beyond the two communities. Weighted UniFrac accounts for abundance
    differences (normalised weighted UniFrac, Lozupone et al. 2007).

    Args:
        tree: Phylogenetic tree as nested dict.
        abundances_a: Dict mapping taxon name to abundance in community A.
        abundances_b: Dict mapping taxon name to abundance in community B.
        weighted: If ``True``, compute weighted UniFrac.

    Returns:
        UniFrac distance in [0, 1]; 0.0 when neither community has a taxon on a branch with positive length.
    """
    all_tips = set(_get_all_tips(tree))
    branches = _edges_with_tips(tree)

    if not weighted:
        # Unweighted UniFrac
        taxa_a = set(t for t, v in abundances_a.items() if v > 0 and t in all_tips)
        taxa_b = set(t for t, v in abundances_b.items() if v > 0 and t in all_tips)

        unique_bl = 0.0
        observed_bl = 0.0

        for bl, tips_below in branches:
            has_a = not tips_below.isdisjoint(taxa_a)
            has_b = not tips_below.isdisjoint(taxa_b)
            if not (has_a or has_b):
                continue  # a branch neither community reaches
            observed_bl += bl
            if has_a != has_b:  # unique to one community
                unique_bl += bl

        if observed_bl == 0:
            return 0.0
        return unique_bl / observed_bl

    else:
        # Weighted UniFrac
        total_a = sum(abundances_a.values()) or 1.0
        total_b = sum(abundances_b.values()) or 1.0

        numerator = 0.0
        denominator = 0.0

        for bl, tips_below in branches:
            prop_a = sum(abundances_a.get(t, 0.0) for t in tips_below) / total_a
            prop_b = sum(abundances_b.get(t, 0.0) for t in tips_below) / total_b

            numerator += bl * abs(prop_a - prop_b)
            denominator += bl * (prop_a + prop_b)

        if denominator == 0:
            return 0.0
        return numerator / denominator


# ---------------------------------------------------------------------------
# Phylogenetic beta diversity
# ---------------------------------------------------------------------------


def phylogenetic_beta_diversity(
    tree: dict,
    communities: list[list[str]],
    metric: str = "unifrac",
) -> dict:
    """Compute phylogenetic beta diversity between multiple communities.

    Args:
        tree: Phylogenetic tree as nested dict.
        communities: List of communities, each a list of taxon names.
        metric: Distance metric: ``"unweighted_unifrac"`` or
            ``"weighted_unifrac"`` (or ``"unifrac"`` for unweighted). Communities are presence/absence lists, so
            weighted UniFrac here weights every listed taxon equally.

    Returns:
        Dictionary with keys:
            - distance_matrix: 2D list of distances between communities.
            - metric: Metric used.
            - n_communities: Number of communities.

    Raises:
        ValueError: If ``metric`` is none of these.
    """
    # "unweighted_unifrac" contains the substring "weighted", so the metric name has to be matched exactly:
    # a substring test computed weighted UniFrac for the metric that asks for unweighted.
    metric_key = metric.lower()
    if metric_key in ("unifrac", "unweighted_unifrac"):
        weighted = False
    elif metric_key == "weighted_unifrac":
        weighted = True
    else:
        raise ValueError(f"Unsupported metric {metric!r}. Use 'unifrac', 'unweighted_unifrac' or 'weighted_unifrac'")
    n = len(communities)
    dist_matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        ab_a = {t: 1.0 for t in communities[i]}
        for j in range(i + 1, n):
            ab_b = {t: 1.0 for t in communities[j]}
            d = compute_unifrac(tree, ab_a, ab_b, weighted=weighted)
            dist_matrix[i][j] = d
            dist_matrix[j][i] = d

    logger.info("Computed phylogenetic beta diversity (%s) for %d communities", metric, n)

    return {
        "distance_matrix": dist_matrix,
        "metric": metric,
        "n_communities": n,
    }


# ---------------------------------------------------------------------------
# NRI / NTI
# ---------------------------------------------------------------------------


def nri_nti(
    tree: dict,
    communities: list[list[str]],
    n_randomizations: int = 999,
    seed: int | None = None,
) -> dict:
    """Compute Net Relatedness Index (NRI) and Nearest Taxon Index (NTI).

    NRI standardises the mean pairwise distance (MPD) relative to a null
    model. NTI standardises the mean nearest taxon distance (MNTD).

    Positive values indicate phylogenetic clustering; negative values
    indicate phylogenetic overdispersion.

    Args:
        tree: Phylogenetic tree as nested dict.
        communities: List of communities (each a list of taxon names).
        n_randomizations: Number of randomizations for null distribution.
        seed: Seed of the random number generator used for the null draws. With ``None`` the draws are
            different on every call; either way they come from a private generator, so the interpreter's
            global ``random`` state is neither read nor changed.

    Returns:
        Dictionary with keys:
            - nri: List of NRI values per community.
            - nti: List of NTI values per community.
            - p_values: Dict with ``"nri"`` and ``"nti"`` lists of one-sided (clustering) p-values,
              ``(1 + number of null values <= observed) / (1 + n_randomizations)``, so none is 0.
            - ses_values: Dict with ``"nri"`` and ``"nti"`` lists of SES values.

        A community with fewer than two taxa on the tree has no pairwise distance to standardise: its NRI, NTI
        and p-values are ``nan``, not 0 and 1. A community whose null distribution has no spread (for instance
        every pairwise distance on the tree is equal) has a zero standard deviation to divide by: its NRI or
        NTI is ``nan`` (the p-value, which does not divide by it, is still reported).
    """
    all_tips = _get_all_tips(tree)
    pw_dist = _pairwise_distances(tree)
    rng = random.Random(seed)

    def _mpd(taxa: list[str]) -> float:
        if len(taxa) < 2:
            return 0.0
        total = 0.0
        count = 0
        for i in range(len(taxa)):
            for j in range(i + 1, len(taxa)):
                d = pw_dist.get((taxa[i], taxa[j]), 0.0)
                total += d
                count += 1
        return total / count if count > 0 else 0.0

    def _mntd(taxa: list[str]) -> float:
        if len(taxa) < 2:
            return 0.0
        min_dists = []
        for i, t1 in enumerate(taxa):
            nearest = float("inf")
            for j, t2 in enumerate(taxa):
                if i != j:
                    d = pw_dist.get((t1, t2), 0.0)
                    nearest = min(nearest, d)
            if nearest < float("inf"):
                min_dists.append(nearest)
        return sum(min_dists) / len(min_dists) if min_dists else 0.0

    nri_vals: list[float] = []
    nti_vals: list[float] = []
    p_nri: list[float] = []
    p_nti: list[float] = []
    undefined = float("nan")

    for comm in communities:
        valid_taxa = [t for t in comm if t in all_tips]
        if len(valid_taxa) < 2:
            nri_vals.append(undefined)
            nti_vals.append(undefined)
            p_nri.append(undefined)
            p_nti.append(undefined)
            continue

        obs_mpd = _mpd(valid_taxa)
        obs_mntd = _mntd(valid_taxa)

        null_mpds = []
        null_mntds = []
        for _ in range(n_randomizations):
            rand_taxa = rng.sample(all_tips, min(len(valid_taxa), len(all_tips)))
            null_mpds.append(_mpd(rand_taxa))
            null_mntds.append(_mntd(rand_taxa))

        mean_null_mpd = sum(null_mpds) / len(null_mpds)
        sd_null_mpd = math.sqrt(sum((x - mean_null_mpd) ** 2 for x in null_mpds) / max(len(null_mpds) - 1, 1))

        mean_null_mntd = sum(null_mntds) / len(null_mntds)
        sd_null_mntd = math.sqrt(sum((x - mean_null_mntd) ** 2 for x in null_mntds) / max(len(null_mntds) - 1, 1))

        # NRI = -(obs - mean_null) / sd_null  (negative sign: clustering = positive); undefined without spread
        nri = -(obs_mpd - mean_null_mpd) / sd_null_mpd if sd_null_mpd > 0 else undefined
        nti = -(obs_mntd - mean_null_mntd) / sd_null_mntd if sd_null_mntd > 0 else undefined

        # One-sided (clustering) p-values: the observed value counts as one more draw from the null, so a
        # finite randomization never yields p = 0
        p_nri_val = (1 + sum(1 for x in null_mpds if x <= obs_mpd)) / (1 + len(null_mpds))
        p_nti_val = (1 + sum(1 for x in null_mntds if x <= obs_mntd)) / (1 + len(null_mntds))

        nri_vals.append(nri)
        nti_vals.append(nti)
        p_nri.append(p_nri_val)
        p_nti.append(p_nti_val)

    logger.info("NRI/NTI computed for %d communities with %d randomizations", len(communities), n_randomizations)

    return {
        "nri": nri_vals,
        "nti": nti_vals,
        "p_values": {"nri": p_nri, "nti": p_nti},
        "ses_values": {"nri": nri_vals, "nti": nti_vals},
    }


# ---------------------------------------------------------------------------
# Phylogenetic signal
# ---------------------------------------------------------------------------


def phylogenetic_signal(
    tree: dict,
    trait_values: dict,
    method: str = "blomberg_k",
    seed: int | None = None,
) -> dict:
    """Test for phylogenetic signal in a continuous trait.

    Both statistics assume the trait evolves by Brownian motion along the tree, whose covariance between two
    tips is the length of the root path they share (the root's own branch is ignored).

    * ``"blomberg_k"``: Blomberg et al.'s (2003) K, the observed ratio of the mean squared error of the trait
      about its arithmetic mean to that about its phylogenetic (generalised least squares) mean, divided by the
      ratio expected under Brownian motion. K = 1 for a Brownian trait, K < 1 means relatives resemble each
      other less than that, K > 1 more. The p-value is a one-sided randomization test (999 shuffles of the
      trait across tips): the share of shuffles whose K is at least the observed one, counting the observed
      value as one more draw.
    * ``"pagel_lambda"``: Pagel's (1999) lambda, the maximum-likelihood value in [0, 1] by which the tree's
      internal (off-diagonal) covariances are multiplied, 0 meaning no signal and 1 Brownian motion on the
      tree as given. The p-value is the likelihood-ratio test against lambda = 0 (chi-square, 1 degree of
      freedom; conservative because 0 is on the boundary of the parameter range).

    Args:
        tree: Phylogenetic tree as nested dict.
        trait_values: Dict mapping taxon name to trait value.
        method: ``"blomberg_k"`` or ``"pagel_lambda"``.
        seed: Seed for the Blomberg's K randomization (a private generator is used; the global ``random``
            state is untouched). Pagel's lambda is deterministic.

    Returns:
        Dictionary with keys:
            - statistic: Blomberg's K, or the maximum-likelihood lambda.
            - p_value: P-value as described above.
            - method: Method used.
            - n_taxa: Number of taxa with trait data.
            - log_likelihood, log_likelihood_lambda0: (``"pagel_lambda"`` only) maximised log-likelihood and its
              value at lambda = 0.

    Raises:
        ValueError: If ``method`` is unknown, fewer than 3 taxa have trait values, the trait does not vary, or
            the tree's covariance matrix is singular (or, for lambda, has no covariance between tips, so lambda is
            not identifiable).
    """
    if method not in ("blomberg_k", "pagel_lambda"):
        raise ValueError(f"Unknown method '{method}'. Use 'blomberg_k' or 'pagel_lambda'.")

    all_tips = _get_all_tips(tree)
    valid_taxa = [t for t in all_tips if t in trait_values]

    if len(valid_taxa) < 3:
        raise ValueError(f"Need at least 3 taxa with trait values, got {len(valid_taxa)}")

    traits = np.asarray([float(trait_values[t]) for t in valid_taxa], dtype=float)
    if not np.all(np.isfinite(traits)):
        raise ValueError("Trait values must be finite")
    if float(np.ptp(traits)) == 0.0:
        raise ValueError("The trait has no variation across taxa, so there is no signal to measure")

    vcv = _phylogenetic_vcv(tree, valid_taxa)
    if method == "blomberg_k":
        return _blomberg_k(traits, vcv, random.Random(seed))
    return _pagel_lambda(traits, vcv)


def _phylogenetic_vcv(tree: dict, taxa: list[str]) -> np.ndarray:
    """Covariance matrix of a Brownian-motion trait at the tips *taxa*.

    ``C[i, j]`` is the length of the root path that tips ``i`` and ``j`` share (their root-to-ancestor distance),
    and ``C[i, i]`` is the root-to-tip distance. Nodes are told apart by identity, not name.
    """
    paths: dict[str, list[tuple[int, float]]] = {}
    counter = [0]

    def walk(node: dict, depth: float, trail: list[tuple[int, float]]) -> None:
        children = node.get("children") or []
        if not children:
            paths[node["name"]] = trail
            return
        for child in children:
            counter[0] += 1
            child_depth = depth + float(child.get("branch_length", 0.0))
            walk(child, child_depth, trail + [(counter[0], child_depth)])

    walk(tree, 0.0, [])

    n = len(taxa)
    vcv = np.zeros((n, n))
    for i in range(n):
        path_i = paths[taxa[i]]
        vcv[i, i] = path_i[-1][1] if path_i else 0.0
        for j in range(i + 1, n):
            shared = 0.0
            for (node_i, depth_i), (node_j, _) in zip(path_i, paths[taxa[j]]):
                if node_i != node_j:
                    break
                shared = depth_i
            vcv[i, j] = vcv[j, i] = shared
    return vcv


def _blomberg_k(traits: np.ndarray, vcv: np.ndarray, rng: random.Random) -> dict:
    """Blomberg's K and its randomization p-value (see :func:`phylogenetic_signal`)."""
    n = len(traits)
    try:
        vcv_inv = np.linalg.inv(vcv)
    except np.linalg.LinAlgError as exc:
        raise ValueError("The tree's covariance matrix is singular (do some tips sit at the same position?)") from exc
    ones = np.ones(n)
    information = float(ones @ vcv_inv @ ones)
    # Expected MSE0 / MSE under Brownian motion: (tr(C) - n / (1' C^-1 1)) / (n - 1)
    expected_ratio = (float(np.trace(vcv)) - n / information) / (n - 1)
    if not expected_ratio > 0:
        raise ValueError("The tree has no branch length separating the tips, so K is not defined")

    def k_statistic(values: np.ndarray) -> float:
        phylogenetic_mean = float(ones @ vcv_inv @ values) / information
        centred = values - phylogenetic_mean
        mse = float(centred @ vcv_inv @ centred) / (n - 1)
        arithmetic = values - values.mean()
        mse0 = float(arithmetic @ arithmetic) / (n - 1)
        return (mse0 / mse) / expected_ratio

    k_obs = k_statistic(traits)

    n_rand = 999
    count_ge = 0
    shuffled = list(traits)
    for _ in range(n_rand):
        rng.shuffle(shuffled)
        if k_statistic(np.asarray(shuffled)) >= k_obs:
            count_ge += 1

    p_value = (count_ge + 1) / (n_rand + 1)

    logger.info("Blomberg's K = %.4f, p = %.4f (n=%d)", k_obs, p_value, n)

    return {
        "statistic": k_obs,
        "p_value": p_value,
        "method": "blomberg_k",
        "n_taxa": n,
    }


def _pagel_lambda(traits: np.ndarray, vcv: np.ndarray) -> dict:
    """Maximum-likelihood Pagel's lambda and its likelihood-ratio p-value (see :func:`phylogenetic_signal`)."""
    n = len(traits)
    diagonal = np.diag(np.diag(vcv))
    off_diagonal = vcv - diagonal
    if not np.any(off_diagonal):
        raise ValueError("Pagel's lambda is not identifiable: no two tips share any ancestry on this tree")
    ones = np.ones(n)

    def log_likelihood(lam: float) -> float:
        # diag(C) + lam * offdiag(C) is positive definite for lam in [0, 1] whenever C is, so a Cholesky factor
        # gives both the log-determinant and the solves (and fails cleanly when C is singular)
        scaled = diagonal + lam * off_diagonal
        try:
            factor = (np.linalg.cholesky(scaled), True)
        except np.linalg.LinAlgError:
            return -math.inf
        log_det = 2.0 * float(np.sum(np.log(np.diag(factor[0]))))
        solved_ones = linalg.cho_solve(factor, ones)
        generalised_mean = float(solved_ones @ traits) / float(solved_ones @ ones)
        residual = traits - generalised_mean
        sigma2 = float(residual @ linalg.cho_solve(factor, residual)) / n
        if not sigma2 > 0:
            return -math.inf
        return -0.5 * (n * math.log(2.0 * math.pi * sigma2) + log_det + n)

    grid = np.linspace(0.0, 1.0, 101)
    grid_values = [log_likelihood(float(lam)) for lam in grid]
    best = int(np.argmax(grid_values))
    if not math.isfinite(grid_values[best]):
        raise ValueError("The likelihood cannot be evaluated: the tree's covariance matrix is singular")
    refined = optimize.minimize_scalar(
        lambda lam: -log_likelihood(float(lam)),
        bounds=(float(grid[max(best - 1, 0)]), float(grid[min(best + 1, len(grid) - 1)])),
        method="bounded",
        options={"xatol": 1e-10},
    )
    if -float(refined.fun) >= grid_values[best]:
        lam_hat, ll_hat = float(refined.x), -float(refined.fun)
    else:
        lam_hat, ll_hat = float(grid[best]), grid_values[best]
    ll_null = grid_values[0]  # lambda = 0

    statistic = max(0.0, 2.0 * (ll_hat - ll_null))
    p_value = float(stats.chi2.sf(statistic, df=1))

    logger.info("Pagel's lambda = %.4f, p = %.4f (n=%d)", lam_hat, p_value, n)

    return {
        "statistic": lam_hat,
        "p_value": p_value,
        "method": "pagel_lambda",
        "n_taxa": n,
        "log_likelihood": ll_hat,
        "log_likelihood_lambda0": ll_null,
    }


# ---------------------------------------------------------------------------
# Tree construction
# ---------------------------------------------------------------------------


def build_simple_tree(
    distance_matrix: list[list[float]],
    taxa_names: list[str],
    method: str = "upgma",
) -> dict:
    """Build a phylogenetic tree from a distance matrix.

    Args:
        distance_matrix: Symmetric pairwise distance matrix (2D list).
        taxa_names: List of taxon names corresponding to matrix rows/columns.
        method: Tree-building method: ``"upgma"`` (default) or ``"nj"``
            (neighbor-joining).

    Returns:
        Tree as nested dict with keys ``name``, ``branch_length``, ``children``.

    Raises:
        ValueError: If matrix dimensions do not match taxa_names.
    """
    n = len(taxa_names)
    if len(distance_matrix) != n:
        raise ValueError(f"Distance matrix has {len(distance_matrix)} rows but {n} taxa")

    if method == "upgma":
        return _upgma(distance_matrix, taxa_names)
    elif method == "nj":
        return _neighbor_joining(distance_matrix, taxa_names)
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'upgma' or 'nj'.")


def _upgma(dist: list[list[float]], names: list[str]) -> dict:
    """UPGMA hierarchical clustering for tree construction."""
    n = len(names)
    # Create mutable distance matrix copy
    d = [row[:] for row in dist]
    # Track clusters: each is a dict representing a subtree
    clusters: list[dict] = [{"name": name, "branch_length": 0.0, "children": []} for name in names]
    cluster_sizes = [1] * n
    heights = [0.0] * n
    active = list(range(n))

    internal_count = 0

    while len(active) > 1:
        # Find minimum distance pair
        min_d = float("inf")
        mi, mj = 0, 1
        for ii in range(len(active)):
            for jj in range(ii + 1, len(active)):
                i_idx = active[ii]
                j_idx = active[jj]
                if d[i_idx][j_idx] < min_d:
                    min_d = d[i_idx][j_idx]
                    mi, mj = ii, jj

        i_idx = active[mi]
        j_idx = active[mj]

        new_height = min_d / 2.0
        bl_i = new_height - heights[i_idx]
        bl_j = new_height - heights[j_idx]

        # Set branch lengths on children
        child_i = clusters[i_idx].copy()
        child_i["branch_length"] = max(0.0, bl_i)
        child_j = clusters[j_idx].copy()
        child_j["branch_length"] = max(0.0, bl_j)

        internal_count += 1
        new_cluster = {
            "name": f"internal_{internal_count}",
            "branch_length": 0.0,
            "children": [child_i, child_j],
        }

        # Update distances using weighted average (UPGMA)
        new_size = cluster_sizes[i_idx] + cluster_sizes[j_idx]
        new_row = [0.0] * (n + internal_count)

        for k_idx in active:
            if k_idx == i_idx or k_idx == j_idx:
                continue
            new_dist = (d[i_idx][k_idx] * cluster_sizes[i_idx] + d[j_idx][k_idx] * cluster_sizes[j_idx]) / new_size
            new_row[k_idx] = new_dist

        # Add new cluster
        new_idx = n + internal_count - 1
        # Extend d matrix
        for row in d:
            row.extend([0.0] * (new_idx + 1 - len(row)))
        while len(d) <= new_idx:
            d.append([0.0] * (new_idx + 1))

        for k in range(len(d)):
            if k < len(new_row):
                d[new_idx][k] = new_row[k]
                d[k][new_idx] = new_row[k]

        while len(clusters) <= new_idx:
            clusters.append({})
        clusters[new_idx] = new_cluster

        while len(cluster_sizes) <= new_idx:
            cluster_sizes.append(0)
        cluster_sizes[new_idx] = new_size

        while len(heights) <= new_idx:
            heights.append(0.0)
        heights[new_idx] = new_height

        # Update active list
        active.remove(i_idx)
        active.remove(j_idx)
        active.append(new_idx)

    root = clusters[active[0]]
    root["name"] = "root"
    logger.info("Built UPGMA tree with %d tips", n)
    return root


def _neighbor_joining(dist: list[list[float]], names: list[str]) -> dict:
    """Neighbor-Joining tree construction."""
    n = len(names)
    d = [row[:] for row in dist]
    nodes: list[dict] = [{"name": name, "branch_length": 0.0, "children": []} for name in names]
    active = list(range(n))
    internal_count = 0

    # Extend matrix as needed
    max_size = 2 * n
    while len(d) < max_size:
        d.append([0.0] * max_size)
    for row in d:
        while len(row) < max_size:
            row.append(0.0)
    while len(nodes) < max_size:
        nodes.append({})

    while len(active) > 2:
        m = len(active)
        # Compute r_i = sum of distances
        r = {}
        for i in active:
            r[i] = sum(d[i][j] for j in active if j != i)

        # Find pair minimising Q
        min_q = float("inf")
        mi, mj = active[0], active[1]
        for ii in range(len(active)):
            for jj in range(ii + 1, len(active)):
                i_idx = active[ii]
                j_idx = active[jj]
                q = (m - 2) * d[i_idx][j_idx] - r[i_idx] - r[j_idx]
                if q < min_q:
                    min_q = q
                    mi, mj = i_idx, j_idx

        # Branch lengths: clamp bl_i first so the sibling branch keeps the
        # original distance sum (bl_i + bl_j == d[mi][mj])
        denom = 2.0 * (m - 2) if m > 2 else 2.0
        bl_i = max(0.0, d[mi][mj] / 2.0 + (r[mi] - r[mj]) / denom)
        bl_i = min(bl_i, d[mi][mj])
        bl_j = max(0.0, d[mi][mj] - bl_i)

        child_i = nodes[mi].copy()
        child_i["branch_length"] = bl_i
        child_j = nodes[mj].copy()
        child_j["branch_length"] = bl_j

        internal_count += 1
        new_idx = n + internal_count - 1
        new_node = {
            "name": f"internal_{internal_count}",
            "branch_length": 0.0,
            "children": [child_i, child_j],
        }
        nodes[new_idx] = new_node

        # Update distances
        for k in active:
            if k == mi or k == mj:
                continue
            new_dist = (d[mi][k] + d[mj][k] - d[mi][mj]) / 2.0
            d[new_idx][k] = max(0.0, new_dist)
            d[k][new_idx] = max(0.0, new_dist)

        active.remove(mi)
        active.remove(mj)
        active.append(new_idx)

    # Join final two
    if len(active) == 2:
        i_idx, j_idx = active
        bl = d[i_idx][j_idx] / 2.0
        child_i = nodes[i_idx].copy()
        child_i["branch_length"] = max(0.0, bl)
        child_j = nodes[j_idx].copy()
        child_j["branch_length"] = max(0.0, d[i_idx][j_idx] - bl)

        root = {
            "name": "root",
            "branch_length": 0.0,
            "children": [child_i, child_j],
        }
    else:
        root = nodes[active[0]]
        root["name"] = "root"

    logger.info("Built NJ tree with %d tips", n)
    return root
