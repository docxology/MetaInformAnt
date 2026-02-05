"""Phylogenetic diversity metrics and community structure analysis.

Provides Faith's Phylogenetic Diversity, UniFrac distances (weighted and
unweighted), Net Relatedness Index (NRI), Nearest Taxon Index (NTI),
phylogenetic signal tests (Blomberg's K, Pagel's lambda), and simple
tree construction (UPGMA/NJ) from distance matrices.

All algorithms are pure Python implementations using nested-dict tree
representations::

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
from typing import Any

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
    total = node.get("branch_length", 0.0)
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


def compute_unifrac(
    tree: dict,
    abundances_a: dict,
    abundances_b: dict,
    weighted: bool = False,
) -> float:
    """Compute UniFrac distance between two communities.

    Unweighted UniFrac measures the fraction of branch length unique to
    either community. Weighted UniFrac accounts for abundance differences.

    Args:
        tree: Phylogenetic tree as nested dict.
        abundances_a: Dict mapping taxon name to abundance in community A.
        abundances_b: Dict mapping taxon name to abundance in community B.
        weighted: If ``True``, compute weighted UniFrac.

    Returns:
        UniFrac distance in [0, 1].
    """
    all_tips = set(_get_all_tips(tree))
    edges = _get_all_edges(tree)

    if not weighted:
        # Unweighted UniFrac
        taxa_a = set(t for t, v in abundances_a.items() if v > 0 and t in all_tips)
        taxa_b = set(t for t, v in abundances_b.items() if v > 0 and t in all_tips)

        unique_bl = 0.0
        total_bl = 0.0

        for _, child_name, bl in edges:
            # Find tips descending from this edge
            # Simple approach: check if child_name is a tip
            # For internal nodes, collect tips below
            tips_below = _tips_below_node(tree, child_name)

            has_a = bool(tips_below.intersection(taxa_a))
            has_b = bool(tips_below.intersection(taxa_b))

            total_bl += bl
            if has_a != has_b:  # unique to one community
                unique_bl += bl

        if total_bl == 0:
            return 0.0
        return unique_bl / total_bl

    else:
        # Weighted UniFrac
        total_a = sum(abundances_a.values()) or 1.0
        total_b = sum(abundances_b.values()) or 1.0

        numerator = 0.0
        denominator = 0.0

        for _, child_name, bl in edges:
            tips_below = _tips_below_node(tree, child_name)

            prop_a = sum(abundances_a.get(t, 0.0) for t in tips_below) / total_a
            prop_b = sum(abundances_b.get(t, 0.0) for t in tips_below) / total_b

            numerator += bl * abs(prop_a - prop_b)
            denominator += bl * (prop_a + prop_b)

        if denominator == 0:
            return 0.0
        return numerator / denominator


def _tips_below_node(tree: dict, node_name: str) -> set[str]:
    """Find all tips below a named node in the tree."""
    result = _find_node(tree, node_name)
    if result is None:
        return set()
    return set(_get_all_tips(result))


def _find_node(tree: dict, name: str) -> dict | None:
    """Find a node by name in the tree."""
    if tree.get("name") == name:
        return tree
    for child in tree.get("children", []):
        result = _find_node(child, name)
        if result is not None:
            return result
    return None


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
            ``"weighted_unifrac"`` (or ``"unifrac"`` for unweighted).

    Returns:
        Dictionary with keys:
            - distance_matrix: 2D list of distances between communities.
            - metric: Metric used.
            - n_communities: Number of communities.
    """
    weighted = "weighted" in metric.lower()
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

    Returns:
        Dictionary with keys:
            - nri: List of NRI values per community.
            - nti: List of NTI values per community.
            - p_values: Dict with ``"nri"`` and ``"nti"`` lists of p-values.
            - ses_values: Dict with ``"nri"`` and ``"nti"`` lists of SES values.
    """
    all_tips = _get_all_tips(tree)
    pw_dist = _pairwise_distances(tree)

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

    for comm in communities:
        valid_taxa = [t for t in comm if t in all_tips]
        if len(valid_taxa) < 2:
            nri_vals.append(0.0)
            nti_vals.append(0.0)
            p_nri.append(1.0)
            p_nti.append(1.0)
            continue

        obs_mpd = _mpd(valid_taxa)
        obs_mntd = _mntd(valid_taxa)

        null_mpds = []
        null_mntds = []
        for _ in range(n_randomizations):
            rand_taxa = random.sample(all_tips, min(len(valid_taxa), len(all_tips)))
            null_mpds.append(_mpd(rand_taxa))
            null_mntds.append(_mntd(rand_taxa))

        mean_null_mpd = sum(null_mpds) / len(null_mpds)
        sd_null_mpd = math.sqrt(sum((x - mean_null_mpd) ** 2 for x in null_mpds) / max(len(null_mpds) - 1, 1))

        mean_null_mntd = sum(null_mntds) / len(null_mntds)
        sd_null_mntd = math.sqrt(sum((x - mean_null_mntd) ** 2 for x in null_mntds) / max(len(null_mntds) - 1, 1))

        # NRI = -(obs - mean_null) / sd_null  (negative sign: clustering = positive)
        nri = -(obs_mpd - mean_null_mpd) / sd_null_mpd if sd_null_mpd > 0 else 0.0
        nti = -(obs_mntd - mean_null_mntd) / sd_null_mntd if sd_null_mntd > 0 else 0.0

        # P-values (proportion of null <= observed for clustering)
        p_nri_val = sum(1 for x in null_mpds if x <= obs_mpd) / len(null_mpds)
        p_nti_val = sum(1 for x in null_mntds if x <= obs_mntd) / len(null_mntds)

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
) -> dict:
    """Test for phylogenetic signal in a continuous trait.

    Args:
        tree: Phylogenetic tree as nested dict.
        trait_values: Dict mapping taxon name to trait value.
        method: ``"blomberg_k"`` or ``"pagel_lambda"``.

    Returns:
        Dictionary with keys:
            - statistic: Test statistic value.
            - p_value: P-value from randomization test.
            - method: Method used.
            - n_taxa: Number of taxa with trait data.
    """
    all_tips = _get_all_tips(tree)
    valid_taxa = [t for t in all_tips if t in trait_values]

    if len(valid_taxa) < 3:
        raise ValueError(f"Need at least 3 taxa with trait values, got {len(valid_taxa)}")

    pw_dist = _pairwise_distances(tree)
    n = len(valid_taxa)
    traits = [float(trait_values[t]) for t in valid_taxa]
    trait_mean = sum(traits) / n

    if method == "blomberg_k":
        return _blomberg_k(valid_taxa, traits, trait_mean, pw_dist, n)
    elif method == "pagel_lambda":
        return _pagel_lambda(valid_taxa, traits, trait_mean, pw_dist, n)
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'blomberg_k' or 'pagel_lambda'.")


def _blomberg_k(
    taxa: list[str],
    traits: list[float],
    trait_mean: float,
    pw_dist: dict[tuple[str, str], float],
    n: int,
) -> dict:
    """Compute Blomberg's K statistic."""
    # MSE0: observed mean squared error
    mse0 = sum((t - trait_mean) ** 2 for t in traits) / (n - 1)
    if mse0 == 0:
        return {"statistic": 0.0, "p_value": 1.0, "method": "blomberg_k", "n_taxa": n}

    # MSE: phylogenetic mean squared error (average of squared differences
    # weighted by inverse distance)
    total_weighted_diff = 0.0
    total_weight = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            d = pw_dist.get((taxa[i], taxa[j]), 1.0)
            if d > 0:
                w = 1.0 / d
                total_weighted_diff += w * (traits[i] - traits[j]) ** 2
                total_weight += w

    mse_phylo = total_weighted_diff / total_weight if total_weight > 0 else mse0

    k_obs = (mse0 / mse_phylo) if mse_phylo > 0 else 0.0

    # Randomization test
    n_rand = 999
    count_ge = 0
    for _ in range(n_rand):
        perm_traits = list(traits)
        random.shuffle(perm_traits)
        perm_mse0 = sum((t - trait_mean) ** 2 for t in perm_traits) / (n - 1)

        tw_diff = 0.0
        tw = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                d = pw_dist.get((taxa[i], taxa[j]), 1.0)
                if d > 0:
                    w = 1.0 / d
                    tw_diff += w * (perm_traits[i] - perm_traits[j]) ** 2
                    tw += w

        perm_mse_phylo = tw_diff / tw if tw > 0 else perm_mse0
        k_perm = (perm_mse0 / perm_mse_phylo) if perm_mse_phylo > 0 else 0.0

        if k_perm >= k_obs:
            count_ge += 1

    p_value = (count_ge + 1) / (n_rand + 1)

    logger.info("Blomberg's K = %.4f, p = %.4f (n=%d)", k_obs, p_value, n)

    return {
        "statistic": k_obs,
        "p_value": p_value,
        "method": "blomberg_k",
        "n_taxa": n,
    }


def _pagel_lambda(
    taxa: list[str],
    traits: list[float],
    trait_mean: float,
    pw_dist: dict[tuple[str, str], float],
    n: int,
) -> dict:
    """Approximate Pagel's lambda using correlation of trait differences with distances."""
    # Compute correlation between |trait_i - trait_j| and phylo distance
    diffs = []
    dists = []
    for i in range(n):
        for j in range(i + 1, n):
            diffs.append(abs(traits[i] - traits[j]))
            dists.append(pw_dist.get((taxa[i], taxa[j]), 0.0))

    if len(diffs) < 3:
        return {"statistic": 0.0, "p_value": 1.0, "method": "pagel_lambda", "n_taxa": n}

    # Pearson correlation
    mean_diff = sum(diffs) / len(diffs)
    mean_dist = sum(dists) / len(dists)
    sd_diff = math.sqrt(sum((d - mean_diff) ** 2 for d in diffs) / (len(diffs) - 1))
    sd_dist = math.sqrt(sum((d - mean_dist) ** 2 for d in dists) / (len(dists) - 1))

    if sd_diff == 0 or sd_dist == 0:
        lam = 0.0
    else:
        r = sum((di - mean_diff) * (dj - mean_dist) for di, dj in zip(diffs, dists)) / (
            (len(diffs) - 1) * sd_diff * sd_dist
        )
        # Lambda approximation: higher correlation -> higher signal -> lambda near 1
        lam = max(0.0, min(1.0, r))

    # Randomization test
    n_rand = 999
    count_ge = 0
    for _ in range(n_rand):
        perm = list(traits)
        random.shuffle(perm)
        perm_diffs = []
        for i in range(n):
            for j in range(i + 1, n):
                perm_diffs.append(abs(perm[i] - perm[j]))

        pm = sum(perm_diffs) / len(perm_diffs)
        ps = math.sqrt(sum((d - pm) ** 2 for d in perm_diffs) / (len(perm_diffs) - 1))
        if ps > 0 and sd_dist > 0:
            r_perm = sum((di - pm) * (dj - mean_dist) for di, dj in zip(perm_diffs, dists)) / (
                (len(perm_diffs) - 1) * ps * sd_dist
            )
            lam_perm = max(0.0, min(1.0, r_perm))
        else:
            lam_perm = 0.0

        if lam_perm >= lam:
            count_ge += 1

    p_value = (count_ge + 1) / (n_rand + 1)

    logger.info("Pagel's lambda = %.4f, p = %.4f (n=%d)", lam, p_value, n)

    return {
        "statistic": lam,
        "p_value": p_value,
        "method": "pagel_lambda",
        "n_taxa": n,
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

        # Branch lengths
        denom = 2.0 * (m - 2) if m > 2 else 2.0
        bl_i = d[mi][mj] / 2.0 + (r[mi] - r[mj]) / denom
        bl_j = d[mi][mj] - bl_i
        bl_i = max(0.0, bl_i)
        bl_j = max(0.0, bl_j)

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
