"""Multi-omic pathway enrichment and active module detection.

Implements methods for combining pathway-level signals across multiple omic
layers, detecting active subnetwork modules, analysing pathway topology
effects, and assessing cross-omic concordance.

Methods:
    - Fisher's combined probability (chi-squared on -2 sum log p).
    - Stouffer's weighted Z (inverse-normal weighted combination).
    - Minimum-p with Bonferroni correction.
    - Greedy seed-and-grow active module detection.
    - Topology-weighted pathway impact analysis.
    - Cross-omic concordance scoring.
"""

from __future__ import annotations

import math
import random
from collections import defaultdict
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
    from scipy import stats as sp_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    sp_stats = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Pure-Python statistical helpers (fallbacks when scipy is absent)
# ---------------------------------------------------------------------------


def _py_chi2_sf(x: float, df: int) -> float:
    """Approximate survival function for chi-squared distribution.

    Uses the incomplete gamma function approximation for even df.
    Falls back to a rough normal approximation for large df.

    Args:
        x: Test statistic.
        df: Degrees of freedom.

    Returns:
        Approximate p-value (survival function).
    """
    if x <= 0:
        return 1.0
    if df <= 0:
        return 0.0

    # For moderate df, use a Wilson-Hilferty normal approximation
    z = ((x / df) ** (1.0 / 3) - (1.0 - 2.0 / (9 * df))) / math.sqrt(2.0 / (9 * df))
    # Standard normal survival
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def _py_norm_sf(z: float) -> float:
    """Standard normal survival function (1 - Phi(z)).

    Args:
        z: Z-score.

    Returns:
        One-sided p-value.
    """
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def _py_norm_ppf(p: float) -> float:
    """Approximate inverse normal CDF (percent point function).

    Uses Abramowitz & Stegun rational approximation (26.2.23).

    Args:
        p: Probability in (0, 1).

    Returns:
        Approximate z-score.
    """
    if p <= 0:
        return -10.0
    if p >= 1:
        return 10.0
    if p < 0.5:
        return -_py_norm_ppf(1.0 - p)

    # Rational approximation for p >= 0.5
    t = math.sqrt(-2.0 * math.log(1.0 - p))
    c0, c1, c2 = 2.515517, 0.802853, 0.010328
    d1, d2, d3 = 1.432788, 0.189269, 0.001308
    return t - (c0 + c1 * t + c2 * t**2) / (1.0 + d1 * t + d2 * t**2 + d3 * t**3)


# ---------------------------------------------------------------------------
# P-value combination methods
# ---------------------------------------------------------------------------


def _fisher_combine(p_values: list[float]) -> float:
    """Fisher's combined probability test.

    Combines p-values using the test statistic -2 * sum(log(p_i)) which
    follows a chi-squared distribution with 2k degrees of freedom.

    Args:
        p_values: List of p-values from individual tests.

    Returns:
        Combined p-value.
    """
    k = len(p_values)
    if k == 0:
        return 1.0
    # Clamp to avoid log(0)
    clamped = [max(p, 1e-300) for p in p_values]
    chi2_stat = -2.0 * sum(math.log(p) for p in clamped)
    df = 2 * k
    if HAS_SCIPY:
        return float(sp_stats.chi2.sf(chi2_stat, df))
    return _py_chi2_sf(chi2_stat, df)


def _stouffer_combine(p_values: list[float], weights: list[float] | None = None) -> float:
    """Stouffer's weighted Z method for combining p-values.

    Converts each p-value to a Z-score, takes a (weighted) sum, and
    converts back to a p-value.

    Args:
        p_values: List of p-values.
        weights: Optional weights for each p-value (default: equal weights).

    Returns:
        Combined p-value.
    """
    k = len(p_values)
    if k == 0:
        return 1.0
    if weights is None:
        weights = [1.0] * k

    z_scores: list[float] = []
    for p in p_values:
        p_clamped = max(min(p, 1.0 - 1e-15), 1e-15)
        if HAS_SCIPY:
            z = float(sp_stats.norm.ppf(1.0 - p_clamped))
        else:
            z = _py_norm_ppf(1.0 - p_clamped)
        z_scores.append(z)

    weighted_z = sum(w * z for w, z in zip(weights, z_scores))
    denom = math.sqrt(sum(w**2 for w in weights))
    if denom == 0:
        return 1.0
    combined_z = weighted_z / denom

    if HAS_SCIPY:
        return float(sp_stats.norm.sf(combined_z))
    return _py_norm_sf(combined_z)


def _min_p_combine(p_values: list[float]) -> float:
    """Minimum-p method with Bonferroni correction.

    Takes the minimum p-value and applies Bonferroni correction.

    Args:
        p_values: List of p-values.

    Returns:
        Bonferroni-corrected minimum p-value.
    """
    if not p_values:
        return 1.0
    return min(1.0, min(p_values) * len(p_values))


# ---------------------------------------------------------------------------
# Multi-omic enrichment
# ---------------------------------------------------------------------------


def multi_omic_enrichment(
    gene_sets: dict[str, list[str]],
    omic_results: dict[str, dict[str, float]],
    method: str = "fisher_combined",
) -> list[dict[str, Any]]:
    """Combine p-values from multiple omic analyses per pathway.

    For each pathway in *gene_sets*, collects the per-omic p-values from
    *omic_results* and combines them using the specified method.

    Args:
        gene_sets: Mapping from pathway_id to list of gene identifiers
            belonging to that pathway.
        omic_results: Mapping from omic name to a dictionary of
            {gene_id: p_value}.  Each gene's p-value represents its
            significance in that omic layer.
        method: P-value combination strategy.  One of ``"fisher_combined"``,
            ``"stouffer"``, or ``"min_p"``.

    Returns:
        List of dictionaries (one per pathway, sorted by combined_p), each
        containing:
            pathway_id: Pathway identifier.
            pathway_name: Same as pathway_id (caller may augment with names).
            combined_p: Combined p-value across omics.
            per_omic_p: Per-omic p-values {omic: p}.
            n_genes: Number of genes in the pathway.
            leading_edge: Genes contributing most to the signal.

    Raises:
        ValueError: If method is unsupported or inputs are empty.
    """
    valid_methods = {"fisher_combined", "stouffer", "min_p"}
    if method not in valid_methods:
        raise ValueError(f"method must be one of {valid_methods}, got '{method}'")
    if not gene_sets:
        raise ValueError("gene_sets must not be empty")
    if not omic_results:
        raise ValueError("omic_results must not be empty")

    logger.info(
        "Running multi-omic enrichment: %d pathways, %d omics, method=%s",
        len(gene_sets),
        len(omic_results),
        method,
    )

    combine_fn = {
        "fisher_combined": _fisher_combine,
        "stouffer": lambda ps: _stouffer_combine(ps),
        "min_p": _min_p_combine,
    }[method]

    results: list[dict[str, Any]] = []

    for pathway_id, genes in gene_sets.items():
        per_omic_p: dict[str, float] = {}
        leading_edge_scores: dict[str, float] = {}

        for omic_name, gene_pvals in omic_results.items():
            # Collect p-values for genes in this pathway
            pathway_pvals: list[float] = []
            for gene in genes:
                if gene in gene_pvals:
                    pathway_pvals.append(gene_pvals[gene])

            if pathway_pvals:
                # Per-omic p-value: Fisher's method within this omic for the pathway
                per_omic_p[omic_name] = _fisher_combine(pathway_pvals)
                # Track the most significant gene per omic for leading edge
                min_gene_p = min(
                    ((g, gene_pvals[g]) for g in genes if g in gene_pvals),
                    key=lambda x: x[1],
                    default=(None, 1.0),
                )
                if min_gene_p[0] is not None:
                    leading_edge_scores[min_gene_p[0]] = min(leading_edge_scores.get(min_gene_p[0], 1.0), min_gene_p[1])
            else:
                per_omic_p[omic_name] = 1.0

        # Combine across omics
        omic_pvals = list(per_omic_p.values())
        combined_p = combine_fn(omic_pvals)

        # Leading edge: genes sorted by combined significance
        leading_edge = sorted(leading_edge_scores, key=lambda g: leading_edge_scores[g])[:10]

        results.append(
            {
                "pathway_id": pathway_id,
                "pathway_name": pathway_id,
                "combined_p": combined_p,
                "per_omic_p": per_omic_p,
                "n_genes": len(genes),
                "leading_edge": leading_edge,
            }
        )

    # Sort by combined p-value
    results.sort(key=lambda r: r["combined_p"])

    logger.info(
        "Multi-omic enrichment complete: %d pathways, top p=%.2e",
        len(results),
        results[0]["combined_p"] if results else 1.0,
    )

    return results


# ---------------------------------------------------------------------------
# Active module detection
# ---------------------------------------------------------------------------


def active_module_detection(
    network: dict[str, list[str]],
    scores: dict[str, float],
    alpha: float = 0.05,
    n_permutations: int = 1000,
) -> list[dict[str, Any]]:
    """Find active subnetwork modules using node scores from multi-omic integration.

    Uses a greedy seed-and-grow approach:
        1. Rank nodes by score (significance).
        2. Seed modules from the most significant nodes.
        3. Iteratively add neighbours that improve the module score.
        4. Assess significance via permutation testing.

    The module score is the sum of transformed node scores:
        score_i = -log(p_i) - threshold, where threshold = -log(alpha).

    Args:
        network: Adjacency list {node: [neighbour1, neighbour2, ...]}.
        scores: Node scores {node: p_value}.  Lower p-values indicate
            more significant nodes.
        alpha: Significance threshold for seed selection and score transform.
        n_permutations: Number of permutations for module significance.

    Returns:
        List of module dictionaries sorted by p-value, each containing:
            module_genes: List of gene/node identifiers in the module.
            module_score: Aggregate score of the module.
            p_value: Permutation-based p-value.
            omic_contributions: Placeholder dict for per-omic contributions.

    Raises:
        ValueError: If network or scores are empty.
    """
    if not network:
        raise ValueError("network must not be empty")
    if not scores:
        raise ValueError("scores must not be empty")

    logger.info(
        "Running active module detection: %d nodes, %d scored, alpha=%.4f",
        len(network),
        len(scores),
        alpha,
    )

    threshold = -math.log(max(alpha, 1e-300))

    def node_score(node: str) -> float:
        p = scores.get(node, 1.0)
        return -math.log(max(p, 1e-300)) - threshold

    # Score all nodes
    node_scores = {node: node_score(node) for node in network}

    # Seeds: nodes with positive score (significant at alpha)
    seeds = [n for n, s in node_scores.items() if s > 0]
    seeds.sort(key=lambda n: node_scores[n], reverse=True)

    if not seeds:
        logger.warning("No seeds found at alpha=%.4f; returning empty modules", alpha)
        return []

    # Greedy seed-and-grow
    used: set[str] = set()
    modules: list[dict[str, Any]] = []

    for seed in seeds:
        if seed in used:
            continue

        module_nodes: set[str] = {seed}
        module_score_val = node_scores[seed]
        used.add(seed)

        # Expand greedily
        improved = True
        while improved:
            improved = False
            best_gain = 0.0
            best_node: str | None = None

            # Candidate neighbours
            candidates: set[str] = set()
            for n in module_nodes:
                for nb in network.get(n, []):
                    if nb not in module_nodes and nb not in used:
                        candidates.add(nb)

            for candidate in candidates:
                gain = node_scores.get(candidate, -threshold)
                if gain > best_gain:
                    best_gain = gain
                    best_node = candidate

            if best_node is not None and best_gain > 0:
                module_nodes.add(best_node)
                used.add(best_node)
                module_score_val += best_gain
                improved = True

        if len(module_nodes) < 2:
            continue

        # Permutation test for module significance
        all_nodes = list(scores.keys())
        module_size = len(module_nodes)
        perm_count = 0
        rng = random.Random(42)

        for _ in range(n_permutations):
            if len(all_nodes) >= module_size:
                perm_nodes = rng.sample(all_nodes, module_size)
            else:
                perm_nodes = all_nodes
            perm_score = sum(node_scores.get(n, -threshold) for n in perm_nodes)
            if perm_score >= module_score_val:
                perm_count += 1

        p_value = (perm_count + 1) / (n_permutations + 1)

        modules.append(
            {
                "module_genes": sorted(module_nodes),
                "module_score": module_score_val,
                "p_value": p_value,
                "omic_contributions": {},
            }
        )

    # Sort by p-value
    modules.sort(key=lambda m: m["p_value"])

    logger.info("Active module detection complete: %d modules found", len(modules))
    return modules


# ---------------------------------------------------------------------------
# Pathway topology analysis
# ---------------------------------------------------------------------------


def pathway_topology_analysis(
    pathway_graph: dict[str, list[str]],
    gene_scores: dict[str, float],
) -> dict[str, Any]:
    """Topology-based pathway analysis incorporating graph structure.

    Weights gene-level perturbation scores by the topological importance
    of each gene in the pathway graph (degree centrality and downstream
    propagation).  Computes an impact factor that accounts for both the
    significance of individual genes and their network position.

    The impact factor is:
        IF = sum_i(score_i * centrality_i) / sum_i(centrality_i)

    where centrality_i is the normalised degree centrality and score_i
    is -log10(p_value) for significant genes.

    Args:
        pathway_graph: Adjacency list {gene: [downstream_genes]}.
        gene_scores: Per-gene scores {gene: p_value}.

    Returns:
        Dictionary with keys:
            impact_factor: Topology-weighted impact score.
            p_value: Approximate p-value for the impact factor (permutation).
            perturbed_genes: List of significantly perturbed genes.
            pathway_perturbation: Per-gene perturbation details
                {gene: {score, centrality, weighted_score}}.

    Raises:
        ValueError: If inputs are empty.
    """
    if not pathway_graph:
        raise ValueError("pathway_graph must not be empty")
    if not gene_scores:
        raise ValueError("gene_scores must not be empty")

    logger.info(
        "Running pathway topology analysis: %d nodes, %d scored",
        len(pathway_graph),
        len(gene_scores),
    )

    # Compute degree centrality (in-degree + out-degree)
    degree: dict[str, int] = defaultdict(int)
    all_nodes: set[str] = set(pathway_graph.keys())
    for node, neighbours in pathway_graph.items():
        degree[node] += len(neighbours)
        for nb in neighbours:
            degree[nb] += 1
            all_nodes.add(nb)

    max_degree = max(degree.values()) if degree else 1
    centrality: dict[str, float] = {}
    for node in all_nodes:
        centrality[node] = (degree.get(node, 0) + 1) / (max_degree + 1)

    # Compute per-gene perturbation
    perturbation: dict[str, dict[str, float]] = {}
    perturbed_genes: list[str] = []
    sum_weighted = 0.0
    sum_centrality = 0.0

    for gene in all_nodes:
        p = gene_scores.get(gene, 1.0)
        score = -math.log10(max(p, 1e-300))
        c = centrality[gene]
        weighted = score * c
        perturbation[gene] = {
            "score": score,
            "centrality": c,
            "weighted_score": weighted,
        }
        sum_weighted += weighted
        sum_centrality += c

        if p < 0.05:
            perturbed_genes.append(gene)

    impact_factor = sum_weighted / sum_centrality if sum_centrality > 0 else 0.0

    # Permutation p-value
    rng = random.Random(42)
    n_perms = 1000
    perm_count = 0
    score_list = list(gene_scores.values())

    for _ in range(n_perms):
        rng.shuffle(score_list)
        perm_weighted = 0.0
        perm_cent = 0.0
        for i, gene in enumerate(all_nodes):
            p = score_list[i % len(score_list)] if score_list else 1.0
            s = -math.log10(max(p, 1e-300))
            c = centrality[gene]
            perm_weighted += s * c
            perm_cent += c
        perm_if = perm_weighted / perm_cent if perm_cent > 0 else 0.0
        if perm_if >= impact_factor:
            perm_count += 1

    p_value = (perm_count + 1) / (n_perms + 1)

    logger.info(
        "Pathway topology analysis complete: IF=%.4f, p=%.4e, %d perturbed genes",
        impact_factor,
        p_value,
        len(perturbed_genes),
    )

    return {
        "impact_factor": impact_factor,
        "p_value": p_value,
        "perturbed_genes": sorted(perturbed_genes),
        "pathway_perturbation": perturbation,
    }


# ---------------------------------------------------------------------------
# Cross-omic pathway concordance
# ---------------------------------------------------------------------------


def cross_omic_pathway_concordance(
    pathway_results: dict[str, dict[str, float]],
) -> dict[str, Any]:
    """Assess concordance of pathway signals across omic layers.

    For each pathway, compares its significance across all omic layers to
    determine whether the signal is concordant (consistently significant or
    non-significant across omics) or discordant (significant in some omics
    but not others).

    Concordance is measured by computing the coefficient of variation (CV)
    of -log10(p) across omics.  Low CV indicates concordance; high CV
    indicates discordance.

    Args:
        pathway_results: Mapping from omic_name to {pathway_id: p_value}.

    Returns:
        Dictionary with keys:
            concordant_pathways: List of pathway IDs with concordant signals
                (CV < median CV).
            discordant_pathways: List of pathway IDs with discordant signals.
            concordance_score: Global concordance score (1 - mean_CV).
            heatmap_data: Dictionary for building a heatmap:
                {pathway_id: {omic: -log10(p)}} sorted by concordance.

    Raises:
        ValueError: If pathway_results is empty or has fewer than 2 omics.
    """
    if not pathway_results:
        raise ValueError("pathway_results must not be empty")
    if len(pathway_results) < 2:
        raise ValueError("Need at least 2 omic layers for concordance analysis")

    logger.info(
        "Assessing cross-omic concordance across %d omic layers",
        len(pathway_results),
    )

    omic_names = list(pathway_results.keys())

    # Collect all pathways
    all_pathways: set[str] = set()
    for pvals in pathway_results.values():
        all_pathways.update(pvals.keys())

    if not all_pathways:
        raise ValueError("No pathways found in any omic layer")

    # Compute -log10(p) per pathway per omic and CV
    heatmap_data: dict[str, dict[str, float]] = {}
    cv_per_pathway: dict[str, float] = {}

    for pathway in all_pathways:
        neg_log_p: dict[str, float] = {}
        values: list[float] = []
        for omic in omic_names:
            p = pathway_results[omic].get(pathway, 1.0)
            val = -math.log10(max(p, 1e-300))
            neg_log_p[omic] = val
            values.append(val)

        heatmap_data[pathway] = neg_log_p

        mean_val = sum(values) / len(values) if values else 0.0
        if mean_val > 0:
            std_val = math.sqrt(sum((v - mean_val) ** 2 for v in values) / len(values))
            cv_per_pathway[pathway] = std_val / mean_val
        else:
            cv_per_pathway[pathway] = 0.0

    # Classify concordant vs discordant
    cv_values = list(cv_per_pathway.values())
    if cv_values:
        cv_values_sorted = sorted(cv_values)
        median_cv = cv_values_sorted[len(cv_values_sorted) // 2]
    else:
        median_cv = 0.0

    concordant: list[str] = []
    discordant: list[str] = []

    for pathway, cv in cv_per_pathway.items():
        if cv <= median_cv:
            concordant.append(pathway)
        else:
            discordant.append(pathway)

    concordant.sort()
    discordant.sort()

    # Global concordance score
    mean_cv = sum(cv_values) / len(cv_values) if cv_values else 0.0
    concordance_score = max(0.0, 1.0 - mean_cv)

    # Sort heatmap data by concordance (low CV first)
    sorted_pathways = sorted(all_pathways, key=lambda pw: cv_per_pathway.get(pw, 0.0))
    sorted_heatmap: dict[str, dict[str, float]] = {pw: heatmap_data[pw] for pw in sorted_pathways}

    logger.info(
        "Concordance analysis complete: %d concordant, %d discordant, score=%.4f",
        len(concordant),
        len(discordant),
        concordance_score,
    )

    return {
        "concordant_pathways": concordant,
        "discordant_pathways": discordant,
        "concordance_score": concordance_score,
        "heatmap_data": sorted_heatmap,
    }
