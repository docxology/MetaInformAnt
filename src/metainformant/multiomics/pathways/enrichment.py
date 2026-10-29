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
import operator
import random
import statistics
from collections import defaultdict
from typing import Any

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

sp_stats: Any  # None when scipy is unavailable (optional dependency)
try:
    from scipy import stats as sp_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    sp_stats = None


# ---------------------------------------------------------------------------
# Pure-Python statistical helpers (fallbacks when scipy is absent)
# ---------------------------------------------------------------------------


def _py_chi2_sf(x: float, df: int) -> float:
    """Survival function of the chi-squared distribution without scipy.

    Exact for even ``df`` (the only kind Fisher's method uses, ``df = 2k``): the closed form
    ``exp(-x/2) * sum_{j<k} (x/2)**j / j!``, summed in log space so large statistics neither overflow nor
    lose their tail. For odd ``df`` it is the Wilson-Hilferty cube-root normal approximation, within about
    0.01 of the exact value in absolute terms (and proportionally less accurate far in the tail).

    Args:
        x: Test statistic.
        df: Degrees of freedom.

    Returns:
        P-value (survival function).
    """
    if x <= 0:
        return 1.0
    if df <= 0:
        return 0.0

    if df % 2 == 0:
        half = x / 2.0
        log_term = -half  # log of the j = 0 term, exp(-x/2)
        total = 0.0
        for j in range(df // 2):
            total += math.exp(log_term)
            log_term += math.log(half) - math.log(j + 1)
        return min(total, 1.0)

    # Odd df: Wilson-Hilferty normal approximation
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
        p: Probability in (0, 1). A value at or beyond the ends is clamped to a z-score of -10.0 or 10.0
            instead of the infinity it stands for.

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
            per_omic_p: Per-omic p-values {omic: p}. An omic whose results contain none of the
                pathway's genes has no evidence about it and is entered as ``p = 1.0``; it is still one
                of the omics being combined (with Fisher's method that adds two degrees of freedom, with
                ``min_p`` it raises the Bonferroni factor).
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

_DEFAULT_PERMUTATION_SEED = 42


def _validated_random_state(random_state: int | None) -> int | None:
    """``random_state`` as an int (or None for an OS-seeded generator); anything else is a TypeError."""
    if random_state is None:
        return None
    try:
        return operator.index(random_state)
    except TypeError:
        raise TypeError(f"random_state must be an int or None, got {random_state!r}") from None


def _validated_n_permutations(n_permutations: int) -> int:
    """``n_permutations`` as an int >= 1.

    With none, the permutation p-value ``(0 + 1) / (0 + 1)`` would report 1.0 for a test that was never run.
    """
    try:
        count = operator.index(n_permutations)
    except TypeError:
        raise ValueError(f"n_permutations must be an integer >= 1, got {n_permutations!r}") from None
    if count < 1:
        raise ValueError(f"n_permutations must be an integer >= 1, got {n_permutations!r}")
    return count


def _check_p_values(p_values: dict[str, float], label: str) -> None:
    """Raise ValueError naming the first entry of ``p_values`` that is not a p-value in [0, 1].

    A negative p-value would be clamped to the most significant score and a NaN one would silently never
    be significant, so neither is allowed through.
    """
    for key, value in p_values.items():
        number = float(value)
        if not 0.0 <= number <= 1.0:  # False for NaN as well
            raise ValueError(f"{label} must be p-values in [0, 1]; {key!r} has {value!r}")


def _active_module_node_score(p_value: float, threshold: float) -> float:
    """Transformed node score used by :func:`active_module_detection`: ``-log(p) - threshold``."""
    return -math.log(max(p_value, 1e-300)) - threshold


def _module_omic_contributions(
    module_nodes: set[str],
    omic_scores: dict[str, dict[str, float]] | None,
    threshold: float,
) -> dict[str, float]:
    """Module score recomputed from each omic layer's p-values alone.

    Uses the same node transform as the module score itself, and treats a node missing from an
    omic as ``p = 1`` exactly as ``scores`` does, so the values are on the scale of
    ``module_score`` and a layer that does not support the module scores below zero.
    """
    if not omic_scores:
        return {}
    members = sorted(module_nodes)
    return {
        omic: sum(_active_module_node_score(p_values.get(node, 1.0), threshold) for node in members)
        for omic, p_values in omic_scores.items()
    }


def active_module_detection(
    network: dict[str, list[str]],
    scores: dict[str, float],
    alpha: float = 0.05,
    n_permutations: int = 1000,
    omic_scores: dict[str, dict[str, float]] | None = None,
    random_state: int | None = _DEFAULT_PERMUTATION_SEED,
) -> list[dict[str, Any]]:
    """Find active subnetwork modules using node scores from multi-omic integration.

    Uses a greedy seed-and-grow approach:
        1. Rank nodes by score (significance).
        2. Seed modules from the most significant nodes.
        3. Iteratively add neighbours that improve the module score.
        4. Assess significance via permutation testing.

    The module score is the sum of transformed node scores:
        score_i = -log(p_i) - threshold, where threshold = -log(alpha).

    The search is deterministic: ties between equally scored neighbours are broken by node
    identifier, and the permutation test of every module draws from its own generator seeded
    with ``random_state``, so the same input and seed always return the same modules and a
    module's p-value does not depend on which other modules were found.

    Args:
        network: Adjacency list {node: [neighbour1, neighbour2, ...]}.
        scores: Node scores {node: p_value}, each in [0, 1].  Lower p-values indicate
            more significant nodes.
        alpha: Significance threshold for seed selection and score transform; in (0, 1].
        n_permutations: Number of permutations for module significance; an integer >= 1.
        omic_scores: Optional per-omic node p-values ``{omic: {node: p_value}}``, the same
            mapping that :func:`multi_omic_enrichment` takes as ``omic_results``. It does not
            change which modules are found or how they are scored; it only supplies the data for
            ``omic_contributions``, because the combined ``scores`` carry no per-omic information.
        random_state: Seed of the permutation test (default 42, the value that used to be fixed).
            Every module starts a generator from this seed, so equally sized modules are tested
            against the same random node sets; ``None`` seeds each module's generator from the
            operating system, which makes the p-values vary from run to run.

    Returns:
        List of module dictionaries sorted by p-value, each containing:
            module_genes: List of gene/node identifiers in the module.
            module_score: Aggregate score of the module.
            p_value: Permutation-based p-value.
            omic_contributions: ``{omic: score}``, the module score recomputed from that omic's
                p-values alone with the same ``-log(p) - threshold`` transform (a node absent from
                an omic counts as ``p = 1``, as in ``scores``). It shows which omic layers support
                the module: a layer that does not scores below zero. Empty when ``omic_scores`` is
                not given.

    Raises:
        ValueError: If network or scores are empty, ``alpha`` is not in (0, 1], ``n_permutations`` is
            below 1 (no permutation test would be run, yet a p-value of 1.0 would be reported), or a
            p-value in ``scores`` or ``omic_scores`` is outside [0, 1] or NaN.
        TypeError: If ``random_state`` is neither an integer nor None.
    """
    if not network:
        raise ValueError("network must not be empty")
    if not scores:
        raise ValueError("scores must not be empty")
    if not 0.0 < alpha <= 1.0:  # False for NaN as well
        raise ValueError(f"alpha must be in (0, 1], got {alpha!r}")
    n_permutations = _validated_n_permutations(n_permutations)
    random_state = _validated_random_state(random_state)
    _check_p_values(scores, "scores")
    for omic, omic_p_values in (omic_scores or {}).items():
        _check_p_values(omic_p_values, f"omic_scores[{omic!r}]")

    logger.info(
        "Running active module detection: %d nodes, %d scored, alpha=%.4f",
        len(network),
        len(scores),
        alpha,
    )

    threshold = -math.log(max(alpha, 1e-300))

    def node_score(node: str) -> float:
        return _active_module_node_score(scores.get(node, 1.0), threshold)

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
    all_nodes = list(scores.keys())

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

            # sorted: set iteration order varies between processes (string hash randomisation), which
            # would make ties between equal gains, and the float sum below, differ from run to run
            for candidate in sorted(candidates):
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
        module_size = len(module_nodes)
        perm_count = 0
        rng = random.Random(random_state)

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
                "omic_contributions": _module_omic_contributions(module_nodes, omic_scores, threshold),
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
    n_permutations: int = 1000,
    random_state: int | None = _DEFAULT_PERMUTATION_SEED,
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
        n_permutations: Number of permutations behind ``p_value``; an integer >= 1.
        random_state: Seed of the permutation test (default 42, the value that used to be fixed);
            ``None`` seeds it from the operating system, which makes ``p_value`` vary from run to run.

    Returns:
        Dictionary with keys:
            impact_factor: Topology-weighted impact score.
            p_value: Approximate p-value for the impact factor (permutation).
            perturbed_genes: List of significantly perturbed genes.
            pathway_perturbation: Per-gene perturbation details
                {gene: {score, centrality, weighted_score}}.

    Raises:
        ValueError: If inputs are empty or ``n_permutations`` is below 1.
        TypeError: If ``random_state`` is neither an integer nor None.
    """
    if not pathway_graph:
        raise ValueError("pathway_graph must not be empty")
    if not gene_scores:
        raise ValueError("gene_scores must not be empty")
    n_perms = _validated_n_permutations(n_permutations)
    random_state = _validated_random_state(random_state)

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

    # A fixed node order: set iteration order varies between processes (string hash randomisation), and
    # the permutation test below pairs shuffled scores with nodes by position, so its p-value (and the
    # order of float sums) would otherwise change from run to run.
    ordered_nodes = sorted(all_nodes)

    max_degree = max(degree.values()) if degree else 1
    centrality: dict[str, float] = {}
    for node in ordered_nodes:
        centrality[node] = (degree.get(node, 0) + 1) / (max_degree + 1)

    # Compute per-gene perturbation
    perturbation: dict[str, dict[str, float]] = {}
    perturbed_genes: list[str] = []
    sum_weighted = 0.0
    sum_centrality = 0.0

    for gene in ordered_nodes:
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
    rng = random.Random(random_state)
    perm_count = 0
    score_list = list(gene_scores.values())

    for _ in range(n_perms):
        rng.shuffle(score_list)
        perm_weighted = 0.0
        perm_cent = 0.0
        for i, gene in enumerate(ordered_nodes):
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
                (CV <= median CV, the true median: the mean of the two middle
                values when the number of pathways is even).
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

    # Compute -log10(p) per pathway per omic and CV. Pathways are visited in identifier order so that
    # the output order of ties, and the float sum of the CVs, do not depend on set iteration order
    # (which varies between processes).
    heatmap_data: dict[str, dict[str, float]] = {}
    cv_per_pathway: dict[str, float] = {}

    for pathway in sorted(all_pathways):
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
    median_cv = statistics.median(cv_values) if cv_values else 0.0

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
    sorted_pathways = sorted(cv_per_pathway, key=lambda pw: (cv_per_pathway[pw], pw))
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
