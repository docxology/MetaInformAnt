"""Pathway enrichment analysis: ORA, GSEA, and pathway networks.

Provides over-representation analysis using the hypergeometric test,
Gene Set Enrichment Analysis with running enrichment score computation,
pathway similarity network construction, and cross-condition enrichment
comparison.

All statistical computations use pure Python implementations with no
external dependencies beyond the standard library.
"""

from __future__ import annotations

import math
import random
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _log_factorial(n: int) -> float:
    """Compute log(n!) using Stirling's approximation for large n."""
    if n <= 1:
        return 0.0
    if n <= 20:
        result = 0.0
        for i in range(2, n + 1):
            result += math.log(i)
        return result
    # Stirling's approximation
    return n * math.log(n) - n + 0.5 * math.log(2 * math.pi * n)


def _log_comb(n: int, k: int) -> float:
    """Compute log(C(n, k))."""
    if k < 0 or k > n:
        return float("-inf")
    if k == 0 or k == n:
        return 0.0
    return _log_factorial(n) - _log_factorial(k) - _log_factorial(n - k)


def _hypergeometric_sf(k: int, M: int, n: int, N: int) -> float:
    """Survival function for hypergeometric distribution.

    P(X >= k) where X ~ Hypergeometric(M, n, N).

    Args:
        k: Observed overlap (successes in sample).
        M: Population size.
        n: Number of successes in population.
        N: Sample size.
    """
    if k <= 0:
        return 1.0

    p_value = 0.0
    max_k = min(n, N)
    for x in range(k, max_k + 1):
        log_p = _log_comb(n, x) + _log_comb(M - n, N - x) - _log_comb(M, N)
        if log_p > -700:
            p_value += math.exp(log_p)

    return min(1.0, p_value)


def _fdr_correction(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR correction.

    Args:
        p_values: List of raw p-values.

    Returns:
        List of adjusted p-values.
    """
    n = len(p_values)
    if n == 0:
        return []

    # Sort by p-value
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0.0] * n

    prev_adj = 1.0
    for rank_idx in range(n - 1, -1, -1):
        orig_idx, p = indexed[rank_idx]
        rank = rank_idx + 1
        adj = min(p * n / rank, prev_adj)
        adj = min(adj, 1.0)
        adjusted[orig_idx] = adj
        prev_adj = adj

    return adjusted


def _bonferroni_correction(p_values: list[float]) -> list[float]:
    """Bonferroni correction."""
    n = len(p_values)
    return [min(1.0, p * n) for p in p_values]


# ---------------------------------------------------------------------------
# Over-Representation Analysis
# ---------------------------------------------------------------------------


def over_representation_analysis(
    gene_list: list[str],
    gene_sets: dict,
    background: list[str] | None = None,
    correction: str = "fdr_bh",
) -> list[dict]:
    """Over-representation analysis using Fisher's exact / hypergeometric test.

    Tests whether genes in the input list are enriched in any gene set
    relative to a background gene universe.

    Args:
        gene_list: List of genes of interest (e.g. differentially expressed).
        gene_sets: Dictionary mapping term ID/name to set or list of genes
            in that pathway/term.
        background: Optional background gene universe. If ``None``, the
            union of all genes in gene_sets is used.
        correction: Multiple testing correction: ``"fdr_bh"`` (default),
            ``"bonferroni"``, or ``"none"``.

    Returns:
        List of result dicts sorted by adjusted p-value, each with:
            - term_id: Gene set identifier.
            - term_name: Gene set name (same as term_id if not provided).
            - p_value: Raw hypergeometric p-value.
            - adjusted_p: Corrected p-value.
            - odds_ratio: Odds ratio.
            - n_genes: Size of the gene set.
            - n_overlap: Number of overlapping genes.
            - overlap_genes: List of overlapping gene names.
    """
    query_set = set(gene_list)

    if background is None:
        bg_set: set[str] = set()
        for gs in gene_sets.values():
            bg_set.update(gs)
        bg_set.update(query_set)
    else:
        bg_set = set(background)

    M = len(bg_set)  # Population size
    N = len(query_set.intersection(bg_set))  # Sample size

    results: list[dict] = []

    for term_id, genes in gene_sets.items():
        term_genes = set(genes).intersection(bg_set)
        n = len(term_genes)  # Successes in population
        if n == 0:
            continue

        overlap = query_set.intersection(term_genes)
        k = len(overlap)  # Successes in sample

        if k == 0:
            p_value = 1.0
        else:
            p_value = _hypergeometric_sf(k, M, n, N)

        # Odds ratio
        a = k
        b = N - k
        c = n - k
        d = M - n - N + k
        if b * c > 0:
            odds_ratio = (a * d) / (b * c)
        elif a > 0:
            odds_ratio = float("inf")
        else:
            odds_ratio = 0.0

        results.append(
            {
                "term_id": term_id,
                "term_name": term_id,
                "p_value": p_value,
                "adjusted_p": p_value,  # Will be updated
                "odds_ratio": odds_ratio,
                "n_genes": n,
                "n_overlap": k,
                "overlap_genes": sorted(overlap),
            }
        )

    # Multiple testing correction
    if results:
        raw_ps = [r["p_value"] for r in results]
        if correction == "fdr_bh":
            adj_ps = _fdr_correction(raw_ps)
        elif correction == "bonferroni":
            adj_ps = _bonferroni_correction(raw_ps)
        else:
            adj_ps = raw_ps

        for i, r in enumerate(results):
            r["adjusted_p"] = adj_ps[i]

    results.sort(key=lambda r: r["adjusted_p"])
    logger.info(
        "ORA complete: %d gene sets tested, %d with overlap",
        len(gene_sets),
        sum(1 for r in results if r["n_overlap"] > 0),
    )

    return results


# ---------------------------------------------------------------------------
# Enrichment Score Computation
# ---------------------------------------------------------------------------


def compute_enrichment_score(
    ranked_list: list[str],
    gene_set: set,
    weighted: bool = True,
    weights: list[float] | None = None,
) -> dict:
    """Compute running enrichment score for GSEA.

    Args:
        ranked_list: Genes ordered by rank (e.g. by t-statistic or fold change).
        gene_set: Set of genes in the pathway.
        weighted: If ``True``, weight hits by their rank metric (|weight|).
        weights: Optional per-gene weights aligned with ranked_list.
            If ``None`` and weighted=True, uses uniform weights.

    Returns:
        Dictionary with keys:
            - es: Enrichment score (maximum deviation from zero).
            - running_es: List of running enrichment scores.
            - hit_indices: Indices of gene set members in ranked list.
            - leading_edge_index: Index of the peak enrichment score.
    """
    n = len(ranked_list)
    gene_set_lower = {g.lower() if isinstance(g, str) else g for g in gene_set}

    if weights is None:
        weights = [1.0] * n

    hits = []
    for i, gene in enumerate(ranked_list):
        gene_key = gene.lower() if isinstance(gene, str) else gene
        if gene_key in gene_set_lower:
            hits.append(i)

    n_hits = len(hits)
    n_miss = n - n_hits

    if n_hits == 0 or n_miss == 0:
        return {
            "es": 0.0,
            "running_es": [0.0] * n,
            "hit_indices": hits,
            "leading_edge_index": 0,
        }

    # Weighted hit scores
    hit_set = set(hits)
    if weighted:
        hit_weights = [abs(weights[i]) for i in hits]
        total_hit_weight = sum(hit_weights)
    else:
        total_hit_weight = float(n_hits)

    miss_penalty = 1.0 / n_miss

    running = []
    score = 0.0
    max_score = 0.0
    min_score = 0.0
    max_idx = 0
    min_idx = 0

    hit_counter = 0
    for i in range(n):
        if i in hit_set:
            if weighted:
                score += abs(weights[i]) / total_hit_weight if total_hit_weight > 0 else 0.0
            else:
                score += 1.0 / total_hit_weight if total_hit_weight > 0 else 0.0
            hit_counter += 1
        else:
            score -= miss_penalty

        running.append(score)

        if score > max_score:
            max_score = score
            max_idx = i
        if score < min_score:
            min_score = score
            min_idx = i

    # ES is the maximum deviation
    if abs(max_score) >= abs(min_score):
        es = max_score
        leading_idx = max_idx
    else:
        es = min_score
        leading_idx = min_idx

    return {
        "es": es,
        "running_es": running,
        "hit_indices": hits,
        "leading_edge_index": leading_idx,
    }


# ---------------------------------------------------------------------------
# GSEA
# ---------------------------------------------------------------------------


def gsea(
    ranked_genes: list[tuple],
    gene_sets: dict,
    n_permutations: int = 1000,
    min_size: int = 15,
    max_size: int = 500,
) -> list[dict]:
    """Gene Set Enrichment Analysis.

    Computes enrichment scores for gene sets against a ranked gene list,
    normalises by permutation, and computes FDR q-values.

    Args:
        ranked_genes: List of ``(gene_name, rank_metric)`` tuples, sorted
            by rank metric (e.g. signed -log10 p-value or t-statistic).
        gene_sets: Dictionary mapping term ID to gene set (list or set).
        n_permutations: Number of permutations for p-value estimation.
        min_size: Minimum gene set size to test.
        max_size: Maximum gene set size to test.

    Returns:
        List of result dicts sorted by FDR, each with:
            - term_id, term_name, es, nes, p_value, fdr,
              leading_edge, leading_edge_genes.
    """
    ranked_list = [g[0] for g in ranked_genes]
    rank_weights = [float(g[1]) for g in ranked_genes]
    n = len(ranked_list)

    results: list[dict] = []

    for term_id, genes in gene_sets.items():
        gs = set(genes)
        gs_in_list = gs.intersection(set(ranked_list))
        if len(gs_in_list) < min_size or len(gs_in_list) > max_size:
            continue

        # Observed enrichment score
        obs = compute_enrichment_score(ranked_list, gs, weighted=True, weights=rank_weights)
        obs_es = obs["es"]

        # Permutation null distribution
        null_es: list[float] = []
        for _ in range(n_permutations):
            perm_weights = list(rank_weights)
            random.shuffle(perm_weights)
            perm_result = compute_enrichment_score(ranked_list, gs, weighted=True, weights=perm_weights)
            null_es.append(perm_result["es"])

        # Normalise: NES = ES / mean(|null_ES|) for same sign
        pos_null = [e for e in null_es if e > 0]
        neg_null = [e for e in null_es if e < 0]

        if obs_es >= 0:
            mean_null = sum(pos_null) / len(pos_null) if pos_null else 1.0
        else:
            mean_null = -sum(neg_null) / len(neg_null) if neg_null else 1.0

        nes = obs_es / mean_null if mean_null != 0 else 0.0

        # P-value
        if obs_es >= 0:
            p_value = (sum(1 for e in null_es if e >= obs_es) + 1) / (n_permutations + 1)
        else:
            p_value = (sum(1 for e in null_es if e <= obs_es) + 1) / (n_permutations + 1)

        # Leading edge genes
        le_idx = obs["leading_edge_index"]
        if obs_es >= 0:
            le_genes = [ranked_list[i] for i in obs["hit_indices"] if i <= le_idx]
        else:
            le_genes = [ranked_list[i] for i in obs["hit_indices"] if i >= le_idx]

        results.append(
            {
                "term_id": term_id,
                "term_name": term_id,
                "es": obs_es,
                "nes": nes,
                "p_value": p_value,
                "fdr": p_value,  # Will be corrected below
                "leading_edge": len(le_genes),
                "leading_edge_genes": le_genes,
            }
        )

    # FDR correction
    if results:
        raw_ps = [r["p_value"] for r in results]
        adj_ps = _fdr_correction(raw_ps)
        for i, r in enumerate(results):
            r["fdr"] = adj_ps[i]

    results.sort(key=lambda r: r["fdr"])

    logger.info(
        "GSEA complete: %d gene sets tested, %d significant (FDR < 0.25)",
        len(results),
        sum(1 for r in results if r["fdr"] < 0.25),
    )

    return results


# ---------------------------------------------------------------------------
# Pathway network
# ---------------------------------------------------------------------------


def pathway_network(
    enrichment_results: list[dict],
    gene_sets: dict,
    similarity_threshold: float = 0.3,
) -> dict:
    """Build a pathway similarity network based on Jaccard coefficient.

    Connects pathways that share a sufficient fraction of their genes.

    Args:
        enrichment_results: List of enrichment result dicts (from ORA or GSEA)
            with at least ``term_id`` key.
        gene_sets: Dictionary mapping term ID to gene set.
        similarity_threshold: Minimum Jaccard similarity for an edge.

    Returns:
        Dictionary with keys:
            - nodes: List of node dicts with term_id and metadata.
            - edges: List of edge dicts with source, target, jaccard.
            - clusters: List of connected component clusters.
    """
    term_ids = [r["term_id"] for r in enrichment_results if r["term_id"] in gene_sets]

    nodes = []
    for r in enrichment_results:
        if r["term_id"] in gene_sets:
            nodes.append(
                {
                    "term_id": r["term_id"],
                    "p_value": r.get("p_value", 1.0),
                    "adjusted_p": r.get("adjusted_p", r.get("fdr", 1.0)),
                    "n_genes": len(gene_sets.get(r["term_id"], [])),
                }
            )

    edges = []
    for i in range(len(term_ids)):
        set_a = set(gene_sets.get(term_ids[i], []))
        for j in range(i + 1, len(term_ids)):
            set_b = set(gene_sets.get(term_ids[j], []))
            intersection = len(set_a.intersection(set_b))
            union = len(set_a.union(set_b))
            jaccard = intersection / union if union > 0 else 0.0

            if jaccard >= similarity_threshold:
                edges.append(
                    {
                        "source": term_ids[i],
                        "target": term_ids[j],
                        "jaccard": jaccard,
                        "n_shared": intersection,
                    }
                )

    # Simple connected components
    parent: dict[str, str] = {tid: tid for tid in term_ids}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for edge in edges:
        union(edge["source"], edge["target"])

    clusters_map: dict[str, list[str]] = {}
    for tid in term_ids:
        root = find(tid)
        clusters_map.setdefault(root, []).append(tid)

    clusters = list(clusters_map.values())

    logger.info(
        "Pathway network: %d nodes, %d edges, %d clusters",
        len(nodes),
        len(edges),
        len(clusters),
    )

    return {
        "nodes": nodes,
        "edges": edges,
        "clusters": clusters,
    }


# ---------------------------------------------------------------------------
# Compare enrichments
# ---------------------------------------------------------------------------


def compare_enrichments(
    results_a: list[dict],
    results_b: list[dict],
) -> dict:
    """Compare enrichment results between two conditions.

    Args:
        results_a: Enrichment results from condition A.
        results_b: Enrichment results from condition B.

    Returns:
        Dictionary with keys:
            - shared_terms: List of terms significant in both conditions.
            - unique_a: List of terms significant only in A.
            - unique_b: List of terms significant only in B.
            - concordance: Fraction of shared terms with same direction.
            - comparison_table: List of dicts for all terms with both results.
    """
    terms_a = {r["term_id"]: r for r in results_a}
    terms_b = {r["term_id"]: r for r in results_b}

    all_terms = set(terms_a.keys()).union(set(terms_b.keys()))

    sig_threshold = 0.05
    sig_a = {
        tid for tid, r in terms_a.items() if r.get("adjusted_p", r.get("fdr", r.get("p_value", 1.0))) < sig_threshold
    }
    sig_b = {
        tid for tid, r in terms_b.items() if r.get("adjusted_p", r.get("fdr", r.get("p_value", 1.0))) < sig_threshold
    }

    shared = sig_a.intersection(sig_b)
    unique_a = sig_a - sig_b
    unique_b = sig_b - sig_a

    # Concordance: for shared terms, check if ES/beta have same sign
    concordant = 0
    for tid in shared:
        es_a = terms_a[tid].get("es", terms_a[tid].get("odds_ratio", 0))
        es_b = terms_b[tid].get("es", terms_b[tid].get("odds_ratio", 0))
        if (es_a > 0 and es_b > 0) or (es_a < 0 and es_b < 0):
            concordant += 1

    concordance = concordant / len(shared) if shared else 0.0

    # Comparison table
    comparison_table = []
    for tid in all_terms:
        entry: dict[str, Any] = {"term_id": tid}
        if tid in terms_a:
            entry["p_a"] = terms_a[tid].get("adjusted_p", terms_a[tid].get("fdr", terms_a[tid].get("p_value", 1.0)))
            entry["es_a"] = terms_a[tid].get("es", terms_a[tid].get("odds_ratio", 0))
        if tid in terms_b:
            entry["p_b"] = terms_b[tid].get("adjusted_p", terms_b[tid].get("fdr", terms_b[tid].get("p_value", 1.0)))
            entry["es_b"] = terms_b[tid].get("es", terms_b[tid].get("odds_ratio", 0))

        in_a = tid in sig_a
        in_b = tid in sig_b
        if in_a and in_b:
            entry["category"] = "shared"
        elif in_a:
            entry["category"] = "unique_a"
        elif in_b:
            entry["category"] = "unique_b"
        else:
            entry["category"] = "not_significant"

        comparison_table.append(entry)

    logger.info(
        "Enrichment comparison: %d shared, %d unique_a, %d unique_b, concordance=%.3f",
        len(shared),
        len(unique_a),
        len(unique_b),
        concordance,
    )

    return {
        "shared_terms": sorted(shared),
        "unique_a": sorted(unique_a),
        "unique_b": sorted(unique_b),
        "concordance": concordance,
        "comparison_table": comparison_table,
    }
