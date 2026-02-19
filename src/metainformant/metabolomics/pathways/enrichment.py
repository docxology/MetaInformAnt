"""Metabolite set enrichment analysis.

Implements over-representation analysis (ORA) and quantitative enrichment
scoring for metabolic pathways based on identified metabolite sets.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb, factorial

import numpy as np


@dataclass
class EnrichmentResult:
    """Result of pathway enrichment analysis.

    Attributes:
        pathway_name: Name of the metabolic pathway.
        pathway_size: Number of metabolites in the pathway.
        overlap: Number of query metabolites in the pathway.
        fold_enrichment: Observed / expected ratio.
        p_value: Fisher's exact test p-value.
        metabolites_in_pathway: List of overlapping metabolite names.
    """

    pathway_name: str
    pathway_size: int
    overlap: int
    fold_enrichment: float
    p_value: float
    metabolites_in_pathway: list[str]


def metabolite_set_enrichment(
    query_metabolites: list[str],
    pathway_db: dict[str, list[str]],
    background_size: int | None = None,
) -> list[EnrichmentResult]:
    """Over-representation analysis for metabolite sets.

    Tests whether query metabolites are enriched in each pathway
    using Fisher's exact test (hypergeometric distribution).

    Args:
        query_metabolites: List of identified metabolite names.
        pathway_db: Dict mapping pathway name to list of metabolite names.
        background_size: Total number of metabolites in the reference universe.
            Defaults to the union of all pathway metabolites + query.

    Returns:
        List of EnrichmentResult, sorted by p-value ascending.
    """
    query_set = set(query_metabolites)

    if background_size is None:
        all_mets = set(query_metabolites)
        for mets in pathway_db.values():
            all_mets.update(mets)
        background_size = len(all_mets)

    n_query = len(query_set)
    results: list[EnrichmentResult] = []

    for pathway_name, pathway_mets in pathway_db.items():
        pathway_set = set(pathway_mets)
        overlap_set = query_set & pathway_set
        k = len(overlap_set)
        K = len(pathway_set)

        expected = n_query * K / background_size if background_size > 0 else 0
        fold = k / expected if expected > 0 else 0.0

        # Hypergeometric p-value
        p_val = _hypergeometric_pvalue(k, K, n_query, background_size)

        results.append(
            EnrichmentResult(
                pathway_name=pathway_name,
                pathway_size=K,
                overlap=k,
                fold_enrichment=fold,
                p_value=p_val,
                metabolites_in_pathway=sorted(overlap_set),
            )
        )

    results.sort(key=lambda r: r.p_value)
    return results


def _hypergeometric_pvalue(k: int, K: int, n: int, N: int) -> float:
    """Compute one-tailed hypergeometric p-value P(X >= k).

    Args:
        k: Number of successes in sample.
        K: Total successes in population.
        n: Sample size.
        N: Population size.

    Returns:
        P-value.
    """
    p = 0.0
    for i in range(k, min(K, n) + 1):
        try:
            p += comb(K, i) * comb(N - K, n - i) / comb(N, n)
        except (ValueError, ZeroDivisionError):
            continue
    return min(p, 1.0)


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Apply Benjamini-Hochberg FDR correction to p-values.

    Args:
        p_values: List of raw p-values.

    Returns:
        List of FDR-adjusted q-values in the same order.
    """
    n = len(p_values)
    if n == 0:
        return []

    # Sort by p-value, keeping track of original indices
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    q_values = [0.0] * n

    # Compute adjusted p-values
    prev_q = 1.0
    for rank_minus_1 in range(n - 1, -1, -1):
        orig_idx, pval = indexed[rank_minus_1]
        rank = rank_minus_1 + 1
        q = min(prev_q, pval * n / rank)
        q = min(q, 1.0)
        q_values[orig_idx] = q
        prev_q = q

    return q_values


def enrichment_with_fdr(
    query_metabolites: list[str],
    pathway_db: dict[str, list[str]],
    background_size: int | None = None,
    fdr_threshold: float = 0.05,
) -> list[EnrichmentResult]:
    """Over-representation analysis with FDR correction.

    Performs metabolite set enrichment followed by Benjamini-Hochberg
    multiple testing correction, returning only significant pathways.

    Args:
        query_metabolites: List of identified metabolite names.
        pathway_db: Dict mapping pathway name to metabolite lists.
        background_size: Size of reference universe. Defaults to auto.
        fdr_threshold: Maximum q-value for significance.

    Returns:
        List of significant EnrichmentResult, sorted by q-value.
    """
    all_results = metabolite_set_enrichment(
        query_metabolites, pathway_db, background_size
    )

    if not all_results:
        return []

    p_vals = [r.p_value for r in all_results]
    q_vals = benjamini_hochberg(p_vals)

    significant = []
    for result, q in zip(all_results, q_vals):
        if q <= fdr_threshold:
            significant.append(result)

    return significant


@dataclass
class PathwayActivityScore:
    """Activity score for a metabolic pathway.

    Attributes:
        pathway_name: Name of the metabolic pathway.
        activity_score: Aggregate activity score.
        n_measured: Number of pathway metabolites measured.
        n_total: Total metabolites in pathway.
        contributing_metabolites: Dict of metabolite name to its contribution.
    """

    pathway_name: str
    activity_score: float
    n_measured: int
    n_total: int
    contributing_metabolites: dict[str, float]


def pathway_activity_scoring(
    metabolite_scores: dict[str, float],
    pathway_db: dict[str, list[str]],
    min_coverage: float = 0.1,
) -> list[PathwayActivityScore]:
    """Compute pathway-level activity scores from metabolite-level scores.

    Aggregates per-metabolite scores (e.g., fold changes or t-statistics)
    to pathway level using the mean of mapped metabolites.

    Args:
        metabolite_scores: Dict mapping metabolite name to numeric score.
        pathway_db: Dict mapping pathway name to metabolite lists.
        min_coverage: Minimum fraction of pathway metabolites that must
            be measured to compute a score.

    Returns:
        List of PathwayActivityScore sorted by absolute activity score.
    """
    results: list[PathwayActivityScore] = []

    for pw_name, pw_mets in pathway_db.items():
        contributions: dict[str, float] = {}
        for met in pw_mets:
            if met in metabolite_scores:
                contributions[met] = metabolite_scores[met]

        coverage = len(contributions) / len(pw_mets) if pw_mets else 0
        if coverage < min_coverage:
            continue

        activity = sum(contributions.values()) / len(contributions) if contributions else 0.0

        results.append(
            PathwayActivityScore(
                pathway_name=pw_name,
                activity_score=activity,
                n_measured=len(contributions),
                n_total=len(pw_mets),
                contributing_metabolites=contributions,
            )
        )

    results.sort(key=lambda r: abs(r.activity_score), reverse=True)
    return results
