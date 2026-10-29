"""Metabolite set enrichment analysis.

Implements over-representation analysis (ORA) and quantitative enrichment
scoring for metabolic pathways based on identified metabolite sets.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb


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

    Raises:
        ValueError: If ``background_size`` is smaller than the query or than a pathway, which
            cannot be a universe containing them.
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

    ``X`` is the number of pathway members among ``n`` metabolites drawn without replacement
    from a universe of ``N`` that contains ``K`` pathway members. The upper tail is summed
    exactly with integers and divided once, so the result is correctly rounded and stays
    accurate for very small tails (it only underflows to 0.0 when the true p-value is below
    about 5e-324).

    Args:
        k: Number of successes in sample.
        K: Total successes in population.
        n: Sample size.
        N: Population size.

    Returns:
        P-value, in [0, 1].

    Raises:
        ValueError: If the arguments do not describe a draw from a population (``K`` or ``n``
            negative or larger than ``N``). Such a call used to return 0.0, the most
            significant p-value possible, for what is really an invalid ``background_size``.
    """
    if N < 0 or not 0 <= K <= N or not 0 <= n <= N:
        raise ValueError(
            f"inconsistent hypergeometric parameters: population N={N}, successes K={K}, draws n={n} "
            "(need 0 <= K <= N and 0 <= n <= N). In metabolite_set_enrichment N is background_size, "
            "K the pathway size and n the number of distinct query metabolites, so background_size "
            "must be at least as large as the query and as every pathway"
        )
    # X takes values from max(0, n + K - N) to min(K, n); P(X >= k) sums that support from k up
    lowest = max(k, 0, n + K - N)
    highest = min(K, n)
    tail = sum(comb(K, i) * comb(N - K, n - i) for i in range(lowest, highest + 1))
    return min(tail / comb(N, n), 1.0)


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
    all_results = metabolite_set_enrichment(query_metabolites, pathway_db, background_size)

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
            be measured to compute a score. A pathway with no measured
            metabolite is never scored, whatever the value: there is nothing
            to average, and a score of 0.0 would be made up.

    Returns:
        List of PathwayActivityScore sorted by absolute activity score.
    """
    results: list[PathwayActivityScore] = []

    for pw_name, pw_mets in pathway_db.items():
        contributions: dict[str, float] = {}
        for met in pw_mets:
            if met in metabolite_scores:
                contributions[met] = metabolite_scores[met]

        if not contributions:
            continue
        coverage = len(contributions) / len(pw_mets)
        if coverage < min_coverage:
            continue

        activity = sum(contributions.values()) / len(contributions)

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
