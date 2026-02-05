"""Protein domain-based classification and clustering utilities.

This module provides tools for classifying proteins into families based on
their domain composition, computing domain architecture similarity, clustering
proteins by shared domain content, and performing domain enrichment analysis.
All methods use pure Python implementations.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any, Dict, List, Optional, Set, Tuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Domain-to-family classification rules
# ---------------------------------------------------------------------------

# Maps domain combinations (frozensets of domain names) to protein families.
# Order-independent matching: presence of domain names is what matters.
_FAMILY_RULES: List[Dict[str, Any]] = [
    {
        "family": "Protein kinase",
        "required_domains": {"Protein kinase domain"},
        "optional_domains": {"SH2", "SH3", "PH domain"},
        "description": "Enzymes that transfer phosphate groups to substrates",
        "go_terms": ["GO:0004672", "GO:0006468"],
    },
    {
        "family": "Serine protease",
        "required_domains": {"Trypsin-like serine protease"},
        "optional_domains": {"Kringle", "EGF-like domain"},
        "description": "Proteases with serine in the catalytic triad",
        "go_terms": ["GO:0004252", "GO:0006508"],
    },
    {
        "family": "Zinc finger transcription factor",
        "required_domains": {"C2H2 zinc finger"},
        "optional_domains": {"KRAB domain", "SCAN domain"},
        "description": "Transcription factors with C2H2 zinc finger DNA binding domains",
        "go_terms": ["GO:0003700", "GO:0006355"],
    },
    {
        "family": "Nuclear receptor",
        "required_domains": {"C4"},
        "optional_domains": {"Ligand-binding domain"},
        "description": "Ligand-activated transcription factors with C4 zinc fingers",
        "go_terms": ["GO:0004879", "GO:0006355"],
    },
    {
        "family": "Small GTPase",
        "required_domains": {"Ras family GTPase"},
        "optional_domains": {"PH domain", "RhoGAP"},
        "description": "GTP-binding molecular switches of the Ras superfamily",
        "go_terms": ["GO:0003924", "GO:0007264"],
    },
    {
        "family": "RNA-binding protein",
        "required_domains": {"RNA recognition motif (RRM)"},
        "optional_domains": {"KH domain", "DEAD-box helicase"},
        "description": "Proteins that bind RNA through RRM domains",
        "go_terms": ["GO:0003723", "GO:0000398"],
    },
    {
        "family": "WD40 repeat protein",
        "required_domains": {"WD40 repeat"},
        "optional_domains": set(),
        "description": "Beta-propeller scaffold proteins with WD40 repeats",
        "go_terms": ["GO:0005515"],
    },
    {
        "family": "Homeodomain protein",
        "required_domains": {"Homeodomain"},
        "optional_domains": {"POU domain", "LIM domain"},
        "description": "Transcription factors with homeodomain DNA binding",
        "go_terms": ["GO:0003700", "GO:0006355"],
    },
    {
        "family": "DEAD-box helicase",
        "required_domains": {"Helicase C-terminal domain"},
        "optional_domains": {"DEAD-box helicase N-terminal", "RNA recognition motif (RRM)"},
        "description": "ATP-dependent RNA helicases",
        "go_terms": ["GO:0004386", "GO:0008026"],
    },
]

# Broader classification by single-domain presence
_SINGLE_DOMAIN_FAMILIES: Dict[str, Dict[str, Any]] = {
    "Protein kinase domain": {
        "superfamily": "Transferase",
        "class": "Kinase",
        "ec_number": "2.7.x.x",
    },
    "Trypsin-like serine protease": {
        "superfamily": "Hydrolase",
        "class": "Protease",
        "ec_number": "3.4.21.x",
    },
    "Ras family GTPase": {
        "superfamily": "Hydrolase",
        "class": "GTPase",
        "ec_number": "3.6.5.x",
    },
    "RNA recognition motif (RRM)": {
        "superfamily": "RNA-binding",
        "class": "RNA-binding protein",
        "ec_number": None,
    },
    "WD40 repeat": {
        "superfamily": "Repeat",
        "class": "Beta-propeller",
        "ec_number": None,
    },
    "Homeodomain": {
        "superfamily": "DNA-binding",
        "class": "Transcription factor",
        "ec_number": None,
    },
    "Helicase C-terminal domain": {
        "superfamily": "Hydrolase",
        "class": "Helicase",
        "ec_number": "3.6.4.x",
    },
}


def classify_protein_family(
    domains: List[Dict[str, Any]],
    sequence: str,
) -> Dict[str, Any]:
    """Classify a protein into a family based on its domain composition.

    Matches the observed domain set against a curated rule base of domain
    combinations that define known protein families (kinases, proteases,
    transcription factors, etc.). Also provides broader superfamily and
    class assignments from individual domains.

    Args:
        domains: List of domain hit dicts, each containing at minimum
            a ``name`` key with the domain name string.
        sequence: Protein amino acid sequence (used for additional
            sequence-based features when domain hits are ambiguous).

    Returns:
        Dict with keys:
            - family: best-matching protein family name, or "Unclassified"
            - confidence: classification confidence (0.0 to 1.0)
            - matched_rule: the rule that matched (or None)
            - superfamily: broader classification
            - enzyme_class: enzyme classification if applicable
            - ec_number: EC number if applicable
            - go_terms: associated GO terms
            - domain_names: list of domain names used for classification
            - alternative_families: other families that partially match
    """
    if not domains:
        return {
            "family": "Unclassified",
            "confidence": 0.0,
            "matched_rule": None,
            "superfamily": "Unknown",
            "enzyme_class": None,
            "ec_number": None,
            "go_terms": [],
            "domain_names": [],
            "alternative_families": [],
        }

    # Extract domain names
    domain_names: Set[str] = {d.get("name", "") for d in domains if d.get("name")}

    # Try to match family rules (best match = most required domains satisfied)
    best_match: Optional[Dict[str, Any]] = None
    best_score = 0.0
    alternatives: List[Dict[str, Any]] = []

    for rule in _FAMILY_RULES:
        required = rule["required_domains"]
        optional = rule.get("optional_domains", set())

        # Count how many required domains are present
        required_hits = required & domain_names
        optional_hits = optional & domain_names

        if not required_hits:
            continue

        # Score: fraction of required domains found + bonus for optional
        req_score = len(required_hits) / len(required)
        opt_score = len(optional_hits) / max(len(optional), 1) * 0.2
        total_score = req_score + opt_score

        if total_score > best_score:
            if best_match is not None:
                alternatives.append(
                    {
                        "family": best_match["family"],
                        "confidence": round(best_score, 4),
                    }
                )
            best_match = rule
            best_score = total_score
        elif total_score > 0.5:
            alternatives.append(
                {
                    "family": rule["family"],
                    "confidence": round(total_score, 4),
                }
            )

    # Determine superfamily / class from individual domains
    superfamily = "Unknown"
    enzyme_class = None
    ec_number = None

    for dname in domain_names:
        if dname in _SINGLE_DOMAIN_FAMILIES:
            info = _SINGLE_DOMAIN_FAMILIES[dname]
            superfamily = info["superfamily"]
            enzyme_class = info["class"]
            ec_number = info["ec_number"]
            break

    if best_match is not None:
        family = best_match["family"]
        confidence = round(min(best_score, 1.0), 4)
        go_terms = best_match.get("go_terms", [])
    else:
        family = "Unclassified"
        confidence = 0.0
        go_terms = []

    return {
        "family": family,
        "confidence": confidence,
        "matched_rule": best_match,
        "superfamily": superfamily,
        "enzyme_class": enzyme_class,
        "ec_number": ec_number,
        "go_terms": go_terms,
        "domain_names": sorted(domain_names),
        "alternative_families": alternatives,
    }


# ---------------------------------------------------------------------------
# Domain architecture similarity
# ---------------------------------------------------------------------------


def compute_domain_similarity(
    domains_a: List[Dict[str, Any]],
    domains_b: List[Dict[str, Any]],
) -> float:
    """Compute Jaccard-like similarity between two domain architectures.

    Compares both domain composition (which domains are present) and domain
    order (the sequential arrangement of domains along the protein). The
    final score is a weighted combination of compositional Jaccard similarity
    (70%) and order similarity via longest common subsequence (30%).

    Args:
        domains_a: Domain list for protein A (each dict must have ``name``).
        domains_b: Domain list for protein B (each dict must have ``name``).

    Returns:
        Similarity score between 0.0 (completely different) and 1.0 (identical).
    """
    if not domains_a and not domains_b:
        return 1.0
    if not domains_a or not domains_b:
        return 0.0

    # Extract ordered domain name lists
    names_a = [d.get("name", d.get("domain_id", "")) for d in sorted(domains_a, key=lambda x: x.get("start", 0))]
    names_b = [d.get("name", d.get("domain_id", "")) for d in sorted(domains_b, key=lambda x: x.get("start", 0))]

    set_a = set(names_a)
    set_b = set(names_b)

    # Jaccard similarity (composition)
    intersection = set_a & set_b
    union = set_a | set_b
    jaccard = len(intersection) / len(union) if union else 1.0

    # Order similarity via longest common subsequence (LCS)
    lcs_len = _longest_common_subsequence_length(names_a, names_b)
    max_len = max(len(names_a), len(names_b))
    order_sim = lcs_len / max_len if max_len > 0 else 1.0

    # Weighted combination
    similarity = 0.7 * jaccard + 0.3 * order_sim

    return round(similarity, 4)


def _longest_common_subsequence_length(seq_a: List[str], seq_b: List[str]) -> int:
    """Compute length of the longest common subsequence.

    Args:
        seq_a: First sequence of strings.
        seq_b: Second sequence of strings.

    Returns:
        Length of the LCS.
    """
    m, n = len(seq_a), len(seq_b)
    if m == 0 or n == 0:
        return 0

    # Space-optimized DP (two rows)
    prev = [0] * (n + 1)
    curr = [0] * (n + 1)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq_a[i - 1] == seq_b[j - 1]:
                curr[j] = prev[j - 1] + 1
            else:
                curr[j] = max(prev[j], curr[j - 1])
        prev, curr = curr, [0] * (n + 1)

    return prev[n]


# ---------------------------------------------------------------------------
# Protein clustering by domain architecture
# ---------------------------------------------------------------------------


def cluster_by_domain(
    proteins: Dict[str, List[Dict[str, Any]]],
    method: str = "jaccard",
    threshold: float = 0.5,
) -> Dict[str, Any]:
    """Cluster proteins by domain architecture similarity.

    Performs single-linkage agglomerative clustering based on pairwise domain
    architecture similarity. Proteins are grouped into clusters where every
    member shares at least ``threshold`` similarity with at least one other
    member.

    Args:
        proteins: Dict mapping protein identifiers to their domain hit lists.
            Each domain hit dict must contain a ``name`` key.
        method: Similarity method. Currently supports ``"jaccard"`` (domain
            composition Jaccard index with order bonus).
        threshold: Minimum similarity to merge two proteins into the same
            cluster (0.0 to 1.0).

    Returns:
        Dict with keys:
            - clusters: dict mapping cluster_id (int) to list of protein IDs
            - cluster_count: number of clusters
            - singleton_count: number of single-protein clusters
            - similarity_matrix: dict of (protein_a, protein_b) -> similarity
            - method: the method used
            - threshold: the threshold used

    Raises:
        ValueError: If method is not supported.
    """
    if method not in ("jaccard",):
        raise ValueError(f"Unsupported clustering method: {method}. Use 'jaccard'.")

    protein_ids = sorted(proteins.keys())
    n = len(protein_ids)

    if n == 0:
        return {
            "clusters": {},
            "cluster_count": 0,
            "singleton_count": 0,
            "similarity_matrix": {},
            "method": method,
            "threshold": threshold,
        }

    # Compute pairwise similarity matrix
    sim_matrix: Dict[Tuple[str, str], float] = {}
    for i in range(n):
        for j in range(i + 1, n):
            pid_a = protein_ids[i]
            pid_b = protein_ids[j]
            sim = compute_domain_similarity(proteins[pid_a], proteins[pid_b])
            sim_matrix[(pid_a, pid_b)] = sim
            sim_matrix[(pid_b, pid_a)] = sim

    # Single-linkage clustering using union-find
    parent: Dict[str, str] = {pid: pid for pid in protein_ids}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: str, y: str) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    # Merge proteins above threshold
    for i in range(n):
        for j in range(i + 1, n):
            pid_a = protein_ids[i]
            pid_b = protein_ids[j]
            if sim_matrix.get((pid_a, pid_b), 0.0) >= threshold:
                union(pid_a, pid_b)

    # Build cluster assignments
    cluster_map: Dict[str, List[str]] = {}
    for pid in protein_ids:
        root = find(pid)
        if root not in cluster_map:
            cluster_map[root] = []
        cluster_map[root].append(pid)

    # Re-number clusters starting from 0
    clusters: Dict[int, List[str]] = {}
    for idx, (_, members) in enumerate(sorted(cluster_map.items())):
        clusters[idx] = sorted(members)

    singleton_count = sum(1 for members in clusters.values() if len(members) == 1)

    logger.debug(
        "Clustered %d proteins into %d clusters (%d singletons)",
        n,
        len(clusters),
        singleton_count,
    )

    return {
        "clusters": clusters,
        "cluster_count": len(clusters),
        "singleton_count": singleton_count,
        "similarity_matrix": sim_matrix,
        "method": method,
        "threshold": threshold,
    }


# ---------------------------------------------------------------------------
# Domain enrichment analysis
# ---------------------------------------------------------------------------


def domain_enrichment(
    protein_set: List[str],
    background: List[str],
    domain_annotations: Dict[str, List[str]],
) -> List[Dict[str, Any]]:
    """Test enrichment of specific domains in a protein set versus background.

    For each domain observed in the protein set, performs a one-sided
    Fisher's exact test (computed via the hypergeometric distribution)
    to determine whether the domain is over-represented compared to the
    background set. Results are sorted by p-value ascending.

    Args:
        protein_set: List of protein identifiers in the foreground set.
        background: List of protein identifiers in the background set.
            Should include the foreground proteins.
        domain_annotations: Dict mapping protein identifiers to lists
            of domain names present in that protein.

    Returns:
        List of dicts sorted by p-value, each with:
            - domain: domain name
            - count_in_set: number of proteins in set with this domain
            - count_in_background: number in background with this domain
            - set_size: total foreground set size
            - background_size: total background size
            - fold_enrichment: ratio of observed to expected frequency
            - p_value: Fisher's exact test p-value (one-sided, over-representation)
            - significant: True if p_value < 0.05

    Raises:
        ValueError: If protein_set is empty or background is empty.
    """
    if not protein_set:
        raise ValueError("protein_set must be non-empty")
    if not background:
        raise ValueError("background must be non-empty")

    set_ids = set(protein_set)
    bg_ids = set(background)
    set_size = len(set_ids)
    bg_size = len(bg_ids)

    # Count domain occurrences in set and background
    set_domain_counts: Counter[str] = Counter()
    bg_domain_counts: Counter[str] = Counter()

    for pid in bg_ids:
        domains = domain_annotations.get(pid, [])
        for domain in domains:
            bg_domain_counts[domain] += 1
            if pid in set_ids:
                set_domain_counts[domain] += 1

    # Also count domains from set proteins not in background
    for pid in set_ids:
        if pid not in bg_ids:
            domains = domain_annotations.get(pid, [])
            for domain in domains:
                set_domain_counts[domain] += 1
                bg_domain_counts[domain] += 1
                bg_size += 1  # conceptually expand background

    # Test each domain
    all_domains = set(set_domain_counts.keys()) | set(bg_domain_counts.keys())
    results: List[Dict[str, Any]] = []

    for domain in sorted(all_domains):
        k = set_domain_counts.get(domain, 0)  # successes in sample
        K = bg_domain_counts.get(domain, 0)  # successes in population
        n = set_size  # sample size
        N = bg_size  # population size

        if k == 0:
            continue

        # Fold enrichment
        expected = (K / N) * n if N > 0 else 0.0
        fold_enrichment = k / expected if expected > 0 else float("inf")

        # One-sided Fisher's exact p-value via hypergeometric CDF
        # P(X >= k) = sum_{i=k}^{min(n,K)} hyper_pmf(i, N, K, n)
        p_value = _hypergeometric_sf(k, N, K, n)

        results.append(
            {
                "domain": domain,
                "count_in_set": k,
                "count_in_background": K,
                "set_size": set_size,
                "background_size": bg_size,
                "fold_enrichment": round(fold_enrichment, 4),
                "p_value": round(p_value, 6),
                "significant": p_value < 0.05,
            }
        )

    # Sort by p-value ascending
    results.sort(key=lambda r: r["p_value"])

    logger.debug(
        "Domain enrichment: tested %d domains, %d significant",
        len(results),
        sum(1 for r in results if r["significant"]),
    )

    return results


def _hypergeometric_sf(k: int, N: int, K: int, n: int) -> float:
    """Compute survival function P(X >= k) for hypergeometric distribution.

    Uses log-space computation to avoid overflow for large values.

    Args:
        k: Number of observed successes.
        N: Population size.
        K: Number of success states in population.
        n: Number of draws (sample size).

    Returns:
        P(X >= k), the probability of observing k or more successes.
    """
    if k <= 0:
        return 1.0
    if k > min(n, K):
        return 0.0

    # P(X >= k) = 1 - P(X <= k-1) = sum_{i=k}^{min(n,K)} pmf(i)
    upper = min(n, K)
    total_prob = 0.0

    for i in range(k, upper + 1):
        log_pmf = _log_hypergeometric_pmf(i, N, K, n)
        if log_pmf > -700:  # avoid underflow
            total_prob += math.exp(log_pmf)

    return min(total_prob, 1.0)


def _log_hypergeometric_pmf(k: int, N: int, K: int, n: int) -> float:
    """Compute log of hypergeometric PMF.

    PMF = C(K, k) * C(N-K, n-k) / C(N, n)

    Args:
        k: Number of successes in draw.
        N: Population size.
        K: Number of success states in population.
        n: Number of draws.

    Returns:
        Natural log of the PMF value.
    """
    if k < 0 or k > K or k > n or (n - k) > (N - K):
        return float("-inf")

    return _log_comb(K, k) + _log_comb(N - K, n - k) - _log_comb(N, n)


def _log_comb(n: int, k: int) -> float:
    """Compute log of binomial coefficient C(n, k) using lgamma.

    Args:
        n: Total items.
        k: Items to choose.

    Returns:
        Natural log of C(n, k).
    """
    if k < 0 or k > n:
        return float("-inf")
    if k == 0 or k == n:
        return 0.0
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
