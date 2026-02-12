"""Semantic information theory measures for biological ontologies.

This module implements semantic information measures including information content,
semantic similarity, and entropy measures for hierarchical biological knowledge
structures like Gene Ontology.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Union

from metainformant.core.data import validation
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def information_content(term_frequencies: Dict[str, int], term: str, total_terms: Optional[int] = None) -> float:
    """Calculate information content for a term based on frequency.

    Information content measures how specific or informative a term is.
    Higher values indicate more specific, less frequent terms.

    Args:
        term_frequencies: Dictionary mapping terms to their frequencies
        term: Term to calculate IC for
        total_terms: Total number of terms (calculated if not provided)

    Returns:
        Information content value (bits)

    Raises:
        ValueError: If term not found in frequencies
    """
    validation.validate_type(term_frequencies, dict, "term_frequencies")
    validation.validate_type(term, str, "term")

    if term not in term_frequencies:
        raise ValueError(f"Term '{term}' not found in frequency dictionary")

    freq = term_frequencies[term]

    if total_terms is None:
        total_terms = sum(term_frequencies.values())

    if total_terms == 0 or freq == 0:
        return 0.0

    # IC = -log2(p(term))
    probability = freq / total_terms
    return -math.log2(probability)


def information_content_from_annotations(annotations: Dict[str, Set[str]], term: str) -> float:
    """Calculate information content from gene-term annotation data.

    Args:
        annotations: Dictionary mapping genes/proteins to sets of terms they are annotated with
        term: Term to calculate IC for

    Returns:
        Information content value

    Raises:
        ValueError: If term not found in annotations
    """
    validation.validate_type(annotations, dict, "annotations")
    validation.validate_type(term, str, "term")

    # Count how many genes are annotated with this term
    term_count = sum(1 for gene_terms in annotations.values() if term in gene_terms)
    total_genes = len(annotations)

    if total_genes == 0:
        return 0.0

    if term_count == 0:
        raise ValueError(f"Term '{term}' not found in annotations")

    # IC = -log2(p(term))
    probability = term_count / total_genes
    return -math.log2(probability)


def semantic_entropy(term_annotations: Dict[str, Set[str]], base: float = 2.0) -> float:
    """Calculate semantic entropy of a set of term annotations.

    This measures the uncertainty or diversity in the semantic content
    of a set of annotations.

    Args:
        term_annotations: Dictionary mapping items to sets of terms
        base: Logarithm base

    Returns:
        Semantic entropy value
    """
    validation.validate_type(term_annotations, dict, "term_annotations")

    if not term_annotations:
        return 0.0

    # Collect all unique terms
    all_terms = set()
    for terms in term_annotations.values():
        all_terms.update(terms)

    if not all_terms:
        return 0.0

    # Calculate term frequencies
    term_freq = defaultdict(int)
    for terms in term_annotations.values():
        for term in terms:
            term_freq[term] += 1

    # Convert to probabilities
    total_annotations = sum(len(terms) for terms in term_annotations.values())
    if total_annotations == 0:
        return 0.0

    probabilities = [freq / total_annotations for freq in term_freq.values()]

    # Calculate entropy
    entropy = 0.0
    for p in probabilities:
        if p > 0:
            entropy -= p * math.log(p) / math.log(base)

    return entropy


def semantic_similarity(
    term1: str, term2: str, term_ic: Dict[str, float], hierarchy: Dict[str, Set[str]], method: str = "resnik"
) -> float:
    """Calculate semantic similarity between two terms.

    Args:
        term1: First term
        term2: Second term
        term_ic: Dictionary mapping terms to their information content
        hierarchy: Dictionary mapping terms to their ancestor terms
        method: Similarity method ('resnik', 'lin', 'jiang_conrath')

    Returns:
        Semantic similarity score (0-1 range, higher = more similar)

    Raises:
        ValueError: If terms not found or invalid method
    """
    validation.validate_type(term1, str, "term1")
    validation.validate_type(term2, str, "term2")
    validation.validate_type(term_ic, dict, "term_ic")
    validation.validate_type(hierarchy, dict, "hierarchy")

    if term1 not in term_ic or term2 not in term_ic:
        raise ValueError("Both terms must be present in term_ic dictionary")

    if method not in ["resnik", "lin", "jiang_conrath"]:
        raise ValueError(f"Unknown similarity method: {method}")

    # Find common ancestors
    ancestors1 = _get_ancestors(term1, hierarchy)
    ancestors2 = _get_ancestors(term2, hierarchy)

    common_ancestors = ancestors1.intersection(ancestors2)

    if not common_ancestors:
        return 0.0

    if method == "resnik":
        # Maximum IC of common ancestors
        max_ic = max(term_ic.get(ca, 0.0) for ca in common_ancestors)
        # Normalize by maximum possible IC
        max_possible_ic = max(term_ic.values()) if term_ic else 1.0
        return max_ic / max_possible_ic if max_possible_ic > 0 else 0.0

    elif method == "lin":
        # Lin similarity: 2 * IC(lcs) / (IC(t1) + IC(t2))
        # Use the most informative common ancestor (highest IC)
        lcs_ic = max(term_ic.get(ca, 0.0) for ca in common_ancestors)
        ic1 = term_ic[term1]
        ic2 = term_ic[term2]

        denominator = ic1 + ic2
        if denominator == 0:
            return 0.0

        return 2.0 * lcs_ic / denominator

    elif method == "jiang_conrath":
        # Jiang-Conrath distance: IC(t1) + IC(t2) - 2 * IC(lcs)
        lcs_ic = max(term_ic.get(ca, 0.0) for ca in common_ancestors)
        ic1 = term_ic[term1]
        ic2 = term_ic[term2]

        distance = ic1 + ic2 - 2.0 * lcs_ic

        # Convert distance to similarity (simple exponential decay)
        # Higher distance = lower similarity
        return math.exp(-distance) if distance >= 0 else 0.0


def semantic_similarity_matrix(
    terms: List[str], term_ic: Dict[str, float], hierarchy: Dict[str, Set[str]], method: str = "resnik"
) -> List[List[float]]:
    """Calculate pairwise semantic similarity matrix for a set of terms.

    Args:
        terms: List of terms to compare
        term_ic: Dictionary mapping terms to their information content
        hierarchy: Dictionary mapping terms to their ancestor terms
        method: Similarity method

    Returns:
        Square similarity matrix (n_terms x n_terms)

    Raises:
        ValueError: If terms list is empty
    """
    validation.validate_type(terms, list, "terms")

    n_terms = len(terms)
    if n_terms == 0:
        raise ValueError("Terms list cannot be empty")

    similarity_matrix = [[0.0] * n_terms for _ in range(n_terms)]

    for i in range(n_terms):
        for j in range(i, n_terms):  # Only compute upper triangle
            if i == j:
                similarity_matrix[i][j] = 1.0  # Self-similarity
            else:
                try:
                    sim = semantic_similarity(terms[i], terms[j], term_ic, hierarchy, method)
                    similarity_matrix[i][j] = sim
                    similarity_matrix[j][i] = sim  # Symmetric
                except ValueError:
                    # If similarity calculation fails, set to 0
                    similarity_matrix[i][j] = 0.0
                    similarity_matrix[j][i] = 0.0

    return similarity_matrix


def _get_ancestors(term: str, hierarchy: Dict[str, Set[str]]) -> Set[str]:
    """Get all ancestors of a term in the hierarchy.

    Args:
        term: Term to find ancestors for
        hierarchy: Dictionary mapping terms to their parent terms

    Returns:
        Set of all ancestor terms including the term itself
    """
    ancestors = {term}
    queue = [term]

    while queue:
        current = queue.pop(0)
        if current in hierarchy:
            parents = hierarchy[current]
            for parent in parents:
                if parent not in ancestors:
                    ancestors.add(parent)
                    queue.append(parent)

    return ancestors


def semantic_distance(
    term1: str, term2: str, term_ic: Dict[str, float], hierarchy: Dict[str, Set[str]], method: str = "jiang_conrath"
) -> float:
    """Calculate semantic distance between two terms.

    Args:
        term1: First term
        term2: Second term
        term_ic: Dictionary mapping terms to information content
        hierarchy: Dictionary mapping terms to ancestors
        method: Distance method ('jiang_conrath', 'lin')

    Returns:
        Semantic distance (higher = more dissimilar)

    Raises:
        ValueError: If invalid method or terms not found
    """
    if method == "jiang_conrath":
        # Jiang-Conrath distance: IC(t1) + IC(t2) - 2 * IC(lcs)
        ancestors1 = _get_ancestors(term1, hierarchy)
        ancestors2 = _get_ancestors(term2, hierarchy)
        common_ancestors = ancestors1.intersection(ancestors2)

        if not common_ancestors:
            return float("inf")

        lcs_ic = max(term_ic.get(ca, 0.0) for ca in common_ancestors)
        ic1 = term_ic[term1]
        ic2 = term_ic[term2]

        return ic1 + ic2 - 2.0 * lcs_ic

    elif method == "lin":
        # Lin distance: 1 - Lin similarity
        similarity = semantic_similarity(term1, term2, term_ic, hierarchy, method="lin")
        return 1.0 - similarity

    else:
        raise ValueError(f"Unknown distance method: {method}")


def term_specificity(term: str, term_ic: Dict[str, float], hierarchy: Dict[str, Set[str]]) -> float:
    """Calculate term specificity based on information content and hierarchy.

    Args:
        term: Term to analyze
        term_ic: Dictionary mapping terms to information content
        hierarchy: Dictionary mapping terms to ancestors

    Returns:
        Specificity score (0-1, higher = more specific)

    Raises:
        ValueError: If term not found
    """
    if term not in term_ic:
        raise ValueError(f"Term '{term}' not found in term_ic dictionary")

    # Get all ancestors
    ancestors = _get_ancestors(term, hierarchy)

    # Calculate average IC of ancestors
    ancestor_ics = [term_ic.get(ancestor, 0.0) for ancestor in ancestors]
    avg_ancestor_ic = sum(ancestor_ics) / len(ancestor_ics) if ancestor_ics else 0.0

    # Specificity = IC(term) / (IC(term) + avg_ancestor_ic)
    term_ic_val = term_ic[term]
    denominator = term_ic_val + avg_ancestor_ic

    if denominator == 0:
        return 0.0

    return term_ic_val / denominator


def ontology_complexity(hierarchy: Dict[str, Set[str]], term_ic: Optional[Dict[str, float]] = None) -> Dict[str, float]:
    """Calculate complexity metrics for an ontology.

    Args:
        hierarchy: Dictionary mapping terms to their parent terms
        term_ic: Optional dictionary of information content values

    Returns:
        Dictionary with complexity metrics
    """
    validation.validate_type(hierarchy, dict, "hierarchy")

    # Basic structure metrics
    n_terms = len(hierarchy)
    n_relationships = sum(len(parents) for parents in hierarchy.values())

    # Calculate depth for each term
    depths = {}
    for term in hierarchy:
        depths[term] = _calculate_term_depth(term, hierarchy)

    avg_depth = sum(depths.values()) / len(depths) if depths else 0.0
    max_depth = max(depths.values()) if depths else 0

    # Calculate information content distribution if provided
    ic_stats = {}
    if term_ic:
        ic_values = list(term_ic.values())
        if ic_values:
            ic_mean = sum(ic_values) / len(ic_values)
            ic_stats = {
                "ic_mean": ic_mean,
                "ic_std": math.sqrt(sum((x - ic_mean) ** 2 for x in ic_values) / len(ic_values)),
                "ic_min": min(ic_values),
                "ic_max": max(ic_values),
                "ic_range": max(ic_values) - min(ic_values),
            }

    return {
        "n_terms": n_terms,
        "n_relationships": n_relationships,
        "avg_depth": avg_depth,
        "max_depth": max_depth,
        "depths": depths,
        **ic_stats,
    }


def _calculate_term_depth(term: str, hierarchy: Dict[str, Set[str]]) -> int:
    """Calculate the depth of a term in the hierarchy.

    Args:
        term: Term to calculate depth for
        hierarchy: Dictionary mapping terms to parents

    Returns:
        Depth (distance from root)
    """
    if term not in hierarchy or not hierarchy[term]:
        return 0  # Root term

    # Find maximum depth among parents
    parent_depths = []
    for parent in hierarchy[term]:
        if parent in hierarchy:  # Avoid infinite recursion
            parent_depths.append(_calculate_term_depth(parent, hierarchy))

    return 1 + max(parent_depths) if parent_depths else 1


def term_redundancy(annotations: Dict[str, Set[str]], term: str) -> float:
    """Calculate term redundancy in annotations.

    This measures how redundant a term is - i.e., how often it's
    accompanied by other terms vs appearing alone.

    Args:
        annotations: Dictionary mapping items to sets of terms
        term: Term to analyze

    Returns:
        Redundancy score (0-1, higher = more redundant)

    Raises:
        ValueError: If term not found in annotations
    """
    validation.validate_type(annotations, dict, "annotations")
    validation.validate_type(term, str, "term")

    term_occurrences = 0
    term_alone_occurrences = 0

    for item_terms in annotations.values():
        if term in item_terms:
            term_occurrences += 1
            if len(item_terms) == 1:
                term_alone_occurrences += 1

    if term_occurrences == 0:
        raise ValueError(f"Term '{term}' not found in annotations")

    # Redundancy = 1 - (times_term_appears_alone / total_times_term_appears)
    return 1.0 - (term_alone_occurrences / term_occurrences)


def annotation_specificity(annotations: Dict[str, Set[str]], term_ic: Dict[str, float]) -> Dict[str, float]:
    """Calculate specificity scores for all annotations.

    Args:
        annotations: Dictionary mapping items to term sets
        term_ic: Dictionary mapping terms to information content

    Returns:
        Dictionary mapping items to their specificity scores
    """
    validation.validate_type(annotations, dict, "annotations")
    validation.validate_type(term_ic, dict, "term_ic")

    specificity_scores = {}

    for item, terms in annotations.items():
        if not terms:
            specificity_scores[item] = 0.0
            continue

        # Average IC of all terms for this item
        ic_values = [term_ic.get(term, 0.0) for term in terms]
        avg_ic = sum(ic_values) / len(ic_values)

        specificity_scores[item] = avg_ic

    return specificity_scores
