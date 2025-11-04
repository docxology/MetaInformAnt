"""Semantic information theory methods.

This module implements semantic information measures including information
content, semantic similarity, and semantic entropy for analyzing biological
data with ontological or hierarchical structure.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any

import numpy as np


def information_content(
    term_frequencies: dict[str, int],
    term: str,
    total_terms: int | None = None
) -> float:
    """Calculate information content of a term.
    
    Information content measures the specificity of a term. More specific
    (rarer) terms have higher information content.
    
    IC(term) = -log(P(term)) where P(term) is the frequency of the term.
    
    Args:
        term_frequencies: Dictionary mapping terms to their frequencies
        term: Term to calculate information content for
        total_terms: Total number of terms (if None, sum of frequencies)
        
    Returns:
        Information content in bits. Returns 0.0 if term not found.
        
    Examples:
        >>> freqs = {"common_term": 100, "rare_term": 1}
        >>> information_content(freqs, "rare_term") > information_content(freqs, "common_term")
        True
    """
    if term not in term_frequencies or term_frequencies[term] == 0:
        return 0.0
    
    if total_terms is None:
        total_terms = sum(term_frequencies.values())
    
    if total_terms == 0:
        return 0.0
    
    p_term = term_frequencies[term] / total_terms
    
    if p_term == 0:
        return 0.0
    
    return -math.log2(p_term)


def semantic_entropy(
    term_annotations: dict[str, set[str]],
    base: float = 2.0
) -> float:
    """Calculate semantic entropy of a set of annotated entities.
    
    Measures the information content considering the semantic structure
    of annotations (e.g., GO terms, pathways).
    
    Args:
        term_annotations: Dictionary mapping entity IDs to sets of terms
        base: Logarithm base
        
    Returns:
        Semantic entropy in bits (or nats/dits depending on base)
    """
    if not term_annotations:
        return 0.0
    
    # Count term frequencies
    term_counts: dict[str, int] = defaultdict(int)
    for terms in term_annotations.values():
        for term in terms:
            term_counts[term] += 1
    
    total_entities = len(term_annotations)
    
    # Calculate entropy using term probabilities
    entropy = 0.0
    for term, count in term_counts.items():
        p_term = count / total_entities
        if p_term > 0:
            entropy -= p_term * math.log(p_term, base)
    
    return entropy


def semantic_similarity(
    term1: str,
    term2: str,
    term_ic: dict[str, float],
    term_hierarchy: dict[str, set[str]] | None = None
) -> float:
    """Calculate semantic similarity between two terms.
    
    Uses information content to measure similarity. Common approaches:
    - Resnik similarity: IC of most informative common ancestor
    - Lin similarity: 2 * IC(MICA) / (IC(term1) + IC(term2))
    
    Args:
        term1: First term
        term2: Second term
        term_ic: Dictionary mapping terms to information content
        term_hierarchy: Optional dictionary mapping terms to parent terms
            (for finding common ancestors)
        
    Returns:
        Semantic similarity score (0 to 1, or based on IC scale)
        
    Examples:
        >>> ic = {"A": 2.0, "B": 2.0, "C": 1.0}  # C is parent of A and B
        >>> hierarchy = {"A": {"C"}, "B": {"C"}}
        >>> sim = semantic_similarity("A", "B", ic, hierarchy)
        >>> sim > 0  # Should have some similarity
        True
    """
    if term1 == term2:
        return 1.0
    
    ic1 = term_ic.get(term1, 0.0)
    ic2 = term_ic.get(term2, 0.0)
    
    if ic1 == 0.0 or ic2 == 0.0:
        return 0.0
    
    # If hierarchy provided, find most informative common ancestor (MICA)
    if term_hierarchy:
        # Find common ancestors
        ancestors1 = set()
        ancestors2 = set()
        
        def get_ancestors(term: str, anc_set: set[str]) -> None:
            if term in term_hierarchy:
                for parent in term_hierarchy[term]:
                    anc_set.add(parent)
                    get_ancestors(parent, anc_set)
        
        get_ancestors(term1, ancestors1)
        get_ancestors(term2, ancestors2)
        
        common_ancestors = ancestors1 & ancestors2
        
        if common_ancestors:
            # Find MICA (term with highest IC)
            mica_ic = max(term_ic.get(anc, 0.0) for anc in common_ancestors)
            # Lin similarity: 2 * IC(MICA) / (IC(term1) + IC(term2))
            return (2.0 * mica_ic) / (ic1 + ic2) if (ic1 + ic2) > 0 else 0.0
    
    # Simple similarity based on shared information
    # If no hierarchy, use a simple distance measure
    min_ic = min(ic1, ic2)
    max_ic = max(ic1, ic2)
    
    if max_ic == 0:
        return 0.0
    
    # Normalized similarity
    return min_ic / max_ic


def semantic_similarity_matrix(
    terms: list[str],
    term_ic: dict[str, float],
    term_hierarchy: dict[str, set[str]] | None = None
) -> np.ndarray:
    """Calculate pairwise semantic similarity matrix.
    
    Args:
        terms: List of terms to compare
        term_ic: Dictionary mapping terms to information content
        term_hierarchy: Optional dictionary mapping terms to parent terms
        
    Returns:
        Symmetric similarity matrix of shape (len(terms), len(terms))
    """
    n = len(terms)
    matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i, n):
            sim = semantic_similarity(terms[i], terms[j], term_ic, term_hierarchy)
            matrix[i, j] = sim
            matrix[j, i] = sim  # Symmetric
    
    return matrix


def information_content_from_annotations(
    term_annotations: dict[str, set[str]]
) -> dict[str, float]:
    """Calculate information content for all terms from annotations.
    
    Args:
        term_annotations: Dictionary mapping entity IDs to sets of terms
        
    Returns:
        Dictionary mapping terms to their information content
    """
    # Count term frequencies
    term_counts: dict[str, int] = defaultdict(int)
    for terms in term_annotations.values():
        for term in terms:
            term_counts[term] += 1
    
    total_entities = len(term_annotations)
    
    # Calculate IC for each term
    term_ic: dict[str, float] = {}
    for term, count in term_counts.items():
        if total_entities > 0:
            p_term = count / total_entities
            if p_term > 0:
                term_ic[term] = -math.log2(p_term)
            else:
                term_ic[term] = 0.0
        else:
            term_ic[term] = 0.0
    
    return term_ic

