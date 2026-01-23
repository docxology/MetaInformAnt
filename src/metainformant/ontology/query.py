"""Ontology querying and traversal.

This module provides functions to query, traverse, and analyze ontologies,
including finding ancestors/descendants, common ancestors, paths, and various
statistics and filters.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Set, Optional, Any, Iterable, Tuple
from metainformant.core import logging, errors, validation
from .types import Ontology, Term, Relationship

logger = logging.get_logger(__name__)

# Global cache for query results
_query_cache: Dict[str, Any] = {}
_cache_enabled = True


def clear_cache() -> None:
    """Clear the ontology query cache."""
    global _query_cache
    _query_cache = {}
    logger.info("Ontology query cache cleared.")


def set_cache_enabled(enabled: bool) -> None:
    """Enable or disable the ontology query cache."""
    global _cache_enabled
    _cache_enabled = enabled
    logger.info(f"Ontology query cache {'enabled' if enabled else 'disabled'}.")


def _get_cache_key(func_name: str, *args, **kwargs) -> str:
    """Generate a cache key for a function call."""
    key_parts = [func_name]
    key_parts.extend(str(arg) for arg in args)
    key_parts.extend(f"{k}={v}" for k, v in sorted(kwargs.items()))
    return "|".join(key_parts)


def ancestors(onto: Ontology, term_id: str, relation_type: str = "is_a") -> Set[str]:
    """Find all ancestors of a term via a specific relationship type.

    Args:
        onto: Ontology to query.
        term_id: ID of the term to find ancestors for.
        relation_type: Type of relationship to traverse (default: "is_a").

    Returns:
        Set of ancestor term IDs, including the term itself.

    Raises:
        errors.TermNotFoundError: If term_id is not in the ontology.

    Examples:
        >>> ancestors(go_ontology, "GO:0008150")
        {'GO:0008150', 'GO:0003674', 'GO:0005575'}
    """
    validation.validate_type(onto, Ontology, "onto")
    validation.validate_not_empty(term_id, "term_id")
    if term_id not in onto.terms:
        raise errors.TermNotFoundError(f"Term '{term_id}' not found in ontology.")

    cache_key = _get_cache_key("ancestors", term_id, relation_type)
    if _cache_enabled and cache_key in _query_cache:
        return set(_query_cache[cache_key])

    ancestors_set = {term_id}  # Include self
    to_process = [term_id]
    visited = set()

    while to_process:
        current_term = to_process.pop()
        if current_term in visited:
            continue
        visited.add(current_term)

        # Find parents via specified relationship type
        for rel in onto.relationships:
            if rel.target == current_term and rel.relation_type == relation_type:
                if rel.source not in ancestors_set:
                    ancestors_set.add(rel.source)
                    to_process.append(rel.source)

    if _cache_enabled:
        _query_cache[cache_key] = list(ancestors_set)

    return ancestors_set


def descendants(onto: Ontology, term_id: str, relation_type: str = "is_a") -> Set[str]:
    """Find all descendants of a term via a specific relationship type.

    Args:
        onto: Ontology to query.
        term_id: ID of the term to find descendants for.
        relation_type: Type of relationship to traverse (default: "is_a").

    Returns:
        Set of descendant term IDs, including the term itself.

    Raises:
        errors.TermNotFoundError: If term_id is not in the ontology.

    Examples:
        >>> descendants(go_ontology, "GO:0003674")
        {'GO:0003674', 'GO:0008150', 'GO:0016740', ...}
    """
    validation.validate_type(onto, Ontology, "onto")
    validation.validate_not_empty(term_id, "term_id")
    if term_id not in onto.terms:
        raise errors.TermNotFoundError(f"Term '{term_id}' not found in ontology.")

    cache_key = _get_cache_key("descendants", term_id, relation_type)
    if _cache_enabled and cache_key in _query_cache:
        return set(_query_cache[cache_key])

    descendants_set = {term_id}  # Include self
    to_process = [term_id]
    visited = set()

    while to_process:
        current_term = to_process.pop()
        if current_term in visited:
            continue
        visited.add(current_term)

        # Find children via specified relationship type
        for rel in onto.relationships:
            if rel.source == current_term and rel.relation_type == relation_type:
                if rel.target not in descendants_set:
                    descendants_set.add(rel.target)
                    to_process.append(rel.target)

    if _cache_enabled:
        _query_cache[cache_key] = list(descendants_set)

    return descendants_set


def common_ancestors(onto: Ontology, term1: str, term2: str, relation_type: str = "is_a") -> Set[str]:
    """Find common ancestors of two terms.

    Args:
        onto: Ontology to query.
        term1: ID of first term.
        term2: ID of second term.
        relation_type: Type of relationship to traverse.

    Returns:
        Set of common ancestor term IDs.

    Examples:
        >>> common_ancestors(go_ontology, "GO:0003674", "GO:0008150")
        {'GO:0003674', 'GO:0005575'}
    """
    ancestors1 = ancestors(onto, term1, relation_type)
    ancestors2 = ancestors(onto, term2, relation_type)
    return ancestors1.intersection(ancestors2)


def most_informative_common_ancestor(
    onto: Ontology, term1: str, term2: str, ic_map: Dict[str, float], relation_type: str = "is_a"
) -> str:
    """Find the most informative common ancestor (MICA) of two terms.

    Args:
        onto: Ontology to query.
        term1: ID of first term.
        term2: ID of second term.
        ic_map: Dictionary mapping term IDs to information content values.
        relation_type: Type of relationship to traverse.

    Returns:
        ID of the most informative common ancestor.

    Raises:
        ValueError: If no common ancestors exist.

    Examples:
        >>> mica = most_informative_common_ancestor(go_ontology, "GO:0003674", "GO:0008150", ic_map)
    """
    common_ancs = common_ancestors(onto, term1, term2, relation_type)
    if not common_ancs:
        raise ValueError(f"No common ancestors found for terms {term1} and {term2}")

    # Find the common ancestor with highest information content
    mica = max(common_ancs, key=lambda t: ic_map.get(t, 0.0))
    return mica


def path_to_root(onto: Ontology, term_id: str, relation_type: str = "is_a") -> List[str]:
    """Find one path from a term to a root term.

    Args:
        onto: Ontology to query.
        term_id: ID of the term to find path for.
        relation_type: Type of relationship to traverse.

    Returns:
        List of term IDs representing a path from term to root, or empty list if no path.

    Examples:
        >>> path_to_root(go_ontology, "GO:0008150")
        ['GO:0008150', 'GO:0003674', 'GO:0005575']
    """
    validation.validate_type(onto, Ontology, "onto")
    validation.validate_not_empty(term_id, "term_id")

    # Find root terms (terms with no parents)
    all_targets = {rel.target for rel in onto.relationships if rel.relation_type == relation_type}
    roots = [t for t in onto.terms.keys() if t not in all_targets]

    if term_id in roots:
        return [term_id]

    # Use BFS to find shortest path to any root
    visited = set()
    queue = [(term_id, [term_id])]  # (current_term, path_to_current)

    while queue:
        current_term, path = queue.pop(0)

        if current_term in visited:
            continue
        visited.add(current_term)

        # Check if current term has parents
        has_parents = any(
            rel.target == current_term and rel.relation_type == relation_type for rel in onto.relationships
        )

        if not has_parents:
            # Found a root
            return path

        # Add parents to queue
        for rel in onto.relationships:
            if rel.target == current_term and rel.relation_type == relation_type:
                if rel.source not in visited:
                    new_path = path + [rel.source]
                    queue.append((rel.source, new_path))

    return []  # No path found


def shortest_path(onto: Ontology, term1: str, term2: str, relation_type: str = "is_a") -> List[str]:
    """Find the shortest path between two terms.

    Args:
        onto: Ontology to query.
        term1: ID of first term.
        term2: ID of second term.
        relation_type: Type of relationship to traverse.

    Returns:
        List of term IDs representing the shortest path, or empty list if no path.

    Examples:
        >>> shortest_path(go_ontology, "GO:0003674", "GO:0008150")
        ['GO:0003674', 'GO:0008150']
    """
    # Build adjacency list
    adj_list = {}
    for rel in onto.relationships:
        if rel.relation_type == relation_type:
            if rel.source not in adj_list:
                adj_list[rel.source] = []
            if rel.target not in adj_list:
                adj_list[rel.target] = []
            adj_list[rel.source].append(rel.target)
            # For undirected search, add reverse edge
            adj_list[rel.target].append(rel.source)

    # BFS to find shortest path
    visited = set()
    queue = [(term1, [term1])]  # (current_term, path_to_current)

    while queue:
        current_term, path = queue.pop(0)

        if current_term == term2:
            return path

        if current_term in visited:
            continue
        visited.add(current_term)

        for neighbor in adj_list.get(current_term, []):
            if neighbor not in visited:
                new_path = path + [neighbor]
                queue.append((neighbor, new_path))

    return []  # No path found


def get_subontology(onto: Ontology, root_terms: Iterable[str], relation_type: str = "is_a") -> Ontology:
    """Extract a sub-ontology rooted at specified terms.

    Args:
        onto: Ontology to extract from.
        root_terms: IDs of root terms for the sub-ontology.
        relation_type: Type of relationship to traverse.

    Returns:
        New Ontology object containing the sub-ontology.

    Examples:
        >>> bp_ontology = get_subontology(go_ontology, ["GO:0008150"])
    """
    validation.validate_type(onto, Ontology, "onto")
    validation.validate_type(root_terms, Iterable, "root_terms")

    all_terms = set()
    for root in root_terms:
        all_terms.update(descendants(onto, root, relation_type))

    # Extract terms and relationships
    sub_terms = {tid: onto.terms[tid] for tid in all_terms if tid in onto.terms}
    sub_relationships = [
        rel
        for rel in onto.relationships
        if rel.source in all_terms and rel.target in all_terms and rel.relation_type == relation_type
    ]

    from .types import create_ontology

    sub_onto = create_ontology(terms=sub_terms, relationships=sub_relationships, **onto.metadata)

    logger.info(f"Extracted sub-ontology with {len(sub_onto)} terms from {len(onto)} total terms")
    return sub_onto


def find_terms_by_name(onto: Ontology, name_pattern: str, case_sensitive: bool = False) -> List[str]:
    """Find terms by name pattern.

    Args:
        onto: Ontology to search.
        name_pattern: Pattern to match in term names.
        case_sensitive: Whether matching should be case-sensitive.

    Returns:
        List of matching term IDs.

    Examples:
        >>> find_terms_by_name(go_ontology, "biological_process")
        ['GO:0008150']
    """
    validation.validate_type(onto, Ontology, "onto")
    validation.validate_not_empty(name_pattern, "name_pattern")

    pattern = name_pattern if case_sensitive else name_pattern.lower()
    matches = []

    for term in onto.terms.values():
        if term.name:
            term_name = term.name if case_sensitive else term.name.lower()
            if pattern in term_name:
                matches.append(term.id)

    return matches


def find_terms_by_namespace(onto: Ontology, namespace: str) -> List[str]:
    """Find all terms in a specific namespace.

    Args:
        onto: Ontology to search.
        namespace: Namespace to filter by.

    Returns:
        List of term IDs in the specified namespace.

    Examples:
        >>> bp_terms = find_terms_by_namespace(go_ontology, "biological_process")
    """
    validation.validate_type(onto, Ontology, "onto")
    validation.validate_not_empty(namespace, "namespace")

    return [term.id for term in onto.terms.values() if term.namespace == namespace]


def get_roots(onto: Ontology, relation_type: str = "is_a") -> Set[str]:
    """Find all root terms (terms with no parents) in the ontology.

    Args:
        onto: Ontology to analyze.
        relation_type: Type of relationship to consider.

    Returns:
        Set of root term IDs.

    Examples:
        >>> roots = get_roots(go_ontology)
        >>> print(roots)  # {'GO:0003674', 'GO:0005575', 'GO:0008150'}
    """
    validation.validate_type(onto, Ontology, "onto")

    # Find all terms that are targets of relationships
    targets = {rel.target for rel in onto.relationships if rel.relation_type == relation_type}
    roots = {term.id for term in onto.terms.values() if term.id not in targets}
    return roots


def get_leaves(onto: Ontology, relation_type: str = "is_a") -> Set[str]:
    """Find all leaf terms (terms with no children) in the ontology.

    Args:
        onto: Ontology to analyze.
        relation_type: Type of relationship to consider.

    Returns:
        Set of leaf term IDs.

    Examples:
        >>> leaves = get_leaves(go_ontology)
    """
    validation.validate_type(onto, Ontology, "onto")

    # Find all terms that are sources of relationships
    sources = {rel.source for rel in onto.relationships if rel.relation_type == relation_type}
    leaves = {term.id for term in onto.terms.values() if term.id not in sources}
    return leaves


def get_subontology_stats(onto: Ontology) -> Dict[str, Any]:
    """Calculate comprehensive statistics for an ontology.

    Args:
        onto: Ontology to analyze.

    Returns:
        Dictionary containing various ontology statistics.

    Examples:
        >>> stats = get_subontology_stats(go_ontology)
        >>> print(f"Terms: {stats['num_terms']}, Relationships: {stats['num_relationships']}")
    """
    validation.validate_type(onto, Ontology, "onto")

    stats = {
        "num_terms": len(onto.terms),
        "num_relationships": len(onto.relationships),
        "num_obsolete_terms": sum(1 for t in onto.terms.values() if t.is_obsolete),
    }

    # Namespace distribution
    namespace_counts = {}
    for term in onto.terms.values():
        ns = term.namespace or "unknown"
        namespace_counts[ns] = namespace_counts.get(ns, 0) + 1
    stats["namespace_distribution"] = namespace_counts

    # Relationship type distribution
    rel_type_counts = {}
    for rel in onto.relationships:
        rel_type_counts[rel.relation_type] = rel_type_counts.get(rel.relation_type, 0) + 1
    stats["relationship_type_distribution"] = rel_type_counts

    # Term depth statistics (using path to root)
    depths = []
    for term_id in onto.terms.keys():
        path = path_to_root(onto, term_id)
        if path:
            depths.append(len(path) - 1)  # Depth is path length minus 1

    if depths:
        stats["depth_stats"] = {
            "mean": sum(depths) / len(depths),
            "max": max(depths),
            "min": min(depths),
        }

    # Root and leaf counts
    roots = get_roots(onto)
    leaves = get_leaves(onto)
    stats["num_roots"] = len(roots)
    stats["num_leaves"] = len(leaves)

    # Terms with definitions
    terms_with_defs = sum(1 for t in onto.terms.values() if t.definition)
    stats["terms_with_definitions"] = terms_with_defs

    return stats


def validate_ontology_integrity(onto: Ontology) -> Tuple[bool, List[str]]:
    """Validate the structural integrity of an ontology.

    Args:
        onto: Ontology to validate.

    Returns:
        Tuple of (is_valid, list_of_error_messages).

    Examples:
        >>> is_valid, errors = validate_ontology_integrity(go_ontology)
        >>> if not is_valid:
        ...     print("Ontology has issues:", errors)
    """
    validation.validate_type(onto, Ontology, "onto")

    errors_list = []

    # Check that all relationship terms exist
    for rel in onto.relationships:
        if rel.source not in onto.terms:
            errors_list.append(f"Relationship source '{rel.source}' not found in terms")
        if rel.target not in onto.terms:
            errors_list.append(f"Relationship target '{rel.target}' not found in terms")

    # Check for cycles (simplified check - full cycle detection would be more complex)
    # This is a basic check for direct self-loops
    for rel in onto.relationships:
        if rel.source == rel.target:
            errors_list.append(f"Self-loop relationship found: {rel.source} -> {rel.target}")

    # Check for duplicate relationships
    rel_set = set()
    for rel in onto.relationships:
        rel_key = (rel.source, rel.target, rel.relation_type)
        if rel_key in rel_set:
            errors_list.append(f"Duplicate relationship: {rel_key}")
        rel_set.add(rel_key)

    is_valid = len(errors_list) == 0
    return is_valid, errors_list


def information_content(onto: Ontology, term_id: str, corpus_size: int | None = None) -> float:
    """Calculate the information content of a term.

    Information content is -log(probability of term occurrence).

    Args:
        onto: Ontology to analyze.
        term_id: ID of the term to calculate IC for.
        corpus_size: Total number of annotations in the corpus. If None, uses number of terms.

    Returns:
        Information content value (higher = more specific term).

    Examples:
        >>> ic = information_content(go_ontology, "GO:0008150", corpus_size=100000)
    """
    validation.validate_type(onto, Ontology, "onto")
    validation.validate_not_empty(term_id, "term_id")

    if term_id not in onto.terms:
        raise errors.TermNotFoundError(f"Term '{term_id}' not found in ontology.")

    # For now, use a simple frequency-based approach
    # In a real implementation, this would use annotation frequencies from a corpus
    descendants_set = descendants(onto, term_id)
    term_frequency = len(descendants_set)

    if corpus_size is None:
        corpus_size = len(onto.terms)

    if term_frequency == 0:
        return 0.0

    import math

    probability = term_frequency / corpus_size
    ic = -math.log(probability) if probability > 0 else 0.0

    return ic


def calculate_ic_map(onto: Ontology, corpus_size: int | None = None) -> Dict[str, float]:
    """Calculate information content for all terms in the ontology.

    Args:
        onto: Ontology to analyze.
        corpus_size: Total corpus size for IC calculation.

    Returns:
        Dictionary mapping term IDs to information content values.

    Examples:
        >>> ic_map = calculate_ic_map(go_ontology)
    """
    ic_map = {}
    for term_id in onto.terms.keys():
        ic_map[term_id] = information_content(onto, term_id, corpus_size)
    return ic_map


def subgraph(onto: Ontology, term_ids: List[str], relation_type: str = "is_a") -> Ontology:
    """Create a subgraph ontology containing only specified terms and their relationships.

    Args:
        onto: Source ontology.
        term_ids: List of term IDs to include.
        relation_type: Type of relationship to consider.

    Returns:
        New Ontology containing only the specified terms and relationships between them.

    Examples:
        >>> subgraph_onto = subgraph(go_ontology, ["GO:0008150", "GO:0003674"])
    """
    validation.validate_type(onto, Ontology, "onto")
    validation.validate_not_empty(term_ids, "term_ids")

    # Create new ontology with selected terms
    selected_terms = {}
    for term_id in term_ids:
        if term_id in onto.terms:
            selected_terms[term_id] = onto.terms[term_id]

    # Find relationships between selected terms
    selected_relationships = []
    for rel in onto.relationships:
        if rel.relation_type == relation_type and rel.source in selected_terms and rel.target in selected_terms:
            selected_relationships.append(rel)

    return Ontology(
        terms=selected_terms,
        relationships=selected_relationships,
        metadata={**onto.metadata, "subgraph_of": onto.metadata.get("id", "unknown")},
    )


def path_to_root(onto: Ontology, term_id: str, relation_type: str = "is_a") -> List[str]:
    """Find the path from a term to the root of the ontology.

    Args:
        onto: Ontology to analyze.
        term_id: Starting term ID.
        relation_type: Type of relationship to traverse.

    Returns:
        List of term IDs from the given term to the root (inclusive).

    Examples:
        >>> path = path_to_root(go_ontology, "GO:0008150")
        >>> print(path)  # ['GO:0008150', 'GO:0003674', 'GO:0005575']
    """
    validation.validate_type(onto, Ontology, "onto")
    validation.validate_not_empty(term_id, "term_id")

    if term_id not in onto.terms:
        raise errors.TermNotFoundError(f"Term '{term_id}' not found in ontology.")

    path = [term_id]
    current = term_id

    # Traverse up to parents until we reach a root
    visited = set()
    while current not in visited:
        visited.add(current)
        parents = [
            rel.source for rel in onto.relationships if rel.target == current and rel.relation_type == relation_type
        ]

        if not parents:
            break  # Reached a root

        # For simplicity, follow the first parent (in real ontologies, there might be multiple paths)
        current = parents[0]
        if current not in path:  # Avoid cycles
            path.append(current)

    return path


def distance(onto: Ontology, term1: str, term2: str, relation_type: str = "is_a") -> int:
    """Calculate the graph distance between two terms.

    Args:
        onto: Ontology to analyze.
        term1: First term ID.
        term2: Second term ID.
        relation_type: Type of relationship to consider.

    Returns:
        Minimum number of edges between the terms, or -1 if not connected.

    Examples:
        >>> dist = distance(go_ontology, "GO:0008150", "GO:0003674")
    """
    validation.validate_type(onto, Ontology, "onto")
    validation.validate_not_empty(term1, "term1")
    validation.validate_not_empty(term2, "term2")

    if term1 not in onto.terms or term2 not in onto.terms:
        return -1

    if term1 == term2:
        return 0

    # Simple BFS to find shortest path
    from collections import deque

    visited = set()
    queue = deque([(term1, 0)])

    while queue:
        current, dist = queue.popleft()

        if current in visited:
            continue
        visited.add(current)

        # Check relationships
        for rel in onto.relationships:
            if rel.relation_type == relation_type:
                # Check both directions (ontology might not be directed)
                if rel.source == current and rel.target not in visited:
                    if rel.target == term2:
                        return dist + 1
                    queue.append((rel.target, dist + 1))
                elif rel.target == current and rel.source not in visited:
                    if rel.source == term2:
                        return dist + 1
                    queue.append((rel.source, dist + 1))

    return -1  # Not connected


def find_term_by_name(onto: Ontology, name_pattern: str, case_sensitive: bool = False) -> List[str]:
    """Find terms by name pattern.

    Args:
        onto: Ontology to search.
        name_pattern: Pattern to match (case-insensitive substring search).
        case_sensitive: Whether to perform case-sensitive matching.

    Returns:
        List of term IDs matching the pattern.

    Examples:
        >>> terms = find_term_by_name(go_ontology, "transport")
    """
    validation.validate_type(onto, Ontology, "onto")

    matches = []
    pattern = name_pattern if case_sensitive else name_pattern.lower()

    for term_id, term in onto.terms.items():
        name = term.name or ""
        if not case_sensitive:
            name = name.lower()

        if pattern in name:
            matches.append(term_id)

    return matches


def filter_by_namespace(onto: Ontology, namespace: str) -> Ontology:
    """Create a new ontology containing only terms from the specified namespace.

    Args:
        onto: Source ontology.
        namespace: Namespace to filter by.

    Returns:
        New Ontology with filtered terms.

    Examples:
        >>> bp_ontology = filter_by_namespace(go_ontology, "biological_process")
    """
    validation.validate_type(onto, Ontology, "onto")

    filtered_terms = {}
    for term_id, term in onto.terms.items():
        if term.namespace == namespace:
            filtered_terms[term_id] = term

    # Include relationships between filtered terms
    filtered_relationships = []
    for rel in onto.relationships:
        if rel.source in filtered_terms and rel.target in filtered_terms:
            filtered_relationships.append(rel)

    return Ontology(
        terms=filtered_terms,
        relationships=filtered_relationships,
        metadata={**onto.metadata, "filtered_by_namespace": namespace},
    )


# Cache management (simple in-memory cache)
_CACHE_ENABLED = True
_QUERY_CACHE = {}


def clear_cache() -> None:
    """Clear the query cache."""
    global _QUERY_CACHE
    _QUERY_CACHE = {}


def set_cache_enabled(enabled: bool) -> None:
    """Enable or disable query caching.

    Args:
        enabled: Whether to enable caching.
    """
    global _CACHE_ENABLED
    _CACHE_ENABLED = enabled

    if not enabled:
        clear_cache()
