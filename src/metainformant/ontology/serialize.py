"""Ontology serialization and file I/O.

This module provides functionality to save and load ontologies in various formats,
including JSON, OBO, and NetworkX graph formats.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Any, Optional
from metainformant.core import logging, io, validation
from .types import Ontology, Term, Relationship

logger = logging.get_logger(__name__)


def save_ontology(onto: Ontology, path: str | Path, format: str = "obo") -> None:
    """Save ontology to a file.

    Args:
        onto: Ontology to save.
        path: Path to save the ontology to.
        format: Format to save in ("obo", "json").

    Raises:
        ValueError: If format is not supported.
        IOError: If writing fails.

    Examples:
        >>> save_ontology(ontology, "output/my_ontology.json", format="json")
    """
    validation.validate_not_none(onto, "onto")
    path = Path(path)

    if format not in ["obo", "json"]:
        raise ValueError(f"Unsupported format: {format}. Use 'obo' or 'json'.")

    logger.info(f"Saving ontology with {len(onto)} terms to {path} in {format} format")

    if format == "obo":
        _save_obo_format(onto, path)
    elif format == "json":
        _save_json_format(onto, path)

    logger.info(f"Ontology saved successfully to {path}")


def load_ontology(path: str | Path, format: str = "obo") -> Ontology:
    """Load ontology from a file.

    Args:
        path: Path to the ontology file.
        format: Format of the file ("obo", "json").

    Returns:
        Loaded Ontology object.

    Raises:
        ValueError: If format is not supported.
        FileNotFoundError: If file doesn't exist.

    Examples:
        >>> ontology = load_ontology("data/go.json", format="json")
    """
    path = validation.validate_path_exists(Path(path))

    if format not in ["obo", "json"]:
        raise ValueError(f"Unsupported format: {format}. Use 'obo' or 'json'.")

    logger.info(f"Loading ontology from {path} in {format} format")

    if format == "obo":
        from .obo import parse_obo
        ontology = parse_obo(path)
    elif format == "json":
        ontology = _load_json_format(path)

    logger.info(f"Ontology loaded successfully with {len(ontology)} terms")
    return ontology


def ontology_to_graph(onto: Ontology) -> Any:
    """Convert ontology to NetworkX directed graph.

    Args:
        onto: Ontology to convert.

    Returns:
        NetworkX DiGraph representing the ontology.

    Raises:
        ImportError: If networkx is not available.

    Examples:
        >>> import networkx as nx
        >>> graph = ontology_to_graph(ontology)
        >>> nx.draw(graph)  # Visualize the ontology
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx is required for ontology_to_graph. Install with: pip install networkx")

    validation.validate_not_none(onto, "onto")

    logger.info(f"Converting ontology with {len(onto)} terms to NetworkX graph")

    # Create directed graph
    graph = nx.DiGraph()

    # Add nodes (terms)
    for term_id, term in onto.terms.items():
        node_attrs = {
            'name': term.name,
            'namespace': term.namespace,
            'definition': term.definition,
            'is_obsolete': term.is_obsolete,
        }
        # Filter out None values
        node_attrs = {k: v for k, v in node_attrs.items() if v is not None}
        graph.add_node(term_id, **node_attrs)

    # Add edges (relationships)
    for rel in onto.relationships:
        graph.add_edge(rel.source, rel.target, relation_type=rel.relation_type)

    logger.info(f"Created graph with {len(graph.nodes)} nodes and {len(graph.edges)} edges")
    return graph


def graph_to_ontology(graph: Any, metadata: Dict[str, Any] | None = None) -> Ontology:
    """Convert NetworkX graph back to Ontology.

    Args:
        graph: NetworkX DiGraph representing an ontology.
        metadata: Optional metadata for the ontology.

    Returns:
        Ontology object.

    Raises:
        ImportError: If networkx is not available.
        ValueError: If graph structure is invalid.

    Examples:
        >>> ontology = graph_to_ontology(my_graph)
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx is required for graph_to_ontology. Install with: pip install networkx")

    validation.validate_not_none(graph, "graph")

    if not isinstance(graph, nx.DiGraph):
        raise ValueError("Graph must be a NetworkX DiGraph")

    logger.info(f"Converting NetworkX graph with {len(graph.nodes)} nodes to Ontology")

    # Extract terms from nodes
    terms = {}
    for node_id, node_attrs in graph.nodes(data=True):
        term = Term(
            id=node_id,
            name=node_attrs.get('name'),
            definition=node_attrs.get('definition'),
            namespace=node_attrs.get('namespace'),
            is_obsolete=node_attrs.get('is_obsolete', False)
        )
        terms[node_id] = term

    # Extract relationships from edges
    relationships = []
    for source, target, edge_attrs in graph.edges(data=True):
        rel_type = edge_attrs.get('relation_type', 'is_a')
        relationship = Relationship(
            source=source,
            target=target,
            relation_type=rel_type
        )
        relationships.append(relationship)

    # Create ontology
    from .types import create_ontology
    ontology = create_ontology(
        terms=terms,
        relationships=relationships,
        **(metadata or {})
    )

    logger.info(f"Converted graph to ontology with {len(ontology)} terms")
    return ontology


def _save_obo_format(onto: Ontology, path: Path) -> None:
    """Save ontology in OBO format."""
    lines = []

    # Write header
    lines.append("format-version: 1.2")
    lines.append(f"data-version: {onto.metadata.get('data-version', 'unknown')}")
    lines.append(f"date: {onto.metadata.get('date', 'unknown')}")
    lines.append(f"saved-by: metainformant")
    lines.append("auto-generated-by: metainformant")

    # Add other metadata
    for key, value in onto.metadata.items():
        if key not in ['data-version', 'date']:
            lines.append(f"{key}: {value}")

    lines.append("")  # Empty line after header

    # Write terms
    for term in onto.terms.values():
        lines.append("[Term]")
        lines.append(f"id: {term.id}")

        if term.name:
            lines.append(f"name: {term.name}")

        if term.namespace:
            lines.append(f"namespace: {term.namespace}")

        if term.definition:
            lines.append(f"def: \"{term.definition}\"")

        if term.synonyms:
            for synonym in term.synonyms:
                lines.append(f"synonym: \"{synonym}\" EXACT []")

        if term.xrefs:
            for xref in term.xrefs:
                lines.append(f"xref: {xref}")

        if term.is_obsolete:
            lines.append("is_obsolete: true")

        # Add relationships from term metadata
        if 'relationships' in term.metadata:
            for rel_type, target in term.metadata['relationships']:
                lines.append(f"{rel_type}: {target}")

        lines.append("")  # Empty line after term

    # Write content to file
    content = "\n".join(lines)
    io.write_text(content, path)


def _save_json_format(onto: Ontology, path: Path) -> None:
    """Save ontology in JSON format."""
    # Convert ontology to serializable dictionary
    onto_dict = {
        'metadata': onto.metadata,
        'terms': {},
        'relationships': []
    }

    # Convert terms
    for term_id, term in onto.terms.items():
        term_dict = {
            'id': term.id,
            'name': term.name,
            'definition': term.definition,
            'namespace': term.namespace,
            'synonyms': term.synonyms,
            'xrefs': term.xrefs,
            'is_obsolete': term.is_obsolete,
            'metadata': term.metadata
        }
        onto_dict['terms'][term_id] = term_dict

    # Convert relationships
    for rel in onto.relationships:
        rel_dict = {
            'source': rel.source,
            'target': rel.target,
            'relation_type': rel.relation_type,
            'metadata': rel.metadata
        }
        onto_dict['relationships'].append(rel_dict)

    # Save as JSON
    io.dump_json(onto_dict, path, indent=2)


def _load_json_format(path: Path) -> Ontology:
    """Load ontology from JSON format."""
    # Load JSON data
    data = io.load_json(path)

    # Extract metadata
    metadata = data.get('metadata', {})

    # Reconstruct terms
    terms = {}
    for term_id, term_dict in data.get('terms', {}).items():
        from .types import create_term
        term = create_term(
            id=term_dict['id'],
            name=term_dict.get('name'),
            definition=term_dict.get('definition'),
            namespace=term_dict.get('namespace'),
            synonyms=term_dict.get('synonyms', []),
            xrefs=term_dict.get('xrefs', []),
            is_obsolete=term_dict.get('is_obsolete', False),
            **term_dict.get('metadata', {})
        )
        terms[term_id] = term

    # Reconstruct relationships
    relationships = []
    for rel_dict in data.get('relationships', []):
        from .types import create_relationship
        relationship = create_relationship(
            source=rel_dict['source'],
            target=rel_dict['target'],
            relation_type=rel_dict['relation_type'],
            **rel_dict.get('metadata', {})
        )
        relationships.append(relationship)

    # Create ontology
    from .types import create_ontology
    ontology = create_ontology(
        terms=terms,
        relationships=relationships,
        **metadata
    )

    return ontology


def export_ontology_stats(onto: Ontology, path: str | Path) -> None:
    """Export comprehensive ontology statistics to a file.

    Args:
        onto: Ontology to analyze.
        path: Path to save the statistics.

    Examples:
        >>> export_ontology_stats(ontology, "output/ontology_stats.json")
    """
    from .query import get_subontology_stats
    stats = get_subontology_stats(onto)

    # Add additional statistics
    stats['exported_by'] = 'metainformant'
    stats['export_timestamp'] = str(Path(path).stat().st_mtime) if path.exists() else None

    # Calculate more detailed stats
    term_name_lengths = [len(term.name) for term in onto.terms.values() if term.name]
    if term_name_lengths:
        stats['term_name_stats'] = {
            'mean_length': sum(term_name_lengths) / len(term_name_lengths),
            'max_length': max(term_name_lengths),
            'min_length': min(term_name_lengths)
        }

    # Relationship type distribution
    rel_type_counts = {}
    for rel in onto.relationships:
        rel_type_counts[rel.relation_type] = rel_type_counts.get(rel.relation_type, 0) + 1
    stats['relationship_type_distribution'] = rel_type_counts

    # Save statistics
    io.dump_json(stats, Path(path), indent=2)
    logger.info(f"Ontology statistics exported to {path}")


def merge_ontologies(*ontologies: Ontology, conflict_resolution: str = "first") -> Ontology:
    """Merge multiple ontologies into one.

    Args:
        *ontologies: Ontology objects to merge.
        conflict_resolution: How to handle term conflicts ("first", "last", "error").

    Returns:
        Merged Ontology object.

    Raises:
        ValueError: If conflict_resolution is invalid or conflicts occur when set to "error".

    Examples:
        >>> merged = merge_ontologies(onto1, onto2, onto3)
    """
    if not ontologies:
        from .types import create_ontology
        return create_ontology()

    if conflict_resolution not in ["first", "last", "error"]:
        raise ValueError("conflict_resolution must be 'first', 'last', or 'error'")

    # Start with first ontology
    merged_terms = dict(ontologies[0].terms)
    merged_relationships = list(ontologies[0].relationships)
    merged_metadata = dict(ontologies[0].metadata)

    # Merge remaining ontologies
    for onto in ontologies[1:]:
        # Merge terms
        for term_id, term in onto.terms.items():
            if term_id in merged_terms:
                if conflict_resolution == "error":
                    raise ValueError(f"Term conflict for {term_id}")
                elif conflict_resolution == "last":
                    merged_terms[term_id] = term
                # For "first", keep existing term
            else:
                merged_terms[term_id] = term

        # Merge relationships (avoid duplicates)
        existing_rels = {(r.source, r.target, r.relation_type) for r in merged_relationships}
        for rel in onto.relationships:
            rel_key = (rel.source, rel.target, rel.relation_type)
            if rel_key not in existing_rels:
                merged_relationships.append(rel)
                existing_rels.add(rel_key)

        # Merge metadata
        merged_metadata.update(onto.metadata)

    # Create merged ontology
    from .types import create_ontology
    merged_onto = create_ontology(
        terms=merged_terms,
        relationships=merged_relationships,
        **merged_metadata
    )

    logger.info(f"Merged {len(ontologies)} ontologies into one with {len(merged_onto)} terms")
    return merged_onto
