"""Pathway analysis and enrichment for biological networks.

This module provides tools for pathway enrichment analysis, pathway topology
analysis, and pathway visualization in biological networks.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Optional dependencies
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    logger.warning("networkx not available, pathway analysis disabled")

try:
    import scipy.stats as stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    logger.warning("scipy not available, statistical tests disabled")


def pathway_enrichment_analysis(
    genes: List[str],
    background_genes: List[str],
    pathways: Dict[str, List[str]],
    method: str = "fisher",
    correction: str = "bonferroni",
    **kwargs: Any
) -> List[Dict[str, Any]]:
    """Perform pathway enrichment analysis.

    Args:
        genes: List of genes of interest
        background_genes: List of background genes
        pathways: Dictionary mapping pathway names to gene lists
        method: Statistical test method ('fisher', 'hypergeometric')
        correction: Multiple testing correction ('bonferroni', 'fdr')
        **kwargs: Additional parameters

    Returns:
        List of enrichment results sorted by p-value

    Raises:
        ImportError: If scipy not available
    """
    if not HAS_SCIPY:
        raise ImportError("scipy required for statistical tests")

    if not HAS_NETWORKX:
        raise ImportError("networkx required for pathway analysis")

    results = []

    # Convert to sets for faster lookup
    gene_set = set(genes)
    background_set = set(background_genes)

    # Validate inputs
    if not gene_set.issubset(background_set):
        logger.warning("Some genes not found in background set")

    for pathway_name, pathway_genes in pathways.items():
        pathway_set = set(pathway_genes)

        # Calculate contingency table
        # [[genes in pathway, genes not in pathway],
        #  [background in pathway, background not in pathway]]
        genes_in_pathway = len(gene_set & pathway_set)
        genes_not_in_pathway = len(gene_set - pathway_set)
        background_in_pathway = len((background_set - gene_set) & pathway_set)
        background_not_in_pathway = len((background_set - gene_set) - pathway_set)

        # Fisher's exact test or hypergeometric test
        if method == "fisher":
            odds_ratio, p_value = stats.fisher_exact([
                [genes_in_pathway, genes_not_in_pathway],
                [background_in_pathway, background_not_in_pathway]
            ])
        elif method == "hypergeometric":
            # Hypergeometric test
            M = len(background_set & pathway_set)  # White balls in urn
            n = len(gene_set)  # Balls drawn
            k = genes_in_pathway  # White balls drawn

            p_value = stats.hypergeom.sf(k - 1, M, n, len(pathway_set))
            odds_ratio = (genes_in_pathway / max(genes_not_in_pathway, 1)) / \
                        (background_in_pathway / max(background_not_in_pathway, 1))
        else:
            raise ValueError(f"Unknown method: {method}")

        # Calculate enrichment metrics
        expected = (len(pathway_set) / len(background_set)) * len(genes)
        enrichment_ratio = genes_in_pathway / expected if expected > 0 else 1.0

        result = {
            'pathway': pathway_name,
            'p_value': p_value,
            'odds_ratio': odds_ratio,
            'enrichment_ratio': enrichment_ratio,
            'genes_in_pathway': genes_in_pathway,
            'pathway_size': len(pathway_set),
            'genes_tested': len(genes),
            'background_size': len(background_set),
            'method': method,
        }

        results.append(result)

    # Apply multiple testing correction
    p_values = [r['p_value'] for r in results]

    if correction == "bonferroni":
        corrected_p_values = [min(p * len(results), 1.0) for p in p_values]
    elif correction == "fdr":
        # Benjamini-Hochberg FDR correction
        sorted_indices = sorted(range(len(p_values)), key=lambda i: p_values[i])
        corrected_p_values = [1.0] * len(p_values)

        for i, idx in enumerate(sorted_indices):
            rank = i + 1
            corrected_p = min(p_values[idx] * len(p_values) / rank, 1.0)
            corrected_p_values[idx] = corrected_p

            # Ensure monotonicity
            if i > 0:
                corrected_p_values[idx] = min(corrected_p_values[idx],
                                            corrected_p_values[sorted_indices[i-1]])
    else:
        corrected_p_values = p_values

    # Add corrected p-values
    for result, corrected_p in zip(results, corrected_p_values):
        result['p_value_corrected'] = corrected_p
        result['significant'] = corrected_p < 0.05

    # Sort by corrected p-value
    results.sort(key=lambda x: x['p_value_corrected'])

    logger.info(f"Pathway enrichment analysis completed: {len(results)} pathways tested")
    return results


def pathway_topology_analysis(
    pathway_graph: Any,
    **kwargs: Any
) -> Dict[str, Any]:
    """Analyze topological properties of pathway graphs.

    Args:
        pathway_graph: NetworkX graph representing a pathway
        **kwargs: Additional analysis parameters

    Returns:
        Dictionary with topology analysis results

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for topology analysis")

    analysis = {
        'n_nodes': len(pathway_graph.nodes()),
        'n_edges': len(pathway_graph.edges()),
        'directed': pathway_graph.is_directed(),
    }

    if analysis['n_nodes'] == 0:
        return analysis

    # Degree analysis
    degrees = [d for n, d in pathway_graph.degree()]
    analysis['degree_stats'] = {
        'mean': float(sum(degrees) / len(degrees)),
        'max': max(degrees),
        'min': min(degrees),
    }

    # Clustering coefficient
    if not pathway_graph.is_directed():
        try:
            clustering = nx.average_clustering(pathway_graph)
            analysis['average_clustering'] = clustering
        except:
            analysis['average_clustering'] = None

    # Centrality measures
    try:
        betweenness = nx.betweenness_centrality(pathway_graph)
        analysis['betweenness_centrality'] = dict(betweenness)
    except:
        analysis['betweenness_centrality'] = {}

    # Connected components
    if not pathway_graph.is_directed():
        components = list(nx.connected_components(pathway_graph))
        analysis['connected_components'] = {
            'n_components': len(components),
            'sizes': [len(c) for c in components],
        }

    # Graph density
    max_edges = analysis['n_nodes'] * (analysis['n_nodes'] - 1)
    if pathway_graph.is_directed():
        max_edges *= 2
    analysis['density'] = analysis['n_edges'] / max_edges if max_edges > 0 else 0

    return analysis


def find_pathway_modules(
    pathway_graph: Any,
    method: str = "louvain",
    **kwargs: Any
) -> List[List[str]]:
    """Identify functional modules within pathways.

    Args:
        pathway_graph: NetworkX pathway graph
        method: Community detection method
        **kwargs: Parameters for community detection

    Returns:
        List of module (community) node lists

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for module detection")

    # Import community detection
    from metainformant.networks.community import detect_communities

    return detect_communities(pathway_graph, method=method, **kwargs)


def pathway_similarity_analysis(
    pathways: Dict[str, Any],
    method: str = "jaccard",
    **kwargs: Any
) -> Dict[str, Any]:
    """Analyze similarities between pathways.

    Args:
        pathways: Dictionary mapping pathway names to pathway data
        method: Similarity method ('jaccard', 'overlap', 'cosine')
        **kwargs: Additional parameters

    Returns:
        Dictionary with pathway similarities

    Raises:
        ValueError: If invalid method
    """
    if method not in ['jaccard', 'overlap', 'cosine']:
        raise ValueError(f"Unknown similarity method: {method}")

    pathway_names = list(pathways.keys())
    n_pathways = len(pathway_names)

    similarity_matrix = {}

    for i in range(n_pathways):
        for j in range(i + 1, n_pathways):
            name1, name2 = pathway_names[i], pathway_names[j]

            # Get gene sets for each pathway
            if isinstance(pathways[name1], dict) and 'genes' in pathways[name1]:
                genes1 = set(pathways[name1]['genes'])
            elif isinstance(pathways[name1], list):
                genes1 = set(pathways[name1])
            else:
                genes1 = set()

            if isinstance(pathways[name2], dict) and 'genes' in pathways[name2]:
                genes2 = set(pathways[name2]['genes'])
            elif isinstance(pathways[name2], list):
                genes2 = set(pathways[name2])
            else:
                genes2 = set()

            # Calculate similarity
            if method == "jaccard":
                intersection = len(genes1 & genes2)
                union = len(genes1 | genes2)
                similarity = intersection / union if union > 0 else 0.0

            elif method == "overlap":
                intersection = len(genes1 & genes2)
                min_size = min(len(genes1), len(genes2))
                similarity = intersection / min_size if min_size > 0 else 0.0

            elif method == "cosine":
                intersection = len(genes1 & genes2)
                similarity = intersection / (math.sqrt(len(genes1)) * math.sqrt(len(genes2)))

            similarity_matrix[f"{name1}_{name2}"] = similarity

    return {
        'similarities': similarity_matrix,
        'method': method,
        'n_pathways': n_pathways,
        'pathway_names': pathway_names,
    }


def pathway_hierarchy_analysis(
    pathways: Dict[str, Any],
    hierarchy_data: Optional[Dict[str, List[str]]] = None,
    **kwargs: Any
) -> Dict[str, Any]:
    """Analyze pathway hierarchies and relationships.

    Args:
        pathways: Dictionary of pathway data
        hierarchy_data: Optional hierarchy relationships
        **kwargs: Additional parameters

    Returns:
        Dictionary with hierarchy analysis
    """
    analysis = {
        'n_pathways': len(pathways),
        'hierarchy_levels': {},
    }

    if hierarchy_data:
        # Build hierarchy graph
        if HAS_NETWORKX:
            hierarchy_graph = nx.DiGraph()

            for parent, children in hierarchy_data.items():
                for child in children:
                    hierarchy_graph.add_edge(parent, child)

            # Calculate hierarchy statistics
            try:
                levels = {}
                for node in nx.topological_sort(hierarchy_graph):
                    predecessors = list(hierarchy_graph.predecessors(node))
                    if not predecessors:
                        level = 0  # Root level
                    else:
                        level = max(levels.get(pred, 0) for pred in predecessors) + 1
                    levels[node] = level

                analysis['hierarchy_levels'] = levels
                analysis['max_depth'] = max(levels.values()) if levels else 0
                analysis['n_roots'] = sum(1 for level in levels.values() if level == 0)

            except nx.NetworkXError:
                logger.warning("Could not perform topological sort on hierarchy")

    return analysis


def create_pathway_network(
    pathways: Dict[str, Any],
    similarity_threshold: float = 0.1,
    **kwargs: Any
) -> Any:
    """Create a network of pathway relationships.

    Args:
        pathways: Dictionary of pathway data
        similarity_threshold: Minimum similarity for edge creation
        **kwargs: Additional parameters

    Returns:
        NetworkX graph of pathway relationships

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for pathway network creation")

    # Calculate pathway similarities
    similarities = pathway_similarity_analysis(pathways, **kwargs)

    # Create network
    pathway_network = nx.Graph()

    # Add pathway nodes
    for pathway_name, pathway_data in pathways.items():
        pathway_network.add_node(pathway_name, **pathway_data)

    # Add similarity edges
    for key, similarity in similarities['similarities'].items():
        if similarity >= similarity_threshold:
            pathway1, pathway2 = key.split('_', 1)
            pathway_network.add_edge(pathway1, pathway2, similarity=similarity)

    logger.info(f"Created pathway network with {len(pathway_network.nodes())} nodes and {len(pathway_network.edges())} edges")
    return pathway_network


def pathway_disease_association(
    pathway_results: List[Dict[str, Any]],
    disease_genes: Dict[str, List[str]],
    **kwargs: Any
) -> Dict[str, Any]:
    """Associate pathways with diseases based on enrichment results.

    Args:
        pathway_results: Results from pathway enrichment analysis
        disease_genes: Dictionary mapping diseases to gene lists
        **kwargs: Additional parameters

    Returns:
        Dictionary with disease-pathway associations
    """
    associations = {}

    for disease, genes in disease_genes.items():
        disease_associations = []

        for pathway_result in pathway_results:
            pathway_name = pathway_result['pathway']
            enriched_genes = []  # Would need pathway gene lists

            # Calculate overlap
            pathway_gene_set = set()  # Would need to populate this
            disease_gene_set = set(genes)
            overlap = len(pathway_gene_set & disease_gene_set)

            if overlap > 0:
                disease_associations.append({
                    'pathway': pathway_name,
                    'overlap': overlap,
                    'p_value': pathway_result.get('p_value_corrected', 1.0),
                })

        if disease_associations:
            # Sort by significance
            disease_associations.sort(key=lambda x: x['p_value'])
            associations[disease] = disease_associations

    return associations


def pathway_visualization_data(
    pathway_graph: Any,
    layout_method: str = "spring",
    **kwargs: Any
) -> Dict[str, Any]:
    """Prepare pathway data for visualization.

    Args:
        pathway_graph: NetworkX pathway graph
        layout_method: Node layout method ('spring', 'circular', 'random')
        **kwargs: Additional visualization parameters

    Returns:
        Dictionary with visualization data

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for visualization data")

    # Calculate layout
    if layout_method == "spring":
        pos = nx.spring_layout(pathway_graph, **kwargs)
    elif layout_method == "circular":
        pos = nx.circular_layout(pathway_graph, **kwargs)
    elif layout_method == "random":
        pos = nx.random_layout(pathway_graph, **kwargs)
    else:
        raise ValueError(f"Unknown layout method: {layout_method}")

    # Prepare node and edge data
    nodes = []
    for node, attrs in pathway_graph.nodes(data=True):
        node_data = {
            'id': node,
            'x': pos[node][0],
            'y': pos[node][1],
        }
        node_data.update(attrs)
        nodes.append(node_data)

    edges = []
    for source, target, attrs in pathway_graph.edges(data=True):
        edge_data = {
            'source': source,
            'target': target,
        }
        edge_data.update(attrs)
        edges.append(edge_data)

    return {
        'nodes': nodes,
        'edges': edges,
        'layout_method': layout_method,
        'n_nodes': len(nodes),
        'n_edges': len(edges),
    }
