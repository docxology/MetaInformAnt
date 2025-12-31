"""Protein-protein interaction network analysis.

This module provides specialized tools for analyzing protein-protein interaction
(PPI) networks, including network construction, analysis, and visualization.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from metainformant.core import logging, io

logger = logging.get_logger(__name__)

# Optional dependencies
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    logger.warning("networkx not available, PPI analysis disabled")


def load_ppi_network(
    ppi_file: Union[str, Path],
    format: str = "tsv",
    **kwargs: Any
) -> Any:
    """Load PPI network from file.

    Args:
        ppi_file: Path to PPI data file
        format: File format ('tsv', 'csv', 'bioplex', 'intact')
        **kwargs: Additional loading parameters

    Returns:
        NetworkX graph representing PPI network

    Raises:
        ImportError: If networkx not available
        ValueError: If format not supported
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for PPI network loading")

    ppi_file = Path(ppi_file)

    if format == "tsv":
        # Assume tab-separated: protein1, protein2, score
        ppi_data = []
        with open(ppi_file, 'r') as f:
            for line_num, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) >= 2:
                    protein1, protein2 = parts[0], parts[1]
                    score = float(parts[2]) if len(parts) > 2 else 1.0
                    ppi_data.append((protein1, protein2, score))

    elif format == "csv":
        # CSV format
        import csv
        ppi_data = []
        with open(ppi_file, 'r') as f:
            reader = csv.reader(f)
            next(reader, None)  # Skip header
            for row in reader:
                if len(row) >= 2:
                    protein1, protein2 = row[0], row[1]
                    score = float(row[2]) if len(row) > 2 else 1.0
                    ppi_data.append((protein1, protein2, score))

    elif format == "bioplex":
        # BioPlex format (specific columns)
        ppi_data = []
        with open(ppi_file, 'r') as f:
            header = next(f).strip().split('\t')
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) > 1:
                    # BioPlex has specific column indices
                    protein1 = parts[0]
                    protein2 = parts[1]
                    score = 1.0  # BioPlex doesn't have scores
                    ppi_data.append((protein1, protein2, score))

    elif format == "intact":
        # IntAct PSI-MI TAB format
        ppi_data = []
        with open(ppi_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 6:
                    # IntAct format: ID interactor A, ID interactor B, etc.
                    protein1_id = parts[0]
                    protein2_id = parts[1]
                    # Extract protein names from IDs
                    protein1 = protein1_id.split(':')[-1] if ':' in protein1_id else protein1_id
                    protein2 = protein2_id.split(':')[-1] if ':' in protein2_id else protein2_id
                    score = 1.0
                    ppi_data.append((protein1, protein2, score))

    else:
        raise ValueError(f"Unsupported PPI format: {format}")

    # Create NetworkX graph
    G = nx.Graph()

    for protein1, protein2, score in ppi_data:
        if G.has_edge(protein1, protein2):
            # Update existing edge with higher score
            current_score = G[protein1][protein2].get('weight', 0)
            if score > current_score:
                G[protein1][protein2]['weight'] = score
        else:
            G.add_edge(protein1, protein2, weight=score)

    logger.info(f"Loaded PPI network: {len(G.nodes())} proteins, {len(G.edges())} interactions")
    return G


def construct_ppi_network_from_interactions(
    interactions: List[Tuple[str, str, float]],
    **kwargs: Any
) -> Any:
    """Construct PPI network from interaction list.

    Args:
        interactions: List of (protein1, protein2, score) tuples
        **kwargs: Additional graph parameters

    Returns:
        NetworkX PPI graph

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for PPI network construction")

    G = nx.Graph(**kwargs)

    for protein1, protein2, score in interactions:
        if G.has_edge(protein1, protein2):
            # Keep higher score
            current_score = G[protein1][protein2].get('weight', 0)
            if score > current_score:
                G[protein1][protein2]['weight'] = score
        else:
            G.add_edge(protein1, protein2, weight=score)

    logger.info(f"Constructed PPI network: {len(G.nodes())} proteins, {len(G.edges())} interactions")
    return G


def ppi_network_analysis(
    ppi_graph: Any,
    **kwargs: Any
) -> Dict[str, Any]:
    """Comprehensive analysis of PPI network properties.

    Args:
        ppi_graph: NetworkX PPI graph
        **kwargs: Additional analysis parameters

    Returns:
        Dictionary with PPI network analysis results

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for PPI analysis")

    analysis = {
        'basic_stats': {
            'n_proteins': len(ppi_graph.nodes()),
            'n_interactions': len(ppi_graph.edges()),
        }
    }

    if analysis['basic_stats']['n_proteins'] == 0:
        return analysis

    # Degree analysis
    degrees = [d for n, d in ppi_graph.degree()]
    analysis['degree_distribution'] = {
        'mean_degree': float(sum(degrees) / len(degrees)),
        'max_degree': max(degrees),
        'min_degree': min(degrees),
        'degrees': degrees,
    }

    # Connected components
    components = list(nx.connected_components(ppi_graph))
    component_sizes = [len(c) for c in components]

    analysis['connectivity'] = {
        'n_components': len(components),
        'largest_component_size': max(component_sizes) if component_sizes else 0,
        'component_sizes': component_sizes,
        'isolated_proteins': len([c for c in components if len(c) == 1]),
    }

    # Clustering coefficient
    try:
        avg_clustering = nx.average_clustering(ppi_graph)
        analysis['clustering'] = {
            'average_clustering_coefficient': avg_clustering,
        }
    except Exception as e:
        logger.warning(f"Could not calculate clustering coefficient: {e}")
        analysis['clustering'] = {'error': str(e)}

    # Network density
    n_nodes = len(ppi_graph.nodes())
    max_edges = n_nodes * (n_nodes - 1) / 2
    density = len(ppi_graph.edges()) / max_edges if max_edges > 0 else 0

    analysis['topology'] = {
        'density': density,
        'average_shortest_path_length': None,  # Too expensive for large networks
    }

    # Centrality measures (sample for large networks)
    max_nodes_for_centrality = kwargs.get('max_nodes_for_centrality', 1000)

    if n_nodes <= max_nodes_for_centrality:
        try:
            degree_cent = nx.degree_centrality(ppi_graph)
            betweenness_cent = nx.betweenness_centrality(ppi_graph)

            analysis['centrality'] = {
                'degree_centrality': dict(list(degree_cent.items())[:10]),  # Top 10
                'betweenness_centrality': dict(list(betweenness_cent.items())[:10]),
            }
        except Exception as e:
            logger.warning(f"Could not calculate centrality: {e}")
            analysis['centrality'] = {'error': str(e)}
    else:
        analysis['centrality'] = {
            'skipped': True,
            'reason': f'Network too large ({n_nodes} > {max_nodes_for_centrality} nodes)',
        }

    return analysis


def find_ppi_hubs(
    ppi_graph: Any,
    degree_threshold: Optional[int] = None,
    percentile: float = 95.0,
    **kwargs: Any
) -> List[Tuple[str, int]]:
    """Identify hub proteins in PPI network.

    Args:
        ppi_graph: NetworkX PPI graph
        degree_threshold: Minimum degree for hub classification
        percentile: Percentile threshold for hub identification
        **kwargs: Additional parameters

    Returns:
        List of (protein, degree) tuples for hub proteins

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for hub identification")

    # Calculate degrees
    degrees = dict(ppi_graph.degree())

    if degree_threshold is None:
        # Use percentile threshold
        degree_values = list(degrees.values())
        degree_threshold = int(np.percentile(degree_values, percentile))

    # Find hubs
    hubs = [(protein, degree) for protein, degree in degrees.items()
            if degree >= degree_threshold]

    # Sort by degree descending
    hubs.sort(key=lambda x: x[1], reverse=True)

    logger.info(f"Identified {len(hubs)} hub proteins (degree ≥ {degree_threshold})")
    return hubs


def ppi_network_clustering(
    ppi_graph: Any,
    method: str = "louvain",
    **kwargs: Any
) -> List[List[str]]:
    """Cluster PPI network into functional modules.

    Args:
        ppi_graph: NetworkX PPI graph
        method: Clustering method
        **kwargs: Additional clustering parameters

    Returns:
        List of protein clusters

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for PPI clustering")

    # Import community detection
    from metainformant.networks.community import detect_communities

    return detect_communities(ppi_graph, method=method, **kwargs)


def analyze_ppi_disease_associations(
    ppi_graph: Any,
    disease_proteins: Dict[str, List[str]],
    **kwargs: Any
) -> Dict[str, Any]:
    """Analyze disease associations in PPI network.

    Args:
        ppi_graph: NetworkX PPI graph
        disease_proteins: Dictionary mapping diseases to protein lists
        **kwargs: Additional analysis parameters

    Returns:
        Dictionary with disease association analysis

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for disease association analysis")

    associations = {}

    for disease, proteins in disease_proteins.items():
        disease_proteins_in_network = [p for p in proteins if p in ppi_graph.nodes()]

        if not disease_proteins_in_network:
            continue

        # Create disease subgraph
        disease_subgraph = ppi_graph.subgraph(disease_proteins_in_network)

        # Analyze disease module properties
        associations[disease] = {
            'proteins_in_network': len(disease_proteins_in_network),
            'total_proteins': len(proteins),
            'subgraph_edges': len(disease_subgraph.edges()),
            'subgraph_density': (
                len(disease_subgraph.edges()) /
                (len(disease_proteins_in_network) * (len(disease_proteins_in_network) - 1) / 2)
                if len(disease_proteins_in_network) > 1 else 0
            ),
        }

        # Find connections to other proteins
        disease_neighbors = set()
        for protein in disease_proteins_in_network:
            disease_neighbors.update(ppi_graph.neighbors(protein))
        disease_neighbors -= set(disease_proteins_in_network)

        associations[disease]['external_connections'] = len(disease_neighbors)

    return associations


def ppi_network_comparison(
    ppi_graph1: Any,
    ppi_graph2: Any,
    **kwargs: Any
) -> Dict[str, Any]:
    """Compare two PPI networks.

    Args:
        ppi_graph1: First PPI network
        ppi_graph2: Second PPI network
        **kwargs: Additional comparison parameters

    Returns:
        Dictionary with network comparison results

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network comparison")

    comparison = {
        'network1_stats': ppi_network_analysis(ppi_graph1),
        'network2_stats': ppi_network_analysis(ppi_graph2),
    }

    # Common proteins
    proteins1 = set(ppi_graph1.nodes())
    proteins2 = set(ppi_graph2.nodes())
    common_proteins = proteins1 & proteins2

    comparison['protein_overlap'] = {
        'common_proteins': len(common_proteins),
        'unique_to_1': len(proteins1 - proteins2),
        'unique_to_2': len(proteins2 - proteins1),
        'jaccard_similarity': len(common_proteins) / len(proteins1 | proteins2),
    }

    # Common interactions
    edges1 = set(ppi_graph1.edges())
    edges2 = set(ppi_graph2.edges())
    common_edges = edges1 & edges2

    comparison['interaction_overlap'] = {
        'common_interactions': len(common_edges),
        'unique_to_1': len(edges1 - edges2),
        'unique_to_2': len(edges2 - edges1),
        'jaccard_similarity': len(common_edges) / len(edges1 | edges2),
    }

    return comparison


def ppi_network_enrichment(
    protein_list: List[str],
    ppi_graph: Any,
    background_proteins: Optional[List[str]] = None,
    **kwargs: Any
) -> Dict[str, Any]:
    """Test for enrichment of protein interactions.

    Args:
        protein_list: List of proteins to test
        ppi_graph: PPI network graph
        background_proteins: Background protein set
        **kwargs: Additional enrichment parameters

    Returns:
        Dictionary with enrichment analysis results

    Raises:
        ImportError: If networkx or scipy not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for PPI enrichment")

    if background_proteins is None:
        background_proteins = list(ppi_graph.nodes())

    # Filter to proteins in network
    test_proteins = [p for p in protein_list if p in ppi_graph.nodes()]
    background_proteins = [p for p in background_proteins if p in ppi_graph.nodes()]

    if not test_proteins:
        return {'error': 'No test proteins found in PPI network'}

    # Count interactions within test set
    test_subgraph = ppi_graph.subgraph(test_proteins)
    observed_interactions = len(test_subgraph.edges())

    # Expected interactions under random model
    n_test = len(test_proteins)
    n_background = len(background_proteins)
    total_possible_interactions = n_background * (n_background - 1) / 2
    observed_background_interactions = len(ppi_graph.edges())

    expected_interactions = (
        observed_background_interactions *
        (n_test * (n_test - 1) / 2) /
        total_possible_interactions
    )

    # Calculate enrichment statistics
    enrichment_ratio = observed_interactions / expected_interactions if expected_interactions > 0 else 1.0

    # Statistical significance (approximate)
    try:
        import scipy.stats as stats
        # Use Poisson approximation for large numbers
        p_value = 1 - stats.poisson.cdf(observed_interactions - 1, expected_interactions)
    except ImportError:
        p_value = None

    return {
        'observed_interactions': observed_interactions,
        'expected_interactions': expected_interactions,
        'enrichment_ratio': enrichment_ratio,
        'p_value': p_value,
        'test_proteins': len(test_proteins),
        'background_proteins': len(background_proteins),
    }


def save_ppi_network(
    ppi_graph: Any,
    output_file: Union[str, Path],
    format: str = "tsv",
    **kwargs: Any
) -> None:
    """Save PPI network to file.

    Args:
        ppi_graph: NetworkX PPI graph
        output_file: Output file path
        format: Output format ('tsv', 'csv', 'json')
        **kwargs: Additional saving parameters

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for PPI saving")

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if format == "tsv":
        with open(output_file, 'w') as f:
            f.write("#protein1\tprotein2\tscore\n")
            for u, v, data in ppi_graph.edges(data=True):
                score = data.get('weight', 1.0)
                f.write(f"{u}\t{v}\t{score}\n")

    elif format == "csv":
        import csv
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['protein1', 'protein2', 'score'])
            for u, v, data in ppi_graph.edges(data=True):
                score = data.get('weight', 1.0)
                writer.writerow([u, v, score])

    elif format == "json":
        # Save as node-link format
        data = nx.node_link_data(ppi_graph)
        io.dump_json(data, output_file)

    else:
        raise ValueError(f"Unsupported output format: {format}")

    logger.info(f"Saved PPI network to {output_file} ({format} format)")
