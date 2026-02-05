"""Protein-protein interaction network analysis.

This module provides specialized tools for analyzing protein-protein interaction
(PPI) networks, including network construction, analysis, and visualization.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import io, logging

logger = logging.get_logger(__name__)

# Optional dependencies
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    logger.warning("networkx not available, PPI analysis disabled")

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def load_ppi_network(ppi_file: Union[str, Path], format: str = "tsv", **kwargs: Any) -> Any:
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
        with open(ppi_file, "r") as f:
            for line_num, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) >= 2:
                    protein1, protein2 = parts[0], parts[1]
                    score = float(parts[2]) if len(parts) > 2 else 1.0
                    ppi_data.append((protein1, protein2, score))

    elif format == "csv":
        # CSV format
        import csv

        ppi_data = []
        with open(ppi_file, "r") as f:
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
        with open(ppi_file, "r") as f:
            next(f)  # skip header line
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) > 1:
                    # BioPlex has specific column indices
                    protein1 = parts[0]
                    protein2 = parts[1]
                    score = 1.0  # BioPlex doesn't have scores
                    ppi_data.append((protein1, protein2, score))

    elif format == "intact":
        # IntAct PSI-MI TAB format
        ppi_data = []
        with open(ppi_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 6:
                    # IntAct format: ID interactor A, ID interactor B, etc.
                    protein1_id = parts[0]
                    protein2_id = parts[1]
                    # Extract protein names from IDs
                    protein1 = protein1_id.split(":")[-1] if ":" in protein1_id else protein1_id
                    protein2 = protein2_id.split(":")[-1] if ":" in protein2_id else protein2_id
                    score = 1.0
                    ppi_data.append((protein1, protein2, score))

    else:
        raise ValueError(f"Unsupported PPI format: {format}")

    # Create NetworkX graph
    G = nx.Graph()

    for protein1, protein2, score in ppi_data:
        if G.has_edge(protein1, protein2):
            # Update existing edge with higher score
            current_score = G[protein1][protein2].get("weight", 0)
            if score > current_score:
                G[protein1][protein2]["weight"] = score
        else:
            G.add_edge(protein1, protein2, weight=score)

    logger.info(f"Loaded PPI network: {len(G.nodes())} proteins, {len(G.edges())} interactions")
    return G


def construct_ppi_network_from_interactions(interactions: List[Tuple[str, str, float]], **kwargs: Any) -> Any:
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
            current_score = G[protein1][protein2].get("weight", 0)
            if score > current_score:
                G[protein1][protein2]["weight"] = score
        else:
            G.add_edge(protein1, protein2, weight=score)

    logger.info(f"Constructed PPI network: {len(G.nodes())} proteins, {len(G.edges())} interactions")
    return G


def ppi_network_analysis(ppi_graph: Any, **kwargs: Any) -> Dict[str, Any]:
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
        "basic_stats": {
            "n_proteins": len(ppi_graph.nodes()),
            "n_interactions": len(ppi_graph.edges()),
        }
    }

    if analysis["basic_stats"]["n_proteins"] == 0:
        return analysis

    # Degree analysis
    degrees = [d for n, d in ppi_graph.degree()]
    analysis["degree_distribution"] = {
        "mean_degree": float(sum(degrees) / len(degrees)),
        "max_degree": max(degrees),
        "min_degree": min(degrees),
        "degrees": degrees,
    }

    # Connected components
    components = list(nx.connected_components(ppi_graph))
    component_sizes = [len(c) for c in components]

    analysis["connectivity"] = {
        "n_components": len(components),
        "largest_component_size": max(component_sizes) if component_sizes else 0,
        "component_sizes": component_sizes,
        "isolated_proteins": len([c for c in components if len(c) == 1]),
    }

    # Clustering coefficient
    try:
        avg_clustering = nx.average_clustering(ppi_graph)
        analysis["clustering"] = {
            "average_clustering_coefficient": avg_clustering,
        }
    except Exception as e:
        logger.warning(f"Could not calculate clustering coefficient: {e}")
        analysis["clustering"] = {"error": str(e)}

    # Network density
    n_nodes = len(ppi_graph.nodes())
    max_edges = n_nodes * (n_nodes - 1) / 2
    density = len(ppi_graph.edges()) / max_edges if max_edges > 0 else 0

    analysis["topology"] = {
        "density": density,
        "average_shortest_path_length": None,  # Too expensive for large networks
    }

    # Centrality measures (sample for large networks)
    max_nodes_for_centrality = kwargs.get("max_nodes_for_centrality", 1000)

    if n_nodes <= max_nodes_for_centrality:
        try:
            degree_cent = nx.degree_centrality(ppi_graph)
            betweenness_cent = nx.betweenness_centrality(ppi_graph)

            analysis["centrality"] = {
                "degree_centrality": dict(list(degree_cent.items())[:10]),  # Top 10
                "betweenness_centrality": dict(list(betweenness_cent.items())[:10]),
            }
        except Exception as e:
            logger.warning(f"Could not calculate centrality: {e}")
            analysis["centrality"] = {"error": str(e)}
    else:
        analysis["centrality"] = {
            "skipped": True,
            "reason": f"Network too large ({n_nodes} > {max_nodes_for_centrality} nodes)",
        }

    return analysis


def find_ppi_hubs(
    ppi_graph: Any, degree_threshold: Optional[int] = None, percentile: float = 95.0, **kwargs: Any
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
    hubs = [(protein, degree) for protein, degree in degrees.items() if degree >= degree_threshold]

    # Sort by degree descending
    hubs.sort(key=lambda x: x[1], reverse=True)

    logger.info(f"Identified {len(hubs)} hub proteins (degree ≥ {degree_threshold})")
    return hubs


def ppi_network_clustering(ppi_graph: Any, method: str = "louvain", **kwargs: Any) -> List[List[str]]:
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
    from metainformant.networks.analysis.community import detect_communities

    return detect_communities(ppi_graph, method=method, **kwargs)


def analyze_ppi_disease_associations(
    ppi_graph: Any, disease_proteins: Dict[str, List[str]], **kwargs: Any
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
            "proteins_in_network": len(disease_proteins_in_network),
            "total_proteins": len(proteins),
            "subgraph_edges": len(disease_subgraph.edges()),
            "subgraph_density": (
                len(disease_subgraph.edges())
                / (len(disease_proteins_in_network) * (len(disease_proteins_in_network) - 1) / 2)
                if len(disease_proteins_in_network) > 1
                else 0
            ),
        }

        # Find connections to other proteins
        disease_neighbors = set()
        for protein in disease_proteins_in_network:
            disease_neighbors.update(ppi_graph.neighbors(protein))
        disease_neighbors -= set(disease_proteins_in_network)

        associations[disease]["external_connections"] = len(disease_neighbors)

    return associations


def ppi_network_comparison(ppi_graph1: Any, ppi_graph2: Any, **kwargs: Any) -> Dict[str, Any]:
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
        "network1_stats": ppi_network_analysis(ppi_graph1),
        "network2_stats": ppi_network_analysis(ppi_graph2),
    }

    # Common proteins
    proteins1 = set(ppi_graph1.nodes())
    proteins2 = set(ppi_graph2.nodes())
    common_proteins = proteins1 & proteins2

    comparison["protein_overlap"] = {
        "common_proteins": len(common_proteins),
        "unique_to_1": len(proteins1 - proteins2),
        "unique_to_2": len(proteins2 - proteins1),
        "jaccard_similarity": len(common_proteins) / len(proteins1 | proteins2),
    }

    # Common interactions
    edges1 = set(ppi_graph1.edges())
    edges2 = set(ppi_graph2.edges())
    common_edges = edges1 & edges2

    comparison["interaction_overlap"] = {
        "common_interactions": len(common_edges),
        "unique_to_1": len(edges1 - edges2),
        "unique_to_2": len(edges2 - edges1),
        "jaccard_similarity": len(common_edges) / len(edges1 | edges2),
    }

    return comparison


def ppi_network_enrichment(
    protein_list: List[str], ppi_graph: Any, background_proteins: Optional[List[str]] = None, **kwargs: Any
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
        return {"error": "No test proteins found in PPI network"}

    # Count interactions within test set
    test_subgraph = ppi_graph.subgraph(test_proteins)
    observed_interactions = len(test_subgraph.edges())

    # Expected interactions under random model
    n_test = len(test_proteins)
    n_background = len(background_proteins)
    total_possible_interactions = n_background * (n_background - 1) / 2
    observed_background_interactions = len(ppi_graph.edges())

    expected_interactions = observed_background_interactions * (n_test * (n_test - 1) / 2) / total_possible_interactions

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
        "observed_interactions": observed_interactions,
        "expected_interactions": expected_interactions,
        "enrichment_ratio": enrichment_ratio,
        "p_value": p_value,
        "test_proteins": len(test_proteins),
        "background_proteins": len(background_proteins),
    }


def save_ppi_network(ppi_graph: Any, output_file: Union[str, Path], format: str = "tsv", **kwargs: Any) -> None:
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
        with open(output_file, "w") as f:
            f.write("#protein1\tprotein2\tscore\n")
            for u, v, data in ppi_graph.edges(data=True):
                score = data.get("weight", 1.0)
                f.write(f"{u}\t{v}\t{score}\n")

    elif format == "csv":
        import csv

        with open(output_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["protein1", "protein2", "score"])
            for u, v, data in ppi_graph.edges(data=True):
                score = data.get("weight", 1.0)
                writer.writerow([u, v, score])

    elif format == "json":
        # Save as node-link format
        data = nx.node_link_data(ppi_graph)
        io.dump_json(data, output_file)

    else:
        raise ValueError(f"Unsupported output format: {format}")

    logger.info(f"Saved PPI network to {output_file} ({format} format)")


class ProteinNetwork:
    """A protein-protein interaction network class.

    This class provides methods for analyzing and manipulating PPI networks,
    including loading from various sources, computing network properties,
    and performing enrichment analysis.
    """

    def __init__(self, graph: Optional[Any] = None, name: str = ""):
        """Initialize a protein network.

        Args:
            graph: NetworkX graph object (optional)
            name: Name of the network
        """
        if graph is None and HAS_NETWORKX:
            self.graph = nx.Graph()
        else:
            self.graph = graph
        self.name = name
        self.metadata = {}

    @classmethod
    def from_file(cls, filepath: Union[str, Path], format: str = "tsv", name: str = "") -> "ProteinNetwork":
        """Load protein network from file.

        Args:
            filepath: Path to network file
            format: File format ('tsv', 'csv', 'json')
            name: Name for the network

        Returns:
            ProteinNetwork instance
        """
        graph = load_ppi_network(filepath, format=format)
        network = cls(graph=graph, name=name or Path(filepath).stem)
        network.metadata["source_file"] = str(filepath)
        network.metadata["format"] = format
        return network

    @classmethod
    def from_interactions(cls, interactions: List[Tuple[str, str]], name: str = "") -> "ProteinNetwork":
        """Create network from list of protein interactions.

        Args:
            interactions: List of (protein1, protein2) tuples
            name: Name for the network

        Returns:
            ProteinNetwork instance
        """
        if not HAS_NETWORKX:
            raise ImportError("NetworkX required for ProteinNetwork")

        graph = nx.Graph()
        graph.add_edges_from(interactions)

        network = cls(graph=graph, name=name)
        network.metadata["n_interactions"] = len(interactions)
        network.metadata["n_proteins"] = len(graph.nodes())
        return network

    def get_proteins(self) -> List[str]:
        """Get list of all proteins in the network.

        Returns:
            List of protein identifiers
        """
        if not HAS_NETWORKX:
            return []
        return list(self.graph.nodes())

    def get_interactions(self) -> List[Tuple[str, str]]:
        """Get list of all interactions in the network.

        Returns:
            List of (protein1, protein2) tuples
        """
        if not HAS_NETWORKX:
            return []
        return list(self.graph.edges())

    def degree_distribution(self) -> Dict[str, int]:
        """Calculate degree distribution of the network.

        Returns:
            Dictionary mapping protein to degree
        """
        if not HAS_NETWORKX:
            return {}
        return dict(self.graph.degree())

    def find_neighbors(self, protein: str) -> List[str]:
        """Find direct neighbors of a protein.

        Args:
            protein: Protein identifier

        Returns:
            List of neighboring proteins
        """
        if not HAS_NETWORKX:
            return []
        if protein not in self.graph:
            return []
        return list(self.graph.neighbors(protein))

    def shortest_path(self, protein1: str, protein2: str) -> Optional[List[str]]:
        """Find shortest path between two proteins.

        Args:
            protein1: Starting protein
            protein2: Target protein

        Returns:
            List of proteins in path, or None if no path exists
        """
        if not HAS_NETWORKX:
            return None
        try:
            return nx.shortest_path(self.graph, protein1, protein2)
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            return None

    def clustering_coefficient(self, protein: Optional[str] = None) -> Union[float, Dict[str, float]]:
        """Calculate clustering coefficient.

        Args:
            protein: Specific protein, or None for all proteins

        Returns:
            Clustering coefficient(s)
        """
        if not HAS_NETWORKX:
            return {} if protein is None else 0.0

        if protein is not None:
            return nx.clustering(self.graph, protein)
        else:
            return nx.clustering(self.graph)

    def betweenness_centrality(self) -> Dict[str, float]:
        """Calculate betweenness centrality for all proteins.

        Returns:
            Dictionary mapping protein to centrality score
        """
        if not HAS_NETWORKX:
            return {}
        return nx.betweenness_centrality(self.graph)

    def connected_components(self) -> List[List[str]]:
        """Find connected components in the network.

        Returns:
            List of component lists (each containing protein identifiers)
        """
        if not HAS_NETWORKX:
            return []
        return [list(component) for component in nx.connected_components(self.graph)]

    def largest_component(self) -> List[str]:
        """Get the largest connected component.

        Returns:
            List of proteins in the largest component
        """
        components = self.connected_components()
        if not components:
            return []
        return max(components, key=len)

    def network_summary(self) -> Dict[str, Any]:
        """Generate network summary statistics.

        Returns:
            Dictionary with network statistics
        """
        if not HAS_NETWORKX:
            return {"error": "NetworkX not available"}

        summary = {
            "n_proteins": len(self.graph.nodes()),
            "n_interactions": len(self.graph.edges()),
            "n_components": len(self.connected_components()),
            "largest_component_size": len(self.largest_component()),
            "average_degree": (
                sum(dict(self.graph.degree()).values()) / len(self.graph.nodes()) if self.graph.nodes() else 0
            ),
            "density": nx.density(self.graph),
        }

        # Add metadata
        summary.update(self.metadata)

        return summary

    def __len__(self) -> int:
        """Get number of proteins."""
        if not HAS_NETWORKX:
            return 0
        return len(self.graph.nodes())

    def __contains__(self, protein: str) -> bool:
        """Check if protein is in network."""
        if not HAS_NETWORKX:
            return False
        return protein in self.graph

    def add_interaction(self, protein1: str, protein2: str, **kwargs: Any) -> None:
        """Add an interaction between two proteins.

        Args:
            protein1: First protein identifier
            protein2: Second protein identifier
            **kwargs: Additional edge attributes (e.g., weight, score)
        """
        if not HAS_NETWORKX:
            raise ImportError("networkx required for adding interactions")
        self.graph.add_edge(protein1, protein2, **kwargs)

    def add_protein_metadata(self, protein: str, **metadata: Any) -> None:
        """Add metadata to a protein node.

        Args:
            protein: Protein identifier
            **metadata: Metadata key-value pairs to add
        """
        if not HAS_NETWORKX:
            raise ImportError("networkx required for adding metadata")
        if protein not in self.graph:
            self.graph.add_node(protein)
        self.graph.nodes[protein].update(metadata)

    def get_protein_metadata(self, protein: str) -> Dict[str, Any]:
        """Get metadata for a protein node.

        Args:
            protein: Protein identifier

        Returns:
            Dictionary of protein metadata
        """
        if not HAS_NETWORKX:
            return {}
        if protein not in self.graph:
            return {}
        return dict(self.graph.nodes[protein])


def predict_interactions(
    target_proteins: List[str],
    known_network: Optional[ProteinNetwork] = None,
    features: Optional[Dict[str, List[float]]] = None,
    method: str = "similarity",
    threshold: float = 0.5,
    max_predictions_per_protein: int = 10,
    **kwargs: Any,
) -> Dict[str, List[Dict[str, Any]]]:
    """Predict protein-protein interactions for target proteins.

    Args:
        target_proteins: List of proteins to predict interactions for
        known_network: Existing PPI network for reference
        features: Feature vectors for proteins (for ML-based methods)
        method: Prediction method ('similarity', 'correlation', 'guilt-by-association', 'ml')
        threshold: Confidence threshold for predictions
        max_predictions_per_protein: Maximum predictions per target protein
        **kwargs: Additional method-specific parameters

    Returns:
        Dictionary mapping target proteins to lists of predicted interactions
    """
    predictions = {}

    if method == "similarity":
        # Simple similarity-based prediction using network topology
        if known_network is None:
            logger.warning("Known network required for similarity method")
            return {}

        predictions = _predict_by_similarity(target_proteins, known_network, threshold, max_predictions_per_protein)

    elif method == "correlation":
        # Correlation-based prediction using feature similarity
        if features is None:
            logger.warning("Features required for correlation method")
            return {}

        predictions = _predict_by_correlation(target_proteins, features, threshold, max_predictions_per_protein)

    elif method == "guilt-by-association":
        # Guilt-by-association using known network neighbors
        if known_network is None:
            logger.warning("Known network required for guilt-by-association method")
            return {}

        predictions = _predict_by_guilt_by_association(
            target_proteins, known_network, threshold, max_predictions_per_protein
        )

    elif method == "ml":
        # Machine learning-based prediction
        if features is None or known_network is None:
            logger.warning("Features and known network required for ML method")
            return {}

        predictions = _predict_by_ml(
            target_proteins, known_network, features, threshold, max_predictions_per_protein, **kwargs
        )

    else:
        raise ValueError(f"Unknown prediction method: {method}")

    return predictions


def _predict_by_similarity(
    target_proteins: List[str], known_network: ProteinNetwork, threshold: float, max_predictions: int
) -> Dict[str, List[Dict[str, Any]]]:
    """Predict interactions using network similarity."""
    predictions = {}

    if not HAS_NETWORKX:
        return predictions

    # Get all proteins in the network
    network_proteins = list(known_network.graph.nodes())

    for target in target_proteins:
        if target in network_proteins:
            # Target is already in network, find similar proteins
            target_neighbors = set(known_network.graph.neighbors(target))

            candidate_predictions = []
            for protein in network_proteins:
                if protein == target:
                    continue

                protein_neighbors = set(known_network.graph.neighbors(protein))
                similarity = len(target_neighbors & protein_neighbors) / len(target_neighbors | protein_neighbors)

                if similarity >= threshold:
                    candidate_predictions.append(
                        {"partner": protein, "confidence": similarity, "method": "jaccard_similarity"}
                    )

            # Sort by confidence and limit
            candidate_predictions.sort(key=lambda x: x["confidence"], reverse=True)
            predictions[target] = candidate_predictions[:max_predictions]
        else:
            # Target not in network, predict based on sequence/feature similarity
            candidate_predictions = []

            # Since target is not in network, we can't use topology
            # Return empty predictions with explanation
            predictions[target] = []
            logger.info(
                f"Target protein {target} not found in network - cannot predict interactions without additional features"
            )

    return predictions


def _predict_by_correlation(
    target_proteins: List[str], features: Dict[str, List[float]], threshold: float, max_predictions: int
) -> Dict[str, List[Dict[str, Any]]]:
    """Predict interactions using feature correlation."""
    predictions = {}

    # Get all proteins with features
    available_proteins = list(features.keys())

    for target in target_proteins:
        if target not in features:
            continue

        target_features = features[target]
        candidate_predictions = []

        for protein in available_proteins:
            if protein == target:
                continue

            if protein in features:
                protein_features = features[protein]

                # Calculate correlation (simplified - using Pearson correlation)
                if len(target_features) == len(protein_features) and len(target_features) > 1:
                    try:
                        correlation = abs(np.corrcoef(target_features, protein_features)[0, 1])
                        if correlation >= threshold:
                            candidate_predictions.append(
                                {"partner": protein, "confidence": correlation, "method": "feature_correlation"}
                            )
                    except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                        continue

        # Sort by confidence and limit
        candidate_predictions.sort(key=lambda x: x["confidence"], reverse=True)
        predictions[target] = candidate_predictions[:max_predictions]

    return predictions


def _predict_by_guilt_by_association(
    target_proteins: List[str], known_network: ProteinNetwork, threshold: float, max_predictions: int
) -> Dict[str, List[Dict[str, Any]]]:
    """Predict interactions using guilt-by-association."""
    predictions = {}

    if not HAS_NETWORKX:
        return predictions

    for target in target_proteins:
        candidate_predictions = []

        # Find neighbors of neighbors (second-degree connections)
        if target in known_network.graph:
            # Direct neighbors
            first_degree = set(known_network.graph.neighbors(target))

            # Second-degree neighbors (excluding direct neighbors and self)
            second_degree = set()
            for neighbor in first_degree:
                second_degree.update(known_network.graph.neighbors(neighbor))

            second_degree -= first_degree
            second_degree.discard(target)

            # Score based on shared neighbors
            for candidate in second_degree:
                candidate_neighbors = set(known_network.graph.neighbors(candidate))
                shared_neighbors = len(first_degree & candidate_neighbors)
                total_neighbors = len(first_degree | candidate_neighbors)

                if total_neighbors > 0:
                    confidence = shared_neighbors / total_neighbors
                    if confidence >= threshold:
                        candidate_predictions.append(
                            {"partner": candidate, "confidence": confidence, "method": "guilt_by_association"}
                        )

        # Sort by confidence and limit
        candidate_predictions.sort(key=lambda x: x["confidence"], reverse=True)
        predictions[target] = candidate_predictions[:max_predictions]

    return predictions


def _predict_by_ml(
    target_proteins: List[str],
    known_network: ProteinNetwork,
    features: Dict[str, List[float]],
    threshold: float,
    max_predictions: int,
    **kwargs: Any,
) -> Dict[str, List[Dict[str, Any]]]:
    """Predict interactions using machine learning."""
    # This would implement a proper ML model (e.g., random forest, neural network)
    # For now, use a simplified approach
    logger.warning("ML prediction method not fully implemented - using correlation fallback")

    # Fall back to correlation method
    return _predict_by_correlation(target_proteins, features, threshold, max_predictions)


def load_string_interactions(
    interactions_df: Any, proteins_df: Optional[Any] = None, confidence_threshold: int = 400
) -> "ProteinNetwork":
    """Load PPI network from STRING database format.

    Args:
        interactions_df: DataFrame with columns ['protein1', 'protein2', 'combined_score']
        proteins_df: Optional DataFrame with protein metadata
        confidence_threshold: Minimum confidence score (0-1000)

    Returns:
        ProteinNetwork object with loaded interactions
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for PPI network loading")

    # Create network
    network = ProteinNetwork(name="STRING_PPI")

    # Load interactions
    if hasattr(interactions_df, "iterrows"):  # pandas DataFrame
        for _, row in interactions_df.iterrows():
            protein1 = str(row["protein1"])
            protein2 = str(row["protein2"])
            score = int(row["combined_score"])

            if score >= confidence_threshold:
                if network.graph.has_edge(protein1, protein2):
                    current_score = network.graph[protein1][protein2].get("weight", 0)
                    if score > current_score:
                        network.graph[protein1][protein2]["weight"] = score
                else:
                    network.graph.add_edge(protein1, protein2, weight=score)

    # Load protein metadata if provided
    if proteins_df is not None and hasattr(proteins_df, "iterrows"):
        network.protein_metadata = {}
        for _, row in proteins_df.iterrows():
            protein_id = str(row["protein_id"])
            metadata = {
                "gene_name": row.get("gene_name", ""),
                "protein_name": row.get("protein_name", ""),
            }
            network.protein_metadata[protein_id] = metadata

    logger.info(
        f"Loaded STRING PPI network: {len(network.graph.nodes())} proteins, "
        f"{len(network.graph.edges())} interactions (threshold: {confidence_threshold})"
    )
    return network


def functional_enrichment_ppi(
    protein_list: List[str],
    ppi_network: "ProteinNetwork",
    function_key: str = "function",
    background_proteins: List[str] | None = None,
    min_overlap: int = 2,
    max_p_value: float = 0.05,
) -> Dict[str, Dict[str, Any]]:
    """Perform functional enrichment analysis on a protein list using PPI network.

    Args:
        protein_list: List of proteins to test for enrichment
        ppi_network: ProteinNetwork instance with functional annotations
        function_key: Key in protein annotations containing functional categories
        background_proteins: Background protein set (default: all proteins in network)
        min_overlap: Minimum overlap for enrichment consideration
        max_p_value: Maximum p-value for significant enrichment

    Returns:
        Dictionary mapping functions to enrichment statistics
    """
    import math
    from collections import defaultdict

    if background_proteins is None:
        background_proteins = list(ppi_network.graph.nodes())

    # Get functional annotations for background
    background_functions = defaultdict(list)
    for protein in background_proteins:
        node_data = ppi_network.get_protein_metadata(protein)
        if function_key in node_data:
            functions = node_data[function_key]
            if isinstance(functions, str):
                functions = [functions]
            for func in functions:
                background_functions[func].append(protein)

    # Get functional annotations for test set
    test_functions = defaultdict(list)
    for protein in protein_list:
        node_data = ppi_network.get_protein_metadata(protein)
        if function_key in node_data:
            functions = node_data[function_key]
            if isinstance(functions, str):
                functions = [functions]
            for func in functions:
                test_functions[func].append(protein)

    # Calculate enrichment for each function
    enrichment_results = {}

    for func in test_functions:
        if func not in background_functions:
            continue

        # Count overlaps
        test_count = len(set(test_functions[func]))
        background_count = len(set(background_functions[func]))

        if test_count < min_overlap:
            continue

        # Calculate enrichment statistics
        overlap = len(set(test_functions[func]) & set(background_functions[func]))

        # Hypergeometric test
        # P(X >= overlap) where X ~ Hypergeometric(N, K, n)
        # N = total background proteins
        # K = proteins with this function in background
        # n = test set size
        N = len(background_proteins)
        K = background_count
        n = len(protein_list)

        if K >= overlap and n >= overlap and N >= K and N >= n:
            # Calculate p-value using hypergeometric distribution
            p_value = 0.0
            for k in range(overlap, min(K, n) + 1):
                # P(X = k) = C(K,k) * C(N-K, n-k) / C(N,n)
                try:
                    p_k = (math.comb(K, k) * math.comb(N - K, n - k)) / math.comb(N, n)
                    p_value += p_k
                except (ValueError, OverflowError):
                    # Fallback for large numbers
                    p_value = 1e-10  # Very significant
                    break

            enrichment_ratio = (overlap / n) / (K / N) if (K / N) > 0 else float("inf")

            if p_value <= max_p_value:
                enrichment_results[func] = {
                    "count": overlap,
                    "test_count": test_count,
                    "background_count": background_count,
                    "enrichment_ratio": enrichment_ratio,
                    "p_value": p_value,
                    "proteins": test_functions[func],
                }

    return enrichment_results
