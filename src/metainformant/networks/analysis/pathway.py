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
    from metainformant.networks.analysis.community import detect_communities

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


class PathwayNetwork:
    """A biological pathway network class for pathway analysis and visualization.

    This class provides methods for loading, analyzing, and visualizing biological
    pathways as network structures.
    """

    def __init__(self, name: str = "", pathways: Optional[Dict[str, List[str]]] = None):
        """Initialize a pathway network.

        Args:
            name: Name of the pathway network
            pathways: Dictionary mapping pathway names to gene lists
        """
        self.name = name
        self.pathways = pathways or {}
        self.metadata = {}

    @classmethod
    def load_from_database(cls, filepath: str | Path, format: str = "json") -> 'PathwayNetwork':
        """Load pathway network from a database file.

        Args:
            filepath: Path to the pathway database file
            format: File format ('json', 'gmt')

        Returns:
            PathwayNetwork instance
        """
        import json
        from pathlib import Path

        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"Pathway database file not found: {filepath}")

        if format.lower() == "json":
            with open(filepath, 'r') as f:
                data = json.load(f)
            pathways = data.get('pathways', {})
            name = data.get('name', filepath.stem)
        elif format.lower() == "gmt":
            # GMT format: pathway_name<TAB>description<TAB>gene1<TAB>gene2...
            pathways = {}
            with open(filepath, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        pathway_name = parts[0]
                        pathways[pathway_name] = parts[2:]  # Skip description
            name = filepath.stem
        else:
            raise ValueError(f"Unsupported format: {format}")

        network = cls(name=name, pathways=pathways)
        network.metadata['source_file'] = str(filepath)
        network.metadata['format'] = format
        return network

    def get_pathway(self, pathway_name: str) -> List[str]:
        """Get genes in a specific pathway.

        Args:
            pathway_name: Name of the pathway

        Returns:
            List of genes in the pathway
        """
        return self.pathways.get(pathway_name, [])

    def find_pathways_containing_gene(self, gene: str) -> List[str]:
        """Find all pathways containing a specific gene.

        Args:
            gene: Gene name

        Returns:
            List of pathway names containing the gene
        """
        return [name for name, genes in self.pathways.items() if gene in genes]

    def get_pathway_sizes(self) -> Dict[str, int]:
        """Get the size (number of genes) for each pathway.

        Returns:
            Dictionary mapping pathway names to sizes
        """
        return {name: len(genes) for name, genes in self.pathways.items()}

    def filter_pathways_by_size(self, min_size: int = 5, max_size: int | None = None) -> 'PathwayNetwork':
        """Filter pathways by size.

        Args:
            min_size: Minimum pathway size
            max_size: Maximum pathway size (None for no limit)

        Returns:
            New PathwayNetwork with filtered pathways
        """
        filtered_pathways = {}
        for name, genes in self.pathways.items():
            size = len(genes)
            if size >= min_size and (max_size is None or size <= max_size):
                filtered_pathways[name] = genes

        new_network = PathwayNetwork(name=f"{self.name}_filtered", pathways=filtered_pathways)
        new_network.metadata = self.metadata.copy()
        return new_network

    def pathway_overlap_matrix(self) -> Dict[str, Dict[str, float]]:
        """Calculate overlap matrix between pathways.

        Returns:
            Dictionary of dictionaries with Jaccard similarity between pathways
        """
        pathway_names = list(self.pathways.keys())
        overlap_matrix = {}

        for i, name1 in enumerate(pathway_names):
            overlap_matrix[name1] = {}
            genes1 = set(self.pathways[name1])

            for j, name2 in enumerate(pathway_names):
                if i == j:
                    overlap_matrix[name1][name2] = 1.0
                elif j > i:
                    genes2 = set(self.pathways[name2])
                    intersection = len(genes1 & genes2)
                    union = len(genes1 | genes2)
                    jaccard = intersection / union if union > 0 else 0.0
                    overlap_matrix[name1][name2] = jaccard
                    overlap_matrix[name2][name1] = jaccard
                else:
                    overlap_matrix[name1][name2] = overlap_matrix[name2][name1]

        return overlap_matrix

    def __len__(self) -> int:
        """Get number of pathways."""
        return len(self.pathways)

    def __contains__(self, pathway_name: str) -> bool:
        """Check if pathway exists."""
        return pathway_name in self.pathways

    def __getitem__(self, pathway_name: str) -> List[str]:
        """Get pathway genes by name."""
        return self.get_pathway(pathway_name)


def pathway_enrichment(
    gene_list: List[str],
    pathway_network: PathwayNetwork,
    background_genes: Optional[List[str]] = None,
    method: str = "fisher",
    correction: str = "bonferroni",
    min_overlap: int = 1
) -> Dict[str, Dict[str, Any]]:
    """Perform pathway enrichment analysis for a gene list against a pathway network.

    Args:
        gene_list: List of genes to test for enrichment
        pathway_network: PathwayNetwork object containing pathway definitions
        background_genes: Background gene set (if None, uses all genes in pathways)
        method: Statistical test method ('fisher', 'hypergeometric')
        correction: Multiple testing correction ('bonferroni', 'fdr', None)
        min_overlap: Minimum overlap size to report

    Returns:
        Dictionary mapping pathway IDs to enrichment results
    """
    if not pathway_network.pathways:
        logger.warning("Empty pathway network")
        return {}

    # Prepare background genes
    if background_genes is None:
        background_genes = set()
        for genes in pathway_network.pathways.values():
            background_genes.update(genes)
        background_genes = list(background_genes)

    # Convert gene list to set for faster lookup
    query_genes = set(gene_list)
    background_set = set(background_genes)

    # Ensure all query genes are in background
    query_genes = query_genes.intersection(background_set)
    if not query_genes:
        logger.warning("No query genes found in background set")
        return {}

    results = {}

    for pathway_id, pathway_genes in pathway_network.pathways.items():
        pathway_set = set(pathway_genes)
        pathway_set = pathway_set.intersection(background_set)  # Only genes in background

        if len(pathway_set) == 0:
            continue

        # Calculate overlap
        overlap = query_genes.intersection(pathway_set)
        overlap_size = len(overlap)

        if overlap_size < min_overlap:
            continue

        # Calculate statistics
        query_size = len(query_genes)
        pathway_size = len(pathway_set)
        background_size = len(background_set)

        # Hypergeometric/Fisher's exact test
        if HAS_SCIPY:
            if method == "fisher":
                # Fisher's exact test
                from scipy.stats import fisher_exact

                # Create contingency table
                # [[in_query_and_pathway, in_query_not_pathway],
                #  [not_in_query_and_pathway, not_in_query_not_pathway]]
                contingency_table = [
                    [overlap_size, query_size - overlap_size],
                    [pathway_size - overlap_size, background_size - pathway_size - (query_size - overlap_size)]
                ]

                odds_ratio, p_value = fisher_exact(contingency_table, alternative='greater')
            else:
                # Hypergeometric test
                from scipy.stats import hypergeom
                p_value = hypergeom.sf(overlap_size - 1, background_size, pathway_size, query_size)
        else:
            # Simple approximation when scipy not available
            logger.warning("SciPy not available, using simplified enrichment calculation")
            expected = (query_size * pathway_size) / background_size
            if overlap_size > expected:
                p_value = max(0.001, 1.0 / (overlap_size / expected))  # Rough approximation
            else:
                p_value = 1.0

        # Calculate enrichment ratio
        expected_overlap = (query_size * pathway_size) / background_size
        enrichment_ratio = overlap_size / expected_overlap if expected_overlap > 0 else float('inf')

        results[pathway_id] = {
            "overlap_size": overlap_size,
            "pathway_size": pathway_size,
            "query_size": query_size,
            "background_size": background_size,
            "p_value": p_value,
            "enrichment_ratio": enrichment_ratio,
            "overlap_genes": list(overlap),
            "expected_overlap": expected_overlap
        }

    # Apply multiple testing correction
    if correction and len(results) > 1:
        p_values = [result["p_value"] for result in results.values()]

        if correction == "bonferroni":
            corrected_p_values = [min(p * len(results), 1.0) for p in p_values]
        elif correction == "fdr" and HAS_SCIPY:
            # Benjamini-Hochberg FDR correction
            from scipy.stats import rankdata
            sorted_indices = sorted(range(len(p_values)), key=lambda i: p_values[i])
            corrected_p_values = [1.0] * len(p_values)

            for i, idx in enumerate(sorted_indices):
                rank = i + 1
                corrected_p = min(p_values[idx] * len(p_values) / rank, 1.0)
                # Ensure monotonicity
                if i > 0:
                    corrected_p = min(corrected_p, corrected_p_values[sorted_indices[i-1]])
                corrected_p_values[idx] = corrected_p
        else:
            corrected_p_values = p_values

        # Add corrected p-values to results
        for i, pathway_id in enumerate(results.keys()):
            results[pathway_id]["corrected_p_value"] = corrected_p_values[i]

    return results


def load_pathway_database(pathway_data: Dict[str, Any], name: str = "loaded_pathways") -> PathwayNetwork:
    """Load pathway database from structured data.

    Args:
        pathway_data: Dictionary containing pathway information
        name: Name for the pathway network

    Returns:
        PathwayNetwork instance
    """
    pathways = {}

    # Handle different data formats
    if "pathways" in pathway_data:
        # Format: {"pathways": {"pathway_id": {"name": "...", "genes": [...], ...}}}
        pathway_dict = pathway_data["pathways"]
        for pathway_id, pathway_info in pathway_dict.items():
            if isinstance(pathway_info, dict):
                genes = pathway_info.get("genes", [])
                if isinstance(genes, list):
                    pathways[pathway_id] = genes
                elif isinstance(genes, str):
                    pathways[pathway_id] = [genes]
                else:
                    pathways[pathway_id] = list(genes)
            elif isinstance(pathway_info, list):
                pathways[pathway_id] = pathway_info
    else:
        # Assume direct format: {"pathway_id": {"name": "...", "genes": [...], ...}}
        for pathway_id, pathway_info in pathway_data.items():
            if isinstance(pathway_info, dict):
                genes = pathway_info.get("genes", [])
                if isinstance(genes, list):
                    pathways[pathway_id] = genes
                elif isinstance(genes, str):
                    pathways[pathway_id] = [genes]
                else:
                    pathways[pathway_id] = list(genes)
            elif isinstance(pathway_info, list):
                pathways[pathway_id] = pathway_info

    network = PathwayNetwork(name=name, pathways=pathways)

    # Add metadata if available
    if "metadata" in pathway_data:
        network.metadata.update(pathway_data["metadata"])
    if "version" in pathway_data:
        network.metadata["version"] = pathway_data["version"]
    if "source" in pathway_data:
        network.metadata["source"] = pathway_data["source"]

    return network


def network_enrichment_analysis(
    gene_list: List[str],
    pathway_network: PathwayNetwork,
    background_genes: Optional[List[str]] = None,
    method: str = "fisher",
    correction: str = "bonferroni",
    min_overlap: int = 1
) -> Dict[str, Dict[str, Any]]:
    """Perform network-based enrichment analysis.

    This is an alias for pathway_enrichment with network-aware parameters.

    Args:
        gene_list: List of genes to test for enrichment
        pathway_network: PathwayNetwork object containing pathway definitions
        background_genes: Background gene set
        method: Statistical test method
        correction: Multiple testing correction
        min_overlap: Minimum overlap size

    Returns:
        Dictionary mapping pathway IDs to enrichment results
    """
    return pathway_enrichment(
        gene_list=gene_list,
        pathway_network=pathway_network,
        background_genes=background_genes,
        method=method,
        correction=correction,
        min_overlap=min_overlap
    )
