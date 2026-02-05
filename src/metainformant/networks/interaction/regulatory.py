"""Regulatory network analysis and modeling.

This module provides tools for analyzing and modeling gene regulatory networks,
including network inference, motif discovery, and dynamic analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Optional dependencies
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    logger.warning("networkx not available, regulatory network analysis disabled")

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def construct_regulatory_network(
    regulators: List[str], targets: List[str], interactions: List[Tuple[str, str, str]], **kwargs: Any
) -> Any:
    """Construct a regulatory network from interaction data.

    Args:
        regulators: List of regulator genes/proteins
        targets: List of target genes/proteins
        interactions: List of (regulator, target, interaction_type) tuples
        **kwargs: Additional network parameters

    Returns:
        NetworkX directed graph representing regulatory network

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for regulatory network construction")

    G = nx.DiGraph(**kwargs)

    # Add nodes
    all_nodes = set(regulators) | set(targets)
    G.add_nodes_from(all_nodes)

    # Add regulatory interactions
    for regulator, target, interaction_type in interactions:
        if regulator in all_nodes and target in all_nodes:
            G.add_edge(regulator, target, interaction_type=interaction_type)

    logger.info(f"Constructed regulatory network: {len(G.nodes())} nodes, {len(G.edges())} interactions")
    return G


def infer_regulatory_network_from_expression(
    expression_data: Any, regulators: List[str], method: str = "correlation", threshold: float = 0.5, **kwargs: Any
) -> Any:
    """Infer regulatory network from gene expression data.

    Args:
        expression_data: pandas DataFrame with genes as rows, samples as columns
        regulators: List of potential regulator genes
        method: Inference method ('correlation', 'mutual_info', 'grnboost2')
        threshold: Edge weight threshold
        **kwargs: Additional inference parameters

    Returns:
        NetworkX regulatory network

    Raises:
        ImportError: If required packages not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for regulatory network inference")

    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas required for expression data handling")

    if not isinstance(expression_data, pd.DataFrame):
        raise ValueError("expression_data must be a pandas DataFrame")

    G = nx.DiGraph()

    # Add regulator nodes
    G.add_nodes_from(regulators)

    # Add target nodes (all genes not in regulators)
    targets = [gene for gene in expression_data.index if gene not in regulators]
    G.add_nodes_from(targets)

    if method == "correlation":
        # Pearson correlation between regulators and targets
        for regulator in regulators:
            if regulator not in expression_data.index:
                continue

            reg_expr = expression_data.loc[regulator]

            for target in targets:
                if target not in expression_data.index:
                    continue

                target_expr = expression_data.loc[target]

                # Calculate correlation
                try:
                    corr = reg_expr.corr(target_expr)
                    if abs(corr) >= threshold:
                        G.add_edge(regulator, target, weight=corr, method="correlation")
                except (ValueError, TypeError):
                    continue

    elif method == "mutual_info":
        # Mutual information between regulators and targets
        try:
            from sklearn.feature_selection import mutual_info_regression

            # Transpose for sklearn format (samples x features)
            expr_T = expression_data.T

            for regulator in regulators:
                if regulator not in expression_data.index:
                    continue

                reg_idx = expression_data.index.get_loc(regulator)
                reg_expr = expr_T.iloc[:, reg_idx]

                for target in targets:
                    if target not in expression_data.index:
                        continue

                    target_idx = expression_data.index.get_loc(target)
                    target_expr = expr_T.iloc[:, target_idx]

                    # Calculate mutual information
                    mi = mutual_info_regression(reg_expr.values.reshape(-1, 1), target_expr.values)[0]

                    if mi >= threshold:
                        G.add_edge(regulator, target, weight=mi, method="mutual_info")

        except ImportError:
            raise ImportError("scikit-learn required for mutual information inference")

    elif method == "grnboost2":
        # GRNBoost2 inference (would require arboreto package)
        raise NotImplementedError("GRNBoost2 inference requires arboreto package")

    else:
        raise ValueError(f"Unknown inference method: {method}")

    # Remove isolated nodes
    isolated = list(nx.isolates(G))
    G.remove_nodes_from(isolated)

    logger.info(f"Inferred regulatory network: {len(G.nodes())} nodes, {len(G.edges())} edges")
    return G


def analyze_regulatory_motifs(
    regulatory_graph: Any, motif_types: Optional[List[str]] = None, **kwargs: Any
) -> Dict[str, Any]:
    """Analyze regulatory motifs in the network.

    Args:
        regulatory_graph: NetworkX regulatory graph
        motif_types: Types of motifs to analyze
        **kwargs: Additional analysis parameters

    Returns:
        Dictionary with motif analysis results

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for motif analysis")

    if motif_types is None:
        motif_types = ["feed_forward", "feedback", "cascade"]

    motifs = {}

    # Feed-forward loop analysis
    if "feed_forward" in motif_types:
        motifs["feed_forward_loops"] = _find_feed_forward_loops(regulatory_graph)

    # Feedback loop analysis
    if "feedback" in motif_types:
        motifs["feedback_loops"] = _find_feedback_loops(regulatory_graph)

    # Cascade analysis
    if "cascade" in motif_types:
        motifs["cascades"] = _find_regulatory_cascades(regulatory_graph)

    return {
        "motifs": motifs,
        "network_size": len(regulatory_graph.nodes()),
        "motif_types_analyzed": motif_types,
    }


def _find_feed_forward_loops(G: Any) -> List[List[str]]:
    """Find feed-forward loop motifs in regulatory network."""
    ffl_motifs = []

    # FFL: A → B, A → C, B → C
    for a in G.nodes():
        for b in G.successors(a):
            for c in G.successors(b):
                if G.has_edge(a, c):
                    # Check edge types (activation/repression)
                    ab_type = G[a][b].get("interaction_type", "activation")
                    ac_type = G[a][c].get("interaction_type", "activation")
                    bc_type = G[b][c].get("interaction_type", "activation")

                    ffl_motifs.append(
                        {
                            "nodes": [a, b, c],
                            "edges": [(a, b, ab_type), (a, c, ac_type), (b, c, bc_type)],
                            "type": "feed_forward_loop",
                        }
                    )

    return ffl_motifs


def _find_feedback_loops(G: Any) -> List[List[str]]:
    """Find feedback loop motifs in regulatory network."""
    feedback_motifs = []

    # Simple feedback: A → B → A
    for a in G.nodes():
        for b in G.successors(a):
            if G.has_edge(b, a):
                ab_type = G[a][b].get("interaction_type", "activation")
                ba_type = G[b][a].get("interaction_type", "activation")

                feedback_motifs.append(
                    {"nodes": [a, b], "edges": [(a, b, ab_type), (b, a, ba_type)], "type": "feedback_loop"}
                )

    return feedback_motifs


def _find_regulatory_cascades(G: Any, max_length: int = 5) -> List[List[str]]:
    """Find regulatory cascades in the network."""
    cascades = []

    # Find paths of length 2 to max_length
    for source in G.nodes():
        for target in G.nodes():
            if source != target:
                try:
                    # Find all simple paths
                    paths = list(nx.all_simple_paths(G, source, target, cutoff=max_length))

                    for path in paths:
                        if len(path) >= 3:  # At least 3 nodes
                            cascades.append(
                                {
                                    "nodes": path,
                                    "length": len(path) - 1,  # Number of edges
                                    "type": "regulatory_cascade",
                                }
                            )

                except nx.NetworkXNoPath:
                    continue

    return cascades


def calculate_regulatory_influence(regulatory_graph: Any, source_nodes: List[str], **kwargs: Any) -> Dict[str, Any]:
    """Calculate regulatory influence propagation in the network.

    Args:
        regulatory_graph: NetworkX regulatory graph
        source_nodes: Source nodes for influence calculation
        **kwargs: Additional calculation parameters

    Returns:
        Dictionary with influence analysis results

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for influence calculation")

    # Filter to valid source nodes
    valid_sources = [n for n in source_nodes if n in regulatory_graph.nodes()]

    if not valid_sources:
        return {"error": "No valid source nodes found"}

    # Calculate influence using shortest paths or random walk
    influence_scores = {}

    for node in regulatory_graph.nodes():
        total_influence = 0

        for source in valid_sources:
            if source == node:
                total_influence += 1.0  # Self-influence
            else:
                try:
                    # Use shortest path distance as influence measure
                    distance = nx.shortest_path_length(regulatory_graph, source, node)
                    # Influence decays with distance
                    influence = 1.0 / (distance + 1)
                    total_influence += influence
                except nx.NetworkXNoPath:
                    # No path, no influence
                    pass

        influence_scores[node] = total_influence

    # Normalize influence scores
    max_influence = max(influence_scores.values()) if influence_scores else 1.0
    if max_influence > 0:
        influence_scores = {node: score / max_influence for node, score in influence_scores.items()}

    return {
        "influence_scores": influence_scores,
        "source_nodes": valid_sources,
        "total_nodes": len(regulatory_graph.nodes()),
        "method": "shortest_path_decay",
    }


def analyze_regulatory_dynamics(
    regulatory_graph: Any, initial_states: Dict[str, float], time_steps: int = 10, **kwargs: Any
) -> Dict[str, Any]:
    """Simulate regulatory network dynamics.

    Args:
        regulatory_graph: NetworkX regulatory graph
        initial_states: Initial expression states for nodes
        time_steps: Number of simulation time steps
        **kwargs: Additional simulation parameters

    Returns:
        Dictionary with dynamics simulation results

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for dynamics simulation")

    if not HAS_NUMPY:
        raise ImportError("numpy required for dynamics simulation")

    nodes = list(regulatory_graph.nodes())
    n_nodes = len(nodes)

    # Initialize state vector
    states = np.zeros((time_steps + 1, n_nodes))

    # Set initial states
    for i, node in enumerate(nodes):
        states[0, i] = initial_states.get(node, 0.0)

    # Simple regulatory dynamics simulation
    # Each node's state is influenced by its regulators
    for t in range(1, time_steps + 1):
        for i, node in enumerate(nodes):
            regulators = list(regulatory_graph.predecessors(node))

            if regulators:
                # Simple averaging of regulator states
                regulator_states = [states[t - 1, nodes.index(reg)] for reg in regulators]
                new_state = np.mean(regulator_states)

                # Add some noise and decay
                noise = np.random.normal(0, 0.1)
                decay = 0.9

                states[t, i] = decay * new_state + (1 - decay) * states[t - 1, i] + noise
            else:
                # No regulators, maintain state with decay
                states[t, i] = 0.9 * states[t - 1, i]

    # Convert to dictionary format
    dynamics = {}
    for i, node in enumerate(nodes):
        dynamics[node] = {
            "initial_state": float(states[0, i]),
            "final_state": float(states[-1, i]),
            "trajectory": states[:, i].tolist(),
            "mean_state": float(np.mean(states[:, i])),
            "state_variance": float(np.var(states[:, i])),
        }

    return {
        "dynamics": dynamics,
        "time_steps": time_steps,
        "simulation_method": "simple_regulatory_model",
        "parameters": kwargs,
    }


def identify_regulatory_hubs(regulatory_graph: Any, direction: str = "out", **kwargs: Any) -> List[Tuple[str, int]]:
    """Identify regulatory hubs in the network.

    Args:
        regulatory_graph: NetworkX regulatory graph
        direction: Direction for hub identification ('in', 'out', 'both')
        **kwargs: Additional hub identification parameters

    Returns:
        List of (node, degree) tuples for hub nodes

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for hub identification")

    if direction == "out":
        degrees = dict(regulatory_graph.out_degree())
    elif direction == "in":
        degrees = dict(regulatory_graph.in_degree())
    elif direction == "both":
        degrees = dict(regulatory_graph.degree())
    else:
        raise ValueError(f"Unknown direction: {direction}")

    # Sort by degree
    sorted_hubs = sorted(degrees.items(), key=lambda x: x[1], reverse=True)

    # Return top hubs (e.g., top 10% or absolute threshold)
    threshold = kwargs.get("threshold", 0.1)  # Top 10%
    if isinstance(threshold, float) and threshold < 1:
        n_hubs = max(1, int(len(sorted_hubs) * threshold))
        return sorted_hubs[:n_hubs]
    else:
        # Absolute threshold
        min_degree = threshold if isinstance(threshold, int) else 5
        return [(node, deg) for node, deg in sorted_hubs if deg >= min_degree]


def regulatory_network_stability_analysis(regulatory_graph: Any, **kwargs: Any) -> Dict[str, Any]:
    """Analyze stability properties of regulatory network.

    Args:
        regulatory_graph: NetworkX regulatory graph
        **kwargs: Additional analysis parameters

    Returns:
        Dictionary with stability analysis results

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for stability analysis")

    analysis = {
        "structural_stability": {},
        "dynamical_stability": {},
    }

    # Structural stability measures
    if nx.is_weakly_connected(regulatory_graph):
        analysis["structural_stability"]["connected"] = True

        # Average path length
        try:
            avg_path = nx.average_shortest_path_length(regulatory_graph.to_undirected())
            analysis["structural_stability"]["average_path_length"] = avg_path
        except (nx.NetworkXError, ZeroDivisionError):
            analysis["structural_stability"]["average_path_length"] = None

    else:
        analysis["structural_stability"]["connected"] = False
        components = list(nx.weakly_connected_components(regulatory_graph))
        analysis["structural_stability"]["n_components"] = len(components)

    # Feedback loops (indicate potential oscillations)
    feedback_loops = _find_feedback_loops(regulatory_graph)
    analysis["dynamical_stability"]["n_feedback_loops"] = len(feedback_loops)

    # Network motifs that affect stability
    motifs = analyze_regulatory_motifs(regulatory_graph)
    analysis["dynamical_stability"]["motifs"] = motifs

    return analysis


def export_regulatory_network(regulatory_graph: Any, output_file: str, format: str = "sif", **kwargs: Any) -> None:
    """Export regulatory network in various formats.

    Args:
        regulatory_graph: NetworkX regulatory graph
        output_file: Output file path
        format: Export format ('sif', 'dot', 'json')
        **kwargs: Additional export parameters

    Raises:
        ImportError: If networkx not available
    """
    if not HAS_NETWORKX:
        raise ImportError("networkx required for network export")

    if format == "sif":
        # Simple Interaction Format
        with open(output_file, "w") as f:
            for u, v, data in regulatory_graph.edges(data=True):
                interaction_type = data.get("interaction_type", "regulates")
                f.write(f"{u}\t{interaction_type}\t{v}\n")

    elif format == "dot":
        # GraphViz DOT format
        from networkx.drawing.nx_pydot import write_dot

        write_dot(regulatory_graph, output_file)

    elif format == "json":
        # Node-link JSON format
        import json

        data = nx.node_link_data(regulatory_graph)
        with open(output_file, "w") as f:
            json.dump(data, f, indent=2)

    else:
        raise ValueError(f"Unsupported export format: {format}")

    logger.info(f"Exported regulatory network to {output_file} ({format} format)")


class GeneRegulatoryNetwork:
    """A gene regulatory network class for analyzing transcriptional regulation.

    This class provides methods for loading, analyzing, and modeling gene regulatory
    networks, including transcription factor-target relationships and regulatory motifs.
    """

    def __init__(self, graph: Optional[Any] = None, name: str = ""):
        """Initialize a gene regulatory network.

        Args:
            graph: NetworkX DiGraph object (optional)
            name: Name of the network
        """
        if graph is None and HAS_NETWORKX:
            self.graph = nx.DiGraph()
        else:
            self.graph = graph
        self.name = name
        self.metadata = {}
        self.tf_targets = {}  # Map of transcription factors to their targets
        self.target_tfs = {}  # Map of targets to their regulators

    @classmethod
    def from_tf_target_file(cls, filepath: Union[str, Path], name: str = "") -> "GeneRegulatoryNetwork":
        """Load regulatory network from TF-target file.

        Args:
            filepath: Path to TF-target file (format: TF<TAB>Target<TAB>Confidence)
            name: Name for the network

        Returns:
            GeneRegulatoryNetwork instance
        """
        if not HAS_NETWORKX:
            raise ImportError("NetworkX required for GeneRegulatoryNetwork")

        network = cls(name=name or Path(filepath).stem)

        with open(filepath, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    tf = parts[0]
                    target = parts[1]
                    confidence = float(parts[2]) if len(parts) > 2 else 1.0

                    # Add edge from TF to target
                    network.graph.add_edge(tf, target, weight=confidence, type="regulates")

                    # Update mappings
                    if tf not in network.tf_targets:
                        network.tf_targets[tf] = []
                    network.tf_targets[tf].append(target)

                    if target not in network.target_tfs:
                        network.target_tfs[target] = []
                    network.target_tfs[target].append(tf)

        network.metadata["source_file"] = str(filepath)
        network.metadata["n_regulations"] = len(network.graph.edges())
        network.metadata["n_tfs"] = len(network.tf_targets)
        network.metadata["n_targets"] = len(network.target_tfs)

        return network

    @classmethod
    def from_interactions(cls, interactions: List[Tuple[str, str, float]], name: str = "") -> "GeneRegulatoryNetwork":
        """Create network from list of TF-target interactions.

        Args:
            interactions: List of (TF, target, confidence) tuples
            name: Name for the network

        Returns:
            GeneRegulatoryNetwork instance
        """
        if not HAS_NETWORKX:
            raise ImportError("NetworkX required for GeneRegulatoryNetwork")

        network = cls(name=name)

        for tf, target, confidence in interactions:
            network.graph.add_edge(tf, target, weight=confidence, type="regulates")

            # Update mappings
            if tf not in network.tf_targets:
                network.tf_targets[tf] = []
            network.tf_targets[tf].append(target)

            if target not in network.target_tfs:
                network.target_tfs[target] = []
            network.target_tfs[target].append(tf)

        network.metadata["n_regulations"] = len(interactions)
        network.metadata["n_tfs"] = len(network.tf_targets)
        network.metadata["n_targets"] = len(network.target_tfs)

        return network

    def get_transcription_factors(self) -> List[str]:
        """Get list of all transcription factors in the network.

        Returns:
            List of TF identifiers
        """
        return list(self.tf_targets.keys())

    def get_targets(self) -> List[str]:
        """Get list of all target genes in the network.

        Returns:
            List of target gene identifiers
        """
        return list(self.target_tfs.keys())

    def get_tf_targets(self, tf: str) -> List[str]:
        """Get targets regulated by a specific transcription factor.

        Args:
            tf: Transcription factor identifier

        Returns:
            List of target genes
        """
        return self.tf_targets.get(tf, [])

    def get_target_regulators(self, target: str) -> List[str]:
        """Get transcription factors that regulate a specific target.

        Args:
            target: Target gene identifier

        Returns:
            List of regulating TFs
        """
        return self.target_tfs.get(target, [])

    def find_common_targets(self, tf1: str, tf2: str) -> List[str]:
        """Find genes regulated by both transcription factors.

        Args:
            tf1: First transcription factor
            tf2: Second transcription factor

        Returns:
            List of genes regulated by both TFs
        """
        targets1 = set(self.get_tf_targets(tf1))
        targets2 = set(self.get_tf_targets(tf2))
        return list(targets1 & targets2)

    def find_shared_regulators(self, gene1: str, gene2: str) -> List[str]:
        """Find transcription factors that regulate both genes.

        Args:
            gene1: First gene
            gene2: Second gene

        Returns:
            List of TFs regulating both genes
        """
        regulators1 = set(self.get_target_regulators(gene1))
        regulators2 = set(self.get_target_regulators(gene2))
        return list(regulators1 & regulators2)

    def regulatory_motifs(self) -> Dict[str, List[List[str]]]:
        """Identify regulatory motifs in the network.

        Returns:
            Dictionary of motif types and their instances
        """
        motifs = {
            "feed_forward": [],  # TF1 -> TF2 -> Target
            "auto_regulation": [],  # TF regulates itself
            "mutual_regulation": [],  # TF1 <-> TF2
            "cascade": [],  # TF1 -> TF2 -> TF3
        }

        # Find auto-regulation
        for tf in self.get_transcription_factors():
            if tf in self.get_tf_targets(tf):
                motifs["auto_regulation"].append([tf])

        # Find feed-forward loops
        for tf1 in self.get_transcription_factors():
            for tf2 in self.get_tf_targets(tf1):
                if tf2 in self.tf_targets:  # TF2 is also a TF
                    for target in self.get_tf_targets(tf2):
                        if target not in [tf1, tf2]:  # Avoid self-loops
                            motifs["feed_forward"].append([tf1, tf2, target])

        # Find mutual regulation
        for tf1 in self.get_transcription_factors():
            for tf2 in self.get_tf_targets(tf1):
                if tf2 in self.tf_targets and tf1 in self.get_tf_targets(tf2):
                    if [tf1, tf2] not in motifs["mutual_regulation"]:
                        motifs["mutual_regulation"].append([tf1, tf2])

        return motifs

    def network_summary(self) -> Dict[str, Any]:
        """Generate network summary statistics.

        Returns:
            Dictionary with network statistics
        """
        if not HAS_NETWORKX:
            return {"error": "NetworkX not available"}

        summary = {
            "name": self.name,
            "n_nodes": len(self.graph.nodes()),
            "n_edges": len(self.graph.edges()),
            "n_tfs": len(self.tf_targets),
            "n_targets": len(self.target_tfs),
            "density": nx.density(self.graph),
            "average_degree": (
                sum(dict(self.graph.degree()).values()) / len(self.graph.nodes()) if self.graph.nodes() else 0
            ),
            "is_directed": self.graph.is_directed(),
        }

        # Add motif information
        motifs = self.regulatory_motifs()
        summary["motifs"] = {motif_type: len(instances) for motif_type, instances in motifs.items()}

        # Add metadata
        summary.update(self.metadata)

        return summary

    def __len__(self) -> int:
        """Get number of regulations."""
        if not HAS_NETWORKX:
            return 0
        return len(self.graph.edges())

    def __contains__(self, gene: str) -> bool:
        """Check if gene is in network."""
        if not HAS_NETWORKX:
            return False
        return gene in self.graph

    def add_regulation(self, tf: str, target: str, **kwargs: Any) -> None:
        """Add a regulatory interaction from TF to target.

        Args:
            tf: Transcription factor identifier
            target: Target gene identifier
            **kwargs: Additional edge attributes (e.g., weight, confidence)
        """
        if not HAS_NETWORKX:
            raise ImportError("networkx required for adding regulation")

        # Add edge
        self.graph.add_edge(tf, target, type="regulates", **kwargs)

        # Update mappings
        if tf not in self.tf_targets:
            self.tf_targets[tf] = []
        if target not in self.tf_targets[tf]:
            self.tf_targets[tf].append(target)

        if target not in self.target_tfs:
            self.target_tfs[target] = []
        if tf not in self.target_tfs[target]:
            self.target_tfs[target].append(tf)

    def add_gene_metadata(self, gene: str, **metadata: Any) -> None:
        """Add metadata to a gene node.

        Args:
            gene: Gene identifier
            **metadata: Metadata key-value pairs to add
        """
        if not HAS_NETWORKX:
            raise ImportError("networkx required for adding metadata")
        if gene not in self.graph:
            self.graph.add_node(gene)
        self.graph.nodes[gene].update(metadata)


def infer_grn(
    expression_data: Any,
    gene_names: List[str],
    method: str = "correlation",
    threshold: float = 0.5,
    tf_genes: Optional[List[str]] = None,
    **kwargs: Any,
) -> GeneRegulatoryNetwork:
    """Infer gene regulatory network from expression data.

    Args:
        expression_data: Expression matrix (samples x genes)
        gene_names: List of gene names
        method: Inference method ('correlation', 'mutual_info', 'granger')
        threshold: Correlation threshold for edge inclusion
        tf_genes: Optional list of known transcription factors
        **kwargs: Additional method-specific parameters

    Returns:
        GeneRegulatoryNetwork instance
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for GRN inference")

    # Convert expression data to numpy array if needed
    if hasattr(expression_data, "values"):
        expression_matrix = expression_data.values
    else:
        expression_matrix = expression_data

    n_samples, n_genes = expression_matrix.shape
    if len(gene_names) != n_genes:
        raise ValueError("gene_names length must match expression_data columns")

    # Create network
    network = GeneRegulatoryNetwork(name=f"GRN_{method}")

    # Infer regulatory relationships based on method
    if method == "correlation":
        # Use correlation to infer regulatory relationships
        import numpy as np

        # Calculate correlation matrix
        corr_matrix = np.corrcoef(expression_matrix.T)

        # Identify potential TFs (if not provided)
        if tf_genes is None:
            # Simple heuristic: genes with high variance are likely TFs
            gene_variance = np.var(expression_matrix, axis=0)
            variance_threshold = np.percentile(gene_variance, 75)  # Top 25%
            tf_indices = np.where(gene_variance >= variance_threshold)[0]
            tf_genes = [gene_names[i] for i in tf_indices]

        # Add regulatory edges
        for i, tf in enumerate(tf_genes):
            if tf not in gene_names:
                continue
            tf_idx = gene_names.index(tf)

            for j, target in enumerate(gene_names):
                if i == j:  # Skip self-regulation
                    continue

                correlation = abs(corr_matrix[tf_idx, j])
                if correlation >= threshold:
                    # Add edge from TF to target
                    network.graph.add_edge(tf, target, weight=correlation, type="regulates")

                    # Update mappings
                    if tf not in network.tf_targets:
                        network.tf_targets[tf] = []
                    network.tf_targets[tf].append(target)

                    if target not in network.target_tfs:
                        network.target_tfs[target] = []
                    network.target_tfs[target].append(tf)

    elif method == "mutual_info":
        # Use mutual information for inference
        from sklearn.feature_selection import mutual_info_regression
        import numpy as np

        if tf_genes is None:
            # Use all genes as potential TFs
            tf_genes = gene_names

        for tf in tf_genes:
            tf_idx = gene_names.index(tf)
            tf_expression = expression_matrix[:, tf_idx]

            for j, target in enumerate(gene_names):
                if gene_names[j] == tf:
                    continue

                target_expression = expression_matrix[:, j]

                # Calculate mutual information
                mi = mutual_info_regression(tf_expression.reshape(-1, 1), target_expression)[0]

                if mi >= threshold:
                    network.graph.add_edge(tf, target, weight=mi, type="regulates")

                    # Update mappings
                    if tf not in network.tf_targets:
                        network.tf_targets[tf] = []
                    network.tf_targets[tf].append(target)

                    if target not in network.target_tfs:
                        network.target_tfs[target] = []
                    network.target_tfs[target].append(tf)

    else:
        raise ValueError(f"Unknown GRN inference method: {method}")

    # Set network attributes
    network.metadata.update(
        {
            "inference_method": method,
            "threshold": threshold,
            "n_samples": n_samples,
            "n_genes": n_genes,
            "n_tfs": len(network.tf_targets),
            "n_targets": len(network.target_tfs),
            "n_regulations": len(network.graph.edges()),
        }
    )

    return network


def regulatory_motifs(grn: GeneRegulatoryNetwork) -> List[Dict[str, Any]]:
    """Identify regulatory motifs in a gene regulatory network.

    This is a standalone function that analyzes regulatory motifs.

    Args:
        grn: GeneRegulatoryNetwork instance

    Returns:
        List of identified motifs with type, genes, and confidence
    """
    motifs = []

    if not HAS_NETWORKX:
        return motifs

    # Find feed-forward loops (TF1 -> TF2 -> Target)
    for tf1 in grn.get_transcription_factors():
        for tf2 in grn.get_tf_targets(tf1):
            if tf2 in grn.tf_targets:  # TF2 is also a TF
                for target in grn.get_tf_targets(tf2):
                    if target not in [tf1, tf2]:  # Avoid self-loops
                        # Check if TF1 regulates TF2 and TF2 regulates target
                        if grn.graph.has_edge(tf1, tf2) and grn.graph.has_edge(tf2, target):
                            motifs.append(
                                {"motif_type": "feed_forward_loop", "genes": [tf1, tf2, target], "confidence": 0.8}
                            )

    # Find auto-regulation (TF regulates itself)
    for tf in grn.get_transcription_factors():
        if grn.graph.has_edge(tf, tf):
            motifs.append({"motif_type": "auto_regulation", "genes": [tf], "confidence": 1.0})

    # Find mutual regulation (TF1 <-> TF2)
    for tf1 in grn.get_transcription_factors():
        for tf2 in grn.get_tf_targets(tf1):
            if tf2 in grn.tf_targets and grn.graph.has_edge(tf2, tf1):
                if [tf1, tf2] not in [
                    [m["genes"][0], m["genes"][1]] for m in motifs if m["motif_type"] == "mutual_regulation"
                ]:
                    motifs.append({"motif_type": "mutual_regulation", "genes": [tf1, tf2], "confidence": 0.9})

    return motifs


def pathway_regulation_analysis(grn: GeneRegulatoryNetwork, pathway_genes: List[str]) -> Dict[str, Any]:
    """Analyze regulation patterns within a pathway.

    Args:
        grn: GeneRegulatoryNetwork instance
        pathway_genes: List of genes in the pathway

    Returns:
        Dictionary with pathway regulation analysis results
    """
    results = {
        "pathway_tfs": [],
        "internal_regulations": 0,
        "external_regulations": 0,
        "regulation_density": 0.0,
        "hub_tfs": [],
        "isolated_genes": [],
    }

    if not HAS_NETWORKX:
        return results

    pathway_set = set(pathway_genes)

    # Find TFs that regulate pathway genes
    pathway_tfs = set()
    for gene in pathway_genes:
        regulators = grn.get_target_regulators(gene)
        pathway_tfs.update(regulators)

    results["pathway_tfs"] = list(pathway_tfs)

    # Count internal vs external regulations
    internal_regulations = 0
    external_regulations = 0

    for tf in pathway_tfs:
        targets = grn.get_tf_targets(tf)
        for target in targets:
            if target in pathway_set:
                if tf in pathway_set:
                    internal_regulations += 1
                else:
                    external_regulations += 1

    results["internal_regulations"] = internal_regulations
    results["external_regulations"] = external_regulations

    # Calculate regulation density
    if len(pathway_set) > 0:
        total_possible_regulations = len(pathway_tfs) * len(pathway_set)
        if total_possible_regulations > 0:
            results["regulation_density"] = (internal_regulations + external_regulations) / total_possible_regulations

    # Identify hub TFs (regulate many pathway genes)
    tf_target_counts = {}
    for tf in pathway_tfs:
        targets = grn.get_tf_targets(tf)
        pathway_targets = [t for t in targets if t in pathway_set]
        tf_target_counts[tf] = len(pathway_targets)

    # TFs regulating 3+ pathway genes are considered hubs
    results["hub_tfs"] = [tf for tf, count in tf_target_counts.items() if count >= 3]

    # Find isolated genes (no regulation within pathway)
    isolated_genes = []
    for gene in pathway_genes:
        regulators = grn.get_target_regulators(gene)
        pathway_regulators = [r for r in regulators if r in pathway_set]
        if not pathway_regulators:
            isolated_genes.append(gene)

    results["isolated_genes"] = isolated_genes

    return results
