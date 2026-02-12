"""Regulatory network core: GeneRegulatoryNetwork class and network operations.

This module provides the GeneRegulatoryNetwork class and supporting functions
for constructing, analyzing motifs, calculating influence, simulating dynamics,
identifying hubs, assessing stability, and exporting regulatory networks.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging

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

    # FFL: A -> B, A -> C, B -> C
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

    # Simple feedback: A -> B -> A
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

    def __init__(self, name_or_graph: Any = "", name: str = "", **kwargs: Any):
        """Initialize a gene regulatory network.

        Args:
            name_or_graph: Name string or NetworkX DiGraph for backward compatibility
            name: Name of the network (used when name_or_graph is a graph)
        """
        if isinstance(name_or_graph, str):
            self.name = name_or_graph or name
            if HAS_NETWORKX:
                self.graph = nx.DiGraph()
            else:
                self.graph = None
        else:
            # Backward compat: first arg is a graph
            self.graph = name_or_graph
            self.name = name
            if self.graph is None and HAS_NETWORKX:
                self.graph = nx.DiGraph()
        self.metadata = {}
        self.tf_targets = {}  # Map of transcription factors to their targets
        self.target_tfs = {}  # Map of targets to their regulators
        self._gene_metadata: Dict[str, Dict[str, Any]] = {}

    @property
    def regulations(self) -> List[Tuple[str, str, Dict[str, Any]]]:
        """Get all regulatory interactions as tuples."""
        result = []
        if HAS_NETWORKX and self.graph is not None:
            for tf, target, data in self.graph.edges(data=True):
                meta = dict(data)
                # Map internal keys to expected keys
                if "regulation_type" in meta:
                    meta["type"] = meta.pop("regulation_type")
                if "weight" in meta:
                    meta["strength"] = meta.get("weight", 0.0)
                result.append((tf, target, meta))
        return result

    @property
    def genes(self) -> List[str]:
        """Get all genes in the network."""
        if HAS_NETWORKX and self.graph is not None:
            return list(self.graph.nodes())
        all_genes: set[str] = set()
        for tf, targets in self.tf_targets.items():
            all_genes.add(tf)
            all_genes.update(targets)
        return list(all_genes)

    @property
    def transcription_factors(self) -> List[str]:
        """Get all transcription factors."""
        return list(self.tf_targets.keys())

    @property
    def gene_metadata(self) -> Dict[str, Dict[str, Any]]:
        """Get gene metadata dictionary."""
        if HAS_NETWORKX and self.graph is not None:
            result: Dict[str, Dict[str, Any]] = {}
            for node, data in self.graph.nodes(data=True):
                if data:
                    result[node] = dict(data)
            result.update(self._gene_metadata)
            return result
        return self._gene_metadata

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

    def get_targets(self, tf: Optional[str] = None) -> List[str]:
        """Get target genes, optionally for a specific TF.

        Args:
            tf: If provided, get targets of this specific TF. Otherwise get all targets.

        Returns:
            List of target gene identifiers
        """
        if tf is not None:
            return self.get_tf_targets(tf)
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

        # Add edge (only set default type if not provided in kwargs)
        if "type" not in kwargs:
            kwargs["type"] = "regulates"
        self.graph.add_edge(tf, target, **kwargs)

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
        if gene not in self._gene_metadata:
            self._gene_metadata[gene] = {}
        self._gene_metadata[gene].update(metadata)

    def add_transcription_factor(self, gene_id: str, **kwargs: Any) -> None:
        """Register a gene as a transcription factor.

        Args:
            gene_id: Gene identifier
            **kwargs: Additional metadata (tf_family, etc.)
        """
        if gene_id not in self.tf_targets:
            self.tf_targets[gene_id] = []
        if HAS_NETWORKX and self.graph is not None:
            if gene_id not in self.graph:
                self.graph.add_node(gene_id)
            self.graph.nodes[gene_id].update(kwargs)
            self.graph.nodes[gene_id]["is_tf"] = True
        if gene_id not in self._gene_metadata:
            self._gene_metadata[gene_id] = {}
        self._gene_metadata[gene_id].update(kwargs)
        self._gene_metadata[gene_id]["is_tf"] = True

    def get_regulators(self, target: str) -> List[str]:
        """Get regulators of a target gene. Alias for get_target_regulators."""
        return self.get_target_regulators(target)

    def filter_by_confidence(self, threshold: float) -> "GeneRegulatoryNetwork":
        """Filter regulations by confidence/weight threshold.

        Args:
            threshold: Minimum confidence value

        Returns:
            New GeneRegulatoryNetwork with filtered regulations
        """
        filtered = GeneRegulatoryNetwork(name=f"{self.name}_filtered")
        if HAS_NETWORKX and self.graph is not None:
            for tf, target, data in self.graph.edges(data=True):
                confidence = data.get("confidence", data.get("weight", 1.0))
                if confidence >= threshold:
                    filtered.add_regulation(tf, target, **data)
        return filtered

    def filter_by_regulation_type(self, reg_type: str) -> "GeneRegulatoryNetwork":
        """Filter regulations by type.

        Args:
            reg_type: Regulation type to keep (e.g., 'activation', 'repression')

        Returns:
            New GeneRegulatoryNetwork with filtered regulations
        """
        filtered = GeneRegulatoryNetwork(name=f"{self.name}_{reg_type}")
        if HAS_NETWORKX and self.graph is not None:
            for tf, target, data in self.graph.edges(data=True):
                edge_type = data.get("regulation_type", data.get("type", ""))
                if edge_type == reg_type:
                    filtered.add_regulation(tf, target, **data)
        return filtered

    def get_network_statistics(self) -> Dict[str, Any]:
        """Get comprehensive network statistics.

        Returns:
            Dictionary with network statistics
        """
        stats: Dict[str, Any] = {
            "n_regulations": len(self),
            "n_genes": len(self.genes),
            "n_transcription_factors": len(self.transcription_factors),
            "n_targets": len(self.target_tfs),
        }
        if HAS_NETWORKX and self.graph is not None:
            stats["density"] = nx.density(self.graph)
            stats["n_nodes"] = len(self.graph.nodes())
            stats["n_edges"] = len(self.graph.edges())
        return stats

    def to_biological_network(self) -> Any:
        """Convert to a BiologicalNetwork instance.

        Returns:
            BiologicalNetwork instance
        """
        from metainformant.networks.analysis.graph import BiologicalNetwork

        bio_net = BiologicalNetwork()
        if HAS_NETWORKX and self.graph is not None:
            for node in self.graph.nodes():
                bio_net.add_node(node)
            for tf, target, data in self.graph.edges(data=True):
                bio_net.add_edge(tf, target, **data)
        return bio_net
