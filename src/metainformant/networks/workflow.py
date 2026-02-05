"""Network analysis workflow orchestration.

Provides a high-level workflow interface for running multi-step network
analyses combining graph construction, community detection, pathway
enrichment, and PPI analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from metainformant.core import logging

from .config import NetworkWorkflowConfig

logger = logging.get_logger(__name__)


class NetworkWorkflow:
    """Orchestrates multi-step network analysis workflows.

    This class provides a high-level interface for running network
    analysis pipelines that combine multiple analysis steps.
    """

    def __init__(self, config: Optional[NetworkWorkflowConfig] = None):
        """Initialize network workflow.

        Args:
            config: Workflow configuration. Uses defaults if None.
        """
        self.config = config or NetworkWorkflowConfig()
        self.results: Dict[str, Any] = {}
        self._steps_completed: List[str] = []

    def build_network(
        self,
        edges: Optional[List[tuple]] = None,
        correlation_matrix: Optional[Any] = None,
        node_names: Optional[List[str]] = None,
        threshold: Optional[float] = None,
    ) -> "NetworkWorkflow":
        """Build a network from input data.

        Args:
            edges: Edge list for direct construction
            correlation_matrix: Correlation matrix for threshold-based construction
            node_names: Node labels (required with correlation_matrix)
            threshold: Edge weight threshold (overrides config)

        Returns:
            self for method chaining
        """
        from .analysis.graph import BiologicalNetwork, create_network, add_edges_from_correlation

        if edges is not None:
            graph = create_network(edges, directed=self.config.network.directed)
            network = BiologicalNetwork(directed=self.config.network.directed)
            network.graph = graph
        elif correlation_matrix is not None:
            import networkx as nx

            if self.config.network.directed:
                g = nx.DiGraph()
            else:
                g = nx.Graph()
            if node_names:
                g.add_nodes_from(node_names)
            network = BiologicalNetwork(directed=self.config.network.directed)
            network.graph = g
            t = threshold if threshold is not None else self.config.network.min_edge_weight
            add_edges_from_correlation(network.graph, correlation_matrix, threshold=t)
        else:
            network = BiologicalNetwork(directed=self.config.network.directed)

        if self.config.network.remove_isolates:
            from .analysis.graph import remove_isolated_nodes

            remove_isolated_nodes(network.graph)

        self.results["network"] = network
        self._steps_completed.append("build_network")
        logger.info(f"Built network: {network.number_of_nodes()} nodes, " f"{network.number_of_edges()} edges")
        return self

    def detect_communities(self) -> "NetworkWorkflow":
        """Run community detection on the built network.

        Returns:
            self for method chaining
        """
        if "network" not in self.results:
            raise RuntimeError("Must call build_network() before detect_communities()")

        from .analysis.community import detect_communities, evaluate_communities, compare_community_methods

        network = self.results["network"]
        cfg = self.config.community

        if cfg.compare_methods:
            comparison = compare_community_methods(network.graph, methods=cfg.methods_to_compare)
            self.results["community_comparison"] = comparison
        else:
            kwargs = {}
            if cfg.resolution != 1.0:
                kwargs["resolution"] = cfg.resolution
            if cfg.n_communities is not None:
                kwargs["k"] = cfg.n_communities

            communities = detect_communities(network.graph, method=cfg.method, **kwargs)
            self.results["communities"] = communities

            # Convert dict to list of lists for evaluation
            comm_lists: Dict[int, List[str]] = {}
            for node, comm_id in communities.items():
                comm_lists.setdefault(comm_id, []).append(node)
            comm_list = list(comm_lists.values())

            evaluation = evaluate_communities(network.graph, comm_list)
            self.results["community_evaluation"] = evaluation

        self._steps_completed.append("detect_communities")
        logger.info("Community detection complete")
        return self

    def analyze_metrics(self) -> "NetworkWorkflow":
        """Compute comprehensive network metrics.

        Returns:
            self for method chaining
        """
        if "network" not in self.results:
            raise RuntimeError("Must call build_network() before analyze_metrics()")

        from .analysis.graph import network_metrics, centrality_measures

        network = self.results["network"]
        self.results["metrics"] = network_metrics(network.graph)
        self.results["centrality"] = centrality_measures(network.graph)

        self._steps_completed.append("analyze_metrics")
        logger.info("Network metrics analysis complete")
        return self

    def run_pathway_enrichment(
        self,
        gene_list: List[str],
        pathway_data: Optional[Dict[str, Any]] = None,
        pathway_network: Optional[Any] = None,
    ) -> "NetworkWorkflow":
        """Run pathway enrichment analysis.

        Args:
            gene_list: Genes to test for enrichment
            pathway_data: Raw pathway data dict
            pathway_network: Pre-built PathwayNetwork

        Returns:
            self for method chaining
        """
        from .analysis.pathway import pathway_enrichment, load_pathway_database

        cfg = self.config.pathway

        if pathway_network is None and pathway_data is not None:
            pathway_network = load_pathway_database(pathway_data)
            if cfg.min_pathway_size > 0 or cfg.max_pathway_size is not None:
                pathway_network = pathway_network.filter_pathways_by_size(
                    min_size=cfg.min_pathway_size, max_size=cfg.max_pathway_size
                )

        if pathway_network is None:
            raise ValueError("Must provide pathway_data or pathway_network")

        results = pathway_enrichment(
            gene_list=gene_list,
            pathway_network=pathway_network,
            background_genes=cfg.background_genes,
            method=cfg.method,
            correction=cfg.correction,
            min_overlap=cfg.min_overlap,
        )

        self.results["pathway_enrichment"] = results
        self._steps_completed.append("pathway_enrichment")
        logger.info(f"Pathway enrichment complete: {len(results)} pathways tested")
        return self

    def export_results(self, output_dir: Optional[str] = None) -> Dict[str, str]:
        """Export analysis results to files.

        Args:
            output_dir: Output directory (overrides config)

        Returns:
            Dictionary mapping result type to file path
        """
        import json

        out_dir = Path(output_dir or self.config.output_dir or "output/networks")
        out_dir.mkdir(parents=True, exist_ok=True)

        exported = {}

        if "network" in self.results:
            from .analysis.graph import export_network

            network_file = out_dir / f"network.{self.config.export_format}"
            export_network(
                self.results["network"],
                str(network_file),
                format=self.config.export_format,
            )
            exported["network"] = str(network_file)

        for key in [
            "metrics",
            "centrality",
            "communities",
            "community_evaluation",
            "community_comparison",
            "pathway_enrichment",
        ]:
            if key in self.results:
                result_file = out_dir / f"{key}.json"
                with open(result_file, "w") as f:
                    json.dump(self.results[key], f, indent=2, default=str)
                exported[key] = str(result_file)

        logger.info(f"Exported {len(exported)} result files to {out_dir}")
        return exported

    def summary(self) -> Dict[str, Any]:
        """Get workflow execution summary.

        Returns:
            Summary of completed steps and key results
        """
        summary: Dict[str, Any] = {
            "steps_completed": self._steps_completed,
            "n_steps": len(self._steps_completed),
        }

        if "network" in self.results:
            net = self.results["network"]
            summary["network"] = {
                "n_nodes": net.number_of_nodes(),
                "n_edges": net.number_of_edges(),
                "directed": net.is_directed(),
            }

        if "metrics" in self.results:
            m = self.results["metrics"]
            summary["metrics"] = {
                "density": m.get("density"),
                "avg_degree": m.get("avg_degree"),
                "num_components": m.get("num_components"),
            }

        if "communities" in self.results:
            comms = self.results["communities"]
            n_comms = len(set(comms.values())) if comms else 0
            summary["communities"] = {"n_communities": n_comms}

        if "pathway_enrichment" in self.results:
            pe = self.results["pathway_enrichment"]
            summary["pathway_enrichment"] = {"n_pathways_tested": len(pe)}

        return summary
