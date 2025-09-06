"""Gene regulatory network inference and analysis."""

from __future__ import annotations

import warnings
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd

from .graph import BiologicalNetwork, create_network


class GeneRegulatoryNetwork:
    """Gene regulatory network representation."""

    def __init__(self, name: str = "grn"):
        """Initialize gene regulatory network.

        Args:
            name: Name of the regulatory network
        """
        self.name = name
        self.regulations: List[Tuple[str, str, Dict[str, Any]]] = []
        self.genes: Set[str] = set()
        self.transcription_factors: Set[str] = set()
        self.gene_metadata: Dict[str, Dict[str, Any]] = {}

    def add_regulation(
        self,
        regulator: str,
        target: str,
        regulation_type: str = "unknown",
        strength: float = 1.0,
        confidence: float = 1.0,
        **metadata,
    ) -> None:
        """Add regulatory interaction.

        Args:
            regulator: Regulator gene/TF identifier
            target: Target gene identifier
            regulation_type: Type of regulation ("activation", "repression", "unknown")
            strength: Regulatory strength
            confidence: Confidence in the interaction
            **metadata: Additional regulatory metadata
        """
        regulation_data = {"type": regulation_type, "strength": strength, "confidence": confidence, **metadata}

        self.regulations.append((regulator, target, regulation_data))
        self.genes.add(regulator)
        self.genes.add(target)

        # Mark regulator as TF if it regulates other genes
        self.transcription_factors.add(regulator)

    def add_gene_metadata(self, gene_id: str, **metadata) -> None:
        """Add metadata for a gene.

        Args:
            gene_id: Gene identifier
            **metadata: Gene metadata (function, location, etc.)
        """
        self.gene_metadata[gene_id] = metadata

    def get_regulators(self, target_gene: str) -> List[Tuple[str, Dict[str, Any]]]:
        """Get regulators of a specific gene.

        Args:
            target_gene: Target gene identifier

        Returns:
            List of (regulator, regulation_data) tuples
        """
        regulators = []

        for reg, target, data in self.regulations:
            if target == target_gene:
                regulators.append((reg, data))

        return regulators

    def get_targets(self, regulator_gene: str) -> List[Tuple[str, Dict[str, Any]]]:
        """Get targets of a specific regulator.

        Args:
            regulator_gene: Regulator gene identifier

        Returns:
            List of (target, regulation_data) tuples
        """
        targets = []

        for reg, target, data in self.regulations:
            if reg == regulator_gene:
                targets.append((target, data))

        return targets

    def create_network(
        self, min_confidence: float = 0.0, regulation_filter: Optional[List[str]] = None
    ) -> BiologicalNetwork:
        """Create BiologicalNetwork from regulatory interactions.

        Args:
            min_confidence: Minimum confidence for including regulations
            regulation_filter: Only include these regulation types

        Returns:
            BiologicalNetwork object (directed)
        """
        # Create directed network
        network = BiologicalNetwork(directed=True)

        # Add all genes as nodes
        for gene in self.genes:
            metadata = self.gene_metadata.get(gene, {})
            is_tf = gene in self.transcription_factors
            network.add_node(gene, is_transcription_factor=is_tf, **metadata)

        # Add regulatory interactions as directed edges
        for regulator, target, data in self.regulations:
            confidence = data.get("confidence", 1.0)
            reg_type = data.get("type", "unknown")
            strength = data.get("strength", 1.0)

            # Apply filters
            if confidence < min_confidence:
                continue

            if regulation_filter and reg_type not in regulation_filter:
                continue

            # Add directed edge from regulator to target
            # Weight combines strength and confidence
            weight = strength * confidence
            network.add_edge(regulator, target, weight=weight)

        return network

    def regulatory_statistics(self) -> Dict[str, Any]:
        """Calculate regulatory network statistics."""
        network = self.create_network()
        from .graph import centrality_measures, network_metrics

        metrics = network_metrics(network)
        centralities = centrality_measures(network)

        # Regulatory-specific statistics
        activation_count = sum(1 for _, _, data in self.regulations if data.get("type") == "activation")
        repression_count = sum(1 for _, _, data in self.regulations if data.get("type") == "repression")

        # Master regulators (high out-degree)
        out_degrees = {}
        for gene in self.genes:
            out_degrees[gene] = len(self.get_targets(gene))

        sorted_regulators = sorted(out_degrees.items(), key=lambda x: x[1], reverse=True)
        master_regulators = [gene for gene, degree in sorted_regulators[:10] if degree > 0]

        return {
            "total_genes": len(self.genes),
            "total_tfs": len(self.transcription_factors),
            "total_regulations": len(self.regulations),
            "activation_regulations": activation_count,
            "repression_regulations": repression_count,
            "master_regulators": master_regulators,
            "network_metrics": metrics,
            "avg_out_degree": np.mean(list(out_degrees.values())),
        }


def infer_grn(
    expression_data: np.ndarray,
    gene_names: List[str],
    method: str = "correlation",
    tf_list: Optional[List[str]] = None,
    threshold: float = 0.7,
) -> GeneRegulatoryNetwork:
    """Infer gene regulatory network from expression data.

    Args:
        expression_data: Expression matrix (samples x genes)
        gene_names: List of gene identifiers
        method: Inference method ("correlation", "mutual_info", "granger")
        tf_list: List of known transcription factors
        threshold: Threshold for regulatory interaction

    Returns:
        GeneRegulatoryNetwork object
    """
    if len(gene_names) != expression_data.shape[1]:
        raise ValueError("Number of gene names must match expression data columns")

    grn = GeneRegulatoryNetwork(name=f"Inferred_GRN_{method}")

    # Set known TFs
    if tf_list:
        for tf in tf_list:
            if tf in gene_names:
                grn.transcription_factors.add(tf)

    if method == "correlation":
        # Pearson correlation-based inference
        corr_matrix = np.corrcoef(expression_data.T)

        for i, reg_gene in enumerate(gene_names):
            for j, target_gene in enumerate(gene_names):
                if i == j:  # Skip self-regulation
                    continue

                correlation = corr_matrix[i, j]

                if abs(correlation) >= threshold:
                    # Determine regulation type based on correlation sign
                    reg_type = "activation" if correlation > 0 else "repression"
                    strength = abs(correlation)

                    # Only add if regulator is a TF (if TF list provided)
                    if tf_list is None or reg_gene in tf_list:
                        grn.add_regulation(
                            reg_gene,
                            target_gene,
                            regulation_type=reg_type,
                            strength=strength,
                            confidence=strength,  # Use correlation as confidence
                            correlation=correlation,
                        )

    elif method == "mutual_info":
        # Mutual information-based inference
        from ..ml.features import _calculate_mi

        # Discretize expression data for MI calculation
        expr_discretized = np.zeros_like(expression_data)
        for j in range(expression_data.shape[1]):
            gene_expr = expression_data[:, j]
            # Simple binning into 3 levels (low, medium, high)
            low_thresh = np.percentile(gene_expr, 33)
            high_thresh = np.percentile(gene_expr, 67)

            expr_discretized[:, j] = np.where(gene_expr <= low_thresh, 0, np.where(gene_expr <= high_thresh, 1, 2))

        for i, reg_gene in enumerate(gene_names):
            for j, target_gene in enumerate(gene_names):
                if i == j:
                    continue

                reg_expr = expr_discretized[:, i].astype(int)
                target_expr = expr_discretized[:, j].astype(int)

                mi = _calculate_mi(reg_expr, target_expr)

                if mi >= threshold:
                    # Determine regulation type (simplified)
                    corr = np.corrcoef(expression_data[:, i], expression_data[:, j])[0, 1]
                    reg_type = "activation" if corr > 0 else "repression"

                    if tf_list is None or reg_gene in tf_list:
                        grn.add_regulation(
                            reg_gene,
                            target_gene,
                            regulation_type=reg_type,
                            strength=mi,
                            confidence=mi,
                            mutual_information=mi,
                        )

    elif method == "granger":
        # Simplified Granger causality (requires time series data)
        warnings.warn("Granger causality inference is simplified - requires proper time series")

        # For each potential regulator-target pair
        for i, reg_gene in enumerate(gene_names):
            for j, target_gene in enumerate(gene_names):
                if i == j:
                    continue

                reg_expr = expression_data[:, i]
                target_expr = expression_data[:, j]

                # Simple lagged correlation as proxy for Granger causality
                if len(reg_expr) > 1:
                    # Lag-1 correlation
                    reg_lagged = reg_expr[:-1]
                    target_current = target_expr[1:]

                    if len(reg_lagged) > 2:
                        lag_corr = np.corrcoef(reg_lagged, target_current)[0, 1]

                        if not np.isnan(lag_corr) and abs(lag_corr) >= threshold:
                            reg_type = "activation" if lag_corr > 0 else "repression"

                            if tf_list is None or reg_gene in tf_list:
                                grn.add_regulation(
                                    reg_gene,
                                    target_gene,
                                    regulation_type=reg_type,
                                    strength=abs(lag_corr),
                                    confidence=abs(lag_corr),
                                    lag_correlation=lag_corr,
                                )

    else:
        raise ValueError(f"Unknown inference method: {method}")

    return grn


def regulatory_motifs(
    grn: GeneRegulatoryNetwork, motif_types: Optional[List[str]] = None
) -> Dict[str, List[Dict[str, Any]]]:
    """Identify regulatory motifs in the network.

    Args:
        grn: Gene regulatory network
        motif_types: Types of motifs to search for

    Returns:
        Dictionary of motif_type -> list of motif instances
    """
    if motif_types is None:
        motif_types = ["feed_forward_loop", "feedback_loop", "bifan"]

    motifs = {motif_type: [] for motif_type in motif_types}
    genes = list(grn.genes)

    # Build adjacency for faster lookup
    adjacency = {}
    for reg, target, data in grn.regulations:
        if reg not in adjacency:
            adjacency[reg] = []
        adjacency[reg].append((target, data))

    # Search for motifs
    for motif_type in motif_types:

        if motif_type == "feed_forward_loop":
            # A -> B, A -> C, B -> C
            for a in genes:
                a_targets = [t for t, _ in adjacency.get(a, [])]

                for b in a_targets:
                    b_targets = [t for t, _ in adjacency.get(b, [])]

                    # Find common targets of A and B
                    common_targets = set(a_targets).intersection(set(b_targets))

                    for c in common_targets:
                        if c != a and c != b:  # Avoid self-loops
                            motifs["feed_forward_loop"].append({"regulator1": a, "intermediate": b, "target": c})

        elif motif_type == "feedback_loop":
            # A -> B -> A (two-gene feedback loop)
            for a in genes:
                a_targets = [t for t, _ in adjacency.get(a, [])]

                for b in a_targets:
                    b_targets = [t for t, _ in adjacency.get(b, [])]

                    if a in b_targets:
                        motifs["feedback_loop"].append({"gene1": a, "gene2": b})

        elif motif_type == "bifan":
            # A -> C, A -> D, B -> C, B -> D
            for a in genes:
                a_targets = set(t for t, _ in adjacency.get(a, []))

                for b in genes:
                    if b == a:
                        continue

                    b_targets = set(t for t, _ in adjacency.get(b, []))

                    # Find common targets
                    common_targets = a_targets.intersection(b_targets)

                    if len(common_targets) >= 2:
                        # Take first two common targets
                        targets_list = list(common_targets)[:2]
                        motifs["bifan"].append(
                            {"regulator1": a, "regulator2": b, "target1": targets_list[0], "target2": targets_list[1]}
                        )

    return motifs


def pathway_regulation_analysis(
    grn: GeneRegulatoryNetwork, pathway_genes: Dict[str, List[str]], min_confidence: float = 0.5
) -> Dict[str, Any]:
    """Analyze regulatory control of biological pathways.

    Args:
        grn: Gene regulatory network
        pathway_genes: Dictionary of pathway_id -> gene_list
        min_confidence: Minimum confidence for regulations

    Returns:
        Pathway regulation analysis results
    """
    results = {}

    for pathway_id, genes in pathway_genes.items():
        pathway_results = {
            "pathway_size": len(genes),
            "regulated_genes": 0,
            "internal_regulations": 0,
            "external_regulators": set(),
            "pathway_regulators": set(),
        }

        genes_set = set(genes)

        # Analyze regulations for each gene in pathway
        for gene in genes:
            regulators = grn.get_regulators(gene)

            # Filter by confidence
            high_conf_regulators = [
                (reg, data) for reg, data in regulators if data.get("confidence", 1.0) >= min_confidence
            ]

            if high_conf_regulators:
                pathway_results["regulated_genes"] += 1

            for reg, data in high_conf_regulators:
                if reg in genes_set:
                    # Internal regulation (within pathway)
                    pathway_results["internal_regulations"] += 1
                    pathway_results["pathway_regulators"].add(reg)
                else:
                    # External regulation (from outside pathway)
                    pathway_results["external_regulators"].add(reg)

        # Convert sets to lists for JSON serialization
        pathway_results["external_regulators"] = list(pathway_results["external_regulators"])
        pathway_results["pathway_regulators"] = list(pathway_results["pathway_regulators"])

        # Calculate regulation statistics
        total_genes = len(genes)
        if total_genes > 0:
            pathway_results["regulation_coverage"] = pathway_results["regulated_genes"] / total_genes
            pathway_results["internal_regulation_ratio"] = pathway_results["internal_regulations"] / max(
                1, pathway_results["regulated_genes"]
            )

        results[pathway_id] = pathway_results

    return results
