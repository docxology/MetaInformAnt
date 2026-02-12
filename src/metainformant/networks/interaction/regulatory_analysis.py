"""Regulatory network analysis: inference, motifs, cascades, validation.

This module provides standalone analysis functions that operate on
GeneRegulatoryNetwork instances or expression data, including GRN inference,
motif detection, cascade discovery, regulation validation, and pathway
regulation analysis.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Optional dependencies
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False

from metainformant.networks.interaction.regulatory_core import GeneRegulatoryNetwork


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

    # Add all genes as nodes so they appear even without edges
    for gene in gene_names:
        network.graph.add_node(gene)

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

        # Pre-register all TF genes
        for tf in tf_genes:
            if tf not in network.tf_targets:
                network.tf_targets[tf] = []

        # Add regulatory edges
        for i, tf in enumerate(tf_genes):
            if tf not in gene_names:
                continue
            tf_idx = gene_names.index(tf)

            for j, target in enumerate(gene_names):
                if tf == target:  # Skip self-regulation
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
        import numpy as np
        from sklearn.feature_selection import mutual_info_regression

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


def detect_regulatory_cascades(
    grn: GeneRegulatoryNetwork,
    max_length: int = 5,
    min_confidence: float = 0.0,
) -> List[Dict[str, Any]]:
    """Detect regulatory cascades in a gene regulatory network.

    Args:
        grn: GeneRegulatoryNetwork instance
        max_length: Maximum cascade length
        min_confidence: Minimum confidence threshold

    Returns:
        List of cascade dictionaries
    """
    cascades: List[Dict[str, Any]] = []
    if not HAS_NETWORKX or grn.graph is None:
        return cascades

    for source in grn.get_transcription_factors():
        for target in grn.graph.nodes():
            if source != target:
                try:
                    paths = list(nx.all_simple_paths(grn.graph, source, target, cutoff=max_length))
                    for path in paths:
                        if len(path) >= 3:
                            # Check confidence along path
                            valid = True
                            if min_confidence > 0:
                                for i in range(len(path) - 1):
                                    edge_data = grn.graph.get_edge_data(path[i], path[i + 1])
                                    conf = (
                                        edge_data.get("confidence", edge_data.get("weight", 1.0)) if edge_data else 0.0
                                    )
                                    if conf < min_confidence:
                                        valid = False
                                        break
                            if valid:
                                cascades.append(
                                    {
                                        "genes": path,
                                        "length": len(path) - 1,
                                        "type": "regulatory_cascade",
                                    }
                                )
                except (nx.NetworkXNoPath, nx.NodeNotFound):
                    continue
    return cascades


def validate_regulation(
    grn: GeneRegulatoryNetwork,
    regulator: str,
    target: str,
    min_confidence: float = 0.0,
) -> Dict[str, Any]:
    """Validate a specific regulatory interaction.

    Args:
        grn: GeneRegulatoryNetwork instance
        regulator: Regulator gene
        target: Target gene
        min_confidence: Minimum confidence threshold

    Returns:
        Validation result dictionary
    """
    result: Dict[str, Any] = {
        "regulator": regulator,
        "target": target,
        "exists": False,
        "confidence": 0.0,
        "valid": False,
    }

    if HAS_NETWORKX and grn.graph is not None and grn.graph.has_edge(regulator, target):
        edge_data = grn.graph.get_edge_data(regulator, target)
        confidence = edge_data.get("confidence", edge_data.get("weight", 1.0))
        result["exists"] = True
        result["confidence"] = confidence
        result["valid"] = confidence >= min_confidence
        result["edge_data"] = dict(edge_data) if edge_data else {}

    return result


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
