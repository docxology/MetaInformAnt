"""Gene regulatory network inference and analysis."""

from __future__ import annotations

import warnings
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np
import pandas as pd

from .graph import BiologicalNetwork, create_network


class GeneRegulatoryNetwork:
    """Gene regulatory network (GRN) representation and analysis.
    
    Manages regulatory interactions between transcription factors and target
    genes, including activation/repression relationships, regulatory strengths,
    and confidence scores. Supports network inference and motif detection.
    
    Attributes:
        name: Name identifier for this regulatory network
        regulations: List of (regulator, target, regulation_data) tuples
        genes: Set of all gene identifiers in the network
        transcription_factors: Set of gene IDs that act as transcription factors
        gene_metadata: Dictionary mapping gene_id -> metadata dictionary
        
    Examples:
        >>> grn = GeneRegulatoryNetwork(name="Inferred_GRN")
        >>> grn.add_regulation(
        ...     regulator="TF1",
        ...     target="GENE1",
        ...     regulation_type="activation",
        ...     strength=0.85,
        ...     confidence=0.9
        ... )
        >>> network = grn.create_network(min_confidence=0.7)
        >>> stats = grn.regulatory_statistics()
        >>> stats["master_regulators"]
        ['TF1', 'TF2', ...]
    """

    def __init__(self, name: str = "grn"):
        """Initialize gene regulatory network.

        Args:
            name: Name identifier for this regulatory network
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
        """Add a regulatory interaction to the network.
        
        Records a directed regulatory relationship where a regulator (e.g.,
        transcription factor) controls a target gene. Automatically marks
        the regulator as a transcription factor.

        Args:
            regulator: Regulator identifier (transcription factor, miRNA, etc.)
            target: Target gene identifier
            regulation_type: Type of regulation:
                - "activation": Positive regulation (upregulation)
                - "repression": Negative regulation (downregulation)
                - "unknown": Unknown or mixed effect (default)
            strength: Regulatory strength in [0, 1] or arbitrary scale.
                Higher values indicate stronger regulatory effect. Default 1.0.
            confidence: Confidence in the interaction [0, 1]. Default 1.0.
            **metadata: Additional keyword arguments for regulatory metadata
                (e.g., condition="stress", timepoint="t0", binding_site="promoter")
                
        Examples:
            >>> grn = GeneRegulatoryNetwork()
            >>> grn.add_regulation(
            ...     "TF1", "GENE1",
            ...     regulation_type="activation",
            ...     strength=0.85,
            ...     confidence=0.9,
            ...     binding_site="promoter"
            ... )
            >>> "TF1" in grn.transcription_factors
            True
            >>> len(grn.regulations)
            1
        """
        regulation_data = {"type": regulation_type, "strength": strength, "confidence": confidence, **metadata}

        self.regulations.append((regulator, target, regulation_data))
        self.genes.add(regulator)
        self.genes.add(target)

        # Mark regulator as TF if it regulates other genes
        self.transcription_factors.add(regulator)

    def add_gene_metadata(self, gene_id: str, **metadata) -> None:
        """Add or update metadata for a gene.
        
        Stores gene annotations such as function, expression patterns, or
        other properties. Useful for enriching regulatory network analysis.

        Args:
            gene_id: Gene identifier (must match IDs used in regulations)
            **metadata: Keyword arguments for gene metadata. Common fields:
                - function: Gene function description
                - expression_level: Expression level category
                - pathway: Pathways the gene participates in
                - chromosome: Chromosomal location
                
        Examples:
            >>> grn = GeneRegulatoryNetwork()
            >>> grn.add_gene_metadata(
            ...     "GENE1",
            ...     function="cell cycle regulation",
            ...     pathway="p53 signaling"
            ... )
            >>> grn.gene_metadata["GENE1"]["function"]
            'cell cycle regulation'
        """
        self.gene_metadata[gene_id] = metadata

    def get_regulators(self, target_gene: str) -> List[Tuple[str, Dict[str, Any]]]:
        """Retrieve all regulators of a specific target gene.
        
        Returns all transcription factors or other regulators that control
        the specified gene.

        Args:
            target_gene: Target gene identifier

        Returns:
            List of tuples, each containing (regulator_id, regulation_data).
            regulation_data is a dictionary with keys:
            - type: Regulation type ("activation", "repression", "unknown")
            - strength: Regulatory strength
            - confidence: Confidence score
            - Additional metadata fields if provided
            
        Examples:
            >>> grn = GeneRegulatoryNetwork()
            >>> grn.add_regulation("TF1", "GENE1", regulation_type="activation")
            >>> grn.add_regulation("TF2", "GENE1", regulation_type="repression")
            >>> regulators = grn.get_regulators("GENE1")
            >>> len(regulators)
            2
            >>> [r[0] for r in regulators]
            ['TF1', 'TF2']
        """
        regulators = []

        for reg, target, data in self.regulations:
            if target == target_gene:
                regulators.append((reg, data))

        return regulators

    def get_targets(self, regulator_gene: str) -> List[Tuple[str, Dict[str, Any]]]:
        """Retrieve all target genes of a specific regulator.
        
        Returns all genes that are regulated by the specified transcription
        factor or other regulator.

        Args:
            regulator_gene: Regulator gene identifier

        Returns:
            List of tuples, each containing (target_gene_id, regulation_data).
            regulation_data structure same as in get_regulators().
            
        Examples:
            >>> grn = GeneRegulatoryNetwork()
            >>> grn.add_regulation("TF1", "GENE1", regulation_type="activation")
            >>> grn.add_regulation("TF1", "GENE2", regulation_type="repression")
            >>> targets = grn.get_targets("TF1")
            >>> len(targets)
            2
            >>> [t[0] for t in targets]
            ['GENE1', 'GENE2']
        """
        targets = []

        for reg, target, data in self.regulations:
            if reg == regulator_gene:
                targets.append((target, data))

        return targets

    def create_network(
        self, min_confidence: float = 0.0, regulation_filter: Optional[List[str]] = None
    ) -> BiologicalNetwork:
        """Convert GeneRegulatoryNetwork to BiologicalNetwork object.
        
        Creates a directed graph representation of regulatory interactions
        with filtering by confidence and regulation type. Edge weights
        combine regulatory strength and confidence.

        Args:
            min_confidence: Minimum confidence score (0-1) for including
                regulations. Default 0.0 includes all regulations.
            regulation_filter: Optional list of regulation types to include.
                If provided, only regulations matching these types are included.
                Valid values: "activation", "repression", "unknown"
                
        Returns:
            BiologicalNetwork object (directed=True) with:
            - Nodes: All genes (marked as TFs if applicable)
            - Edges: Filtered regulatory interactions
            - Edge weights: strength × confidence
            - Node attributes: Gene metadata if available
            
        Examples:
            >>> grn = GeneRegulatoryNetwork()
            >>> grn.add_regulation("TF1", "GENE1", regulation_type="activation", confidence=0.8)
            >>> grn.add_regulation("TF1", "GENE2", regulation_type="repression", confidence=0.3)
            >>> network = grn.create_network(
            ...     min_confidence=0.5,
            ...     regulation_filter=["activation"]
            ... )
            >>> network.num_edges()
            1  # Only high-confidence activation included
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
        """Calculate comprehensive regulatory network statistics.
        
        Computes network topology, regulatory type counts, and identifies
        master regulators (genes with high out-degree).
        
        Returns:
            Dictionary containing:
            - total_genes: Number of unique genes in network
            - total_tfs: Number of transcription factors
            - total_regulations: Number of regulatory interactions
            - activation_regulations: Count of activation interactions
            - repression_regulations: Count of repression interactions
            - master_regulators: List of top 10 regulators by out-degree
            - network_metrics: Basic network topology metrics
            - avg_out_degree: Average number of targets per regulator
            
        Examples:
            >>> stats = grn.regulatory_statistics()
            >>> stats["total_tfs"]
            50
            >>> stats["master_regulators"]
            ['TF1', 'TF2', 'TF3']
        """
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
    """Infer gene regulatory network from gene expression data.
    
    Constructs a regulatory network by identifying relationships between
    regulators and target genes based on expression patterns. Multiple
    inference methods are available.
    
    Args:
        expression_data: Expression matrix of shape (n_samples, n_genes).
            Each column corresponds to one gene's expression profile.
        gene_names: List of gene identifiers corresponding to columns
            in expression_data. Must have length equal to expression_data.shape[1]
        method: Inference method:
            - "correlation": Pearson correlation between regulator and target
              expression profiles
            - "mutual_info": Mutual information between discretized expressions
            - "granger": Simplified Granger causality (requires time series data)
        tf_list: Optional list of known transcription factor gene names.
            If provided, only these genes can act as regulators. If None,
            all genes can be regulators.
        threshold: Minimum correlation/mutual information for regulatory
            interaction. Interactions with |correlation| >= threshold are included.
            
    Returns:
        GeneRegulatoryNetwork object with inferred regulatory interactions.
        Regulation type is determined by correlation sign (positive=activation,
        negative=repression). Confidence equals correlation/mutual info strength.
        
    Raises:
        ValueError: If number of gene_names doesn't match expression data columns
        ValueError: If method is unknown
        
    Examples:
        >>> expression = np.random.randn(100, 200)  # 100 samples, 200 genes
        >>> gene_names = [f"GENE_{i}" for i in range(200)]
        >>> tf_list = [f"GENE_{i}" for i in range(50)]  # First 50 are TFs
        >>> grn = infer_grn(
        ...     expression, gene_names,
        ...     method="correlation",
        ...     tf_list=tf_list,
        ...     threshold=0.75
        ... )
        >>> grn.regulatory_statistics()["total_regulations"]
        320
        
    References:
        Bansal, M., et al. (2007). How to infer gene networks from expression
        profiles. Molecular Systems Biology, 3(1), 78.
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
        try:
            from ..information.syntactic import mutual_information
        except ImportError:
            try:
                from ..ml.features import _calculate_mi as mutual_information
            except ImportError:
                raise ImportError(
                    "Mutual information calculation requires either information or ml.features module. "
                    "Please ensure one is available."
                )

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

                reg_expr = expr_discretized[:, i].astype(int).tolist()
                target_expr = expr_discretized[:, j].astype(int).tolist()

                mi = mutual_information(reg_expr, target_expr)

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
    """Identify regulatory network motifs (small recurring patterns).
    
    Regulatory motifs are small subnetworks that occur frequently and
    often have specific biological functions. Common motifs include
    feed-forward loops and feedback loops.
    
    Args:
        grn: Gene regulatory network to analyze
        motif_types: List of motif types to search for. If None, searches
            for all supported types. Supported types:
            - "feed_forward_loop": A -> B -> C, A -> C (common in GRNs)
            - "feedback_loop": A -> B -> A (two-gene feedback)
            - "bifan": Two regulators both target two common targets
            
    Returns:
        Dictionary mapping motif_type to list of motif instances.
        Each motif instance is a dictionary with relevant node identifiers.
        
    Examples:
        >>> grn = GeneRegulatoryNetwork()
        >>> grn.add_regulation("TF1", "TF2", regulation_type="activation")
        >>> grn.add_regulation("TF2", "GENE1", regulation_type="activation")
        >>> grn.add_regulation("TF1", "GENE1", regulation_type="activation")
        >>> motifs = regulatory_motifs(grn, motif_types=["feed_forward_loop"])
        >>> len(motifs["feed_forward_loop"])
        1
        
    References:
        Alon, U. (2007). Network motifs: theory and experimental approaches.
        Nature Reviews Genetics, 8(6), 450-461.
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
    """Analyze regulatory control structure of biological pathways.
    
    Examines how pathways are regulated by transcription factors,
    distinguishing between internal regulation (within pathway) and
    external regulation (from outside pathway).
    
    Args:
        grn: Gene regulatory network
        pathway_genes: Dictionary mapping pathway_id to list of gene IDs
            in that pathway
        min_confidence: Minimum confidence threshold for regulatory
            interactions to include in analysis
            
    Returns:
        Dictionary mapping pathway_id to analysis results. Each result
        contains:
        - pathway_size: Number of genes in pathway
        - regulated_genes: Number of genes with regulatory inputs
        - internal_regulations: Count of regulations within pathway
        - external_regulators: List of regulators outside pathway
        - pathway_regulators: List of regulators within pathway
        - regulation_coverage: Fraction of genes with regulations
        - internal_regulation_ratio: Ratio of internal to total regulations
        
    Examples:
        >>> pathway_genes = {
        ...     "pathway1": ["GENE1", "GENE2", "GENE3"],
        ...     "pathway2": ["GENE4", "GENE5"]
        ... }
        >>> analysis = pathway_regulation_analysis(grn, pathway_genes, min_confidence=0.7)
        >>> analysis["pathway1"]["regulation_coverage"]
        0.666...
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


def detect_regulatory_cascades(
    grn: GeneRegulatoryNetwork, max_length: int = 5, min_confidence: float = 0.5
) -> List[Dict[str, Any]]:
    """Detect regulatory cascades (multi-level regulation chains).
    
    Identifies chains of regulation where one regulator controls another,
    which in turn controls downstream targets. Cascades reveal hierarchical
    regulatory structure.
    
    Args:
        grn: Gene regulatory network
        max_length: Maximum cascade length to detect (default 5)
        min_confidence: Minimum confidence for regulatory interactions
        
    Returns:
        List of cascade dictionaries, each containing:
        - cascade_id: Unique identifier
        - length: Number of regulatory steps
        - regulators: List of regulator IDs in cascade order
        - target: Final target gene
        - avg_confidence: Average confidence across cascade
        
    Examples:
        >>> grn = GeneRegulatoryNetwork()
        >>> grn.add_regulation("TF1", "TF2", confidence=0.8)
        >>> grn.add_regulation("TF2", "GENE1", confidence=0.9)
        >>> cascades = detect_regulatory_cascades(grn)
        >>> len(cascades) > 0
        True
    """
    cascades = []
    cascade_id = 0
    
    # Build adjacency list for efficient traversal
    adjacency = {}
    for reg, target, data in grn.regulations:
        if data.get("confidence", 1.0) < min_confidence:
            continue
        if reg not in adjacency:
            adjacency[reg] = []
        adjacency[reg].append((target, data))
    
    # Find cascades starting from each TF
    for start_tf in grn.transcription_factors:
        if start_tf not in adjacency:
            continue
        
        # BFS to find cascades
        queue = [(start_tf, [start_tf], 0.0)]  # (current, path, avg_conf)
        
        while queue:
            current, path, avg_conf = queue.pop(0)
            
            if len(path) > max_length:
                continue
            
            if current in adjacency:
                for target, data in adjacency[current]:
                    if target in path:  # Avoid cycles
                        continue
                    
                    conf = data.get("confidence", 1.0)
                    new_avg = (avg_conf * (len(path) - 1) + conf) / len(path)
                    new_path = path + [target]
                    
                    # Record cascade
                    if len(new_path) >= 2:
                        cascades.append(
                            {
                                "cascade_id": f"cascade_{cascade_id}",
                                "length": len(new_path) - 1,
                                "regulators": new_path[:-1],
                                "target": target,
                                "avg_confidence": new_avg,
                            }
                        )
                        cascade_id += 1
                    
                    # Continue if target is also a regulator
                    if target in grn.transcription_factors and len(new_path) < max_length:
                        queue.append((target, new_path, new_avg))
    
    return cascades


def validate_regulation(
    grn: GeneRegulatoryNetwork, regulator: str, target: str, min_confidence: float = 0.5
) -> Dict[str, Any]:
    """Validate a regulatory interaction and provide evidence.
    
    Checks if a regulation exists and provides validation metrics
    including confidence, evidence types, and network context.
    
    Args:
        grn: Gene regulatory network
        regulator: Regulator identifier
        target: Target gene identifier
        min_confidence: Minimum confidence threshold
        
    Returns:
        Dictionary containing:
        - exists: Whether regulation exists
        - confidence: Confidence score if exists
        - regulation_type: Type of regulation (activation/repression)
        - strength: Regulatory strength
        - evidence: List of evidence types
        - network_support: Whether supported by network structure
        
    Examples:
        >>> grn = GeneRegulatoryNetwork()
        >>> grn.add_regulation("TF1", "GENE1", confidence=0.8)
        >>> validation = validate_regulation(grn, "TF1", "GENE1")
        >>> validation["exists"]
        True
    """
    # Find regulation
    regulation = None
    for reg, tgt, data in grn.regulations:
        if reg == regulator and tgt == target:
            regulation = data
            break
    
    if not regulation:
        return {"exists": False}
    
    confidence = regulation.get("confidence", 0.0)
    if confidence < min_confidence:
        return {"exists": True, "confidence": confidence, "below_threshold": True}
    
    # Check network support (common neighbors, etc.)
    reg_targets = set(tgt for reg, tgt, _ in grn.regulations if reg == regulator)
    target_regulators = set(reg for reg, tgt, _ in grn.regulations if tgt == target)
    common_neighbors = reg_targets.intersection(target_regulators)
    
    return {
        "exists": True,
        "confidence": confidence,
        "regulation_type": regulation.get("type", "unknown"),
        "strength": regulation.get("strength", 1.0),
        "evidence": regulation.get("evidence_types", []),
        "network_support": len(common_neighbors) > 0,
        "common_neighbors": len(common_neighbors),
    }
