"""Pathway network analysis for biological systems."""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from .graph import BiologicalNetwork, create_network


class PathwayNetwork:
    """Biological pathway network representation and analysis.
    
    Manages pathways (gene sets) and their relationships, enabling pathway
    overlap analysis, enrichment analysis, and network construction based
    on pathway similarity.
    
    Attributes:
        name: Name identifier for this pathway network
        pathways: Dictionary mapping pathway_id -> set of gene IDs
        gene_pathways: Dictionary mapping gene_id -> set of pathway IDs
        pathway_metadata: Dictionary mapping pathway_id -> metadata dictionary
        
    Examples:
        >>> pn = PathwayNetwork(name="KEGG_pathways")
        >>> pn.add_pathway("path:00010", ["GENE1", "GENE2", "GENE3"],
        ...                metadata={"name": "Glycolysis"})
        >>> pn.get_pathway_genes("path:00010")
        {'GENE1', 'GENE2', 'GENE3'}
    """

    def __init__(self, name: str = "pathway_network"):
        """Initialize pathway network.

        Args:
            name: Name identifier for this pathway network
        """
        self.name = name
        self.pathways: Dict[str, Set[str]] = {}  # pathway_id -> gene_set
        self.gene_pathways: Dict[str, Set[str]] = defaultdict(set)  # gene -> pathway_ids
        self.pathway_metadata: Dict[str, Dict[str, Any]] = {}

    def add_pathway(self, pathway_id: str, genes: List[str], metadata: Optional[Dict[str, Any]] = None) -> None:
        """Add a pathway with its associated genes.
        
        Adds a biological pathway to the network and automatically updates
        gene-pathway mappings. Pathway metadata can include name, description,
        or other annotations.

        Args:
            pathway_id: Unique pathway identifier (e.g., "KEGG_PATH:00010")
            genes: List of gene identifiers in the pathway
            metadata: Optional dictionary of pathway metadata (e.g.,
                {"name": "Glycolysis", "source": "KEGG", "category": "Metabolism"})
                
        Examples:
            >>> pn = PathwayNetwork()
            >>> pn.add_pathway(
            ...     "path:00010",
            ...     ["GENE1", "GENE2", "GENE3"],
            ...     metadata={"name": "Glycolysis", "source": "KEGG"}
            ... )
            >>> len(pn.get_pathway_genes("path:00010"))
            3
        """
        gene_set = set(genes)
        self.pathways[pathway_id] = gene_set

        # Update gene-pathway mappings
        for gene in genes:
            self.gene_pathways[gene].add(pathway_id)

        # Store metadata
        if metadata:
            self.pathway_metadata[pathway_id] = metadata

    def get_pathway_genes(self, pathway_id: str) -> Set[str]:
        """Retrieve all genes in a specific pathway.
        
        Args:
            pathway_id: Pathway identifier
            
        Returns:
            Set of gene identifiers in the pathway. Returns empty set if
            pathway doesn't exist.
            
        Examples:
            >>> pn = PathwayNetwork()
            >>> pn.add_pathway("path:00010", ["GENE1", "GENE2"])
            >>> pn.get_pathway_genes("path:00010")
            {'GENE1', 'GENE2'}
        """
        return self.pathways.get(pathway_id, set())

    def get_gene_pathways(self, gene: str) -> Set[str]:
        """Retrieve all pathways containing a specific gene.
        
        Args:
            gene: Gene identifier
            
        Returns:
            Set of pathway identifiers that include this gene. Returns empty
            set if gene is not in any pathway.
            
        Examples:
            >>> pn = PathwayNetwork()
            >>> pn.add_pathway("path:00010", ["GENE1", "GENE2"])
            >>> pn.add_pathway("path:00020", ["GENE1", "GENE3"])
            >>> pn.get_gene_pathways("GENE1")
            {'path:00010', 'path:00020'}
        """
        return self.gene_pathways.get(gene, set())

    def pathway_overlap(self, pathway1: str, pathway2: str) -> Tuple[Set[str], float]:
        """Calculate gene overlap and similarity between two pathways.
        
        Computes both the set of overlapping genes and the Jaccard similarity
        coefficient, which measures pathway similarity as the ratio of
        intersection to union.

        Args:
            pathway1: First pathway identifier
            pathway2: Second pathway identifier

        Returns:
            Tuple containing:
            - overlapping_genes: Set of gene IDs present in both pathways
            - jaccard_index: Jaccard similarity coefficient in [0, 1].
              Formula: |A ∩ B| / |A ∪ B|
              
        Examples:
            >>> pn = PathwayNetwork()
            >>> pn.add_pathway("path1", ["GENE1", "GENE2", "GENE3"])
            >>> pn.add_pathway("path2", ["GENE2", "GENE3", "GENE4"])
            >>> overlap, jaccard = pn.pathway_overlap("path1", "path2")
            >>> overlap
            {'GENE2', 'GENE3'}
            >>> jaccard
            0.5  # 2 overlapping / 4 total unique genes
        """
        genes1 = self.get_pathway_genes(pathway1)
        genes2 = self.get_pathway_genes(pathway2)

        if not genes1 or not genes2:
            return set(), 0.0

        overlap = genes1.intersection(genes2)
        union = genes1.union(genes2)

        jaccard = len(overlap) / len(union) if union else 0.0

        return overlap, jaccard

    def create_pathway_network(self, min_overlap: int = 2, min_jaccard: float = 0.1) -> BiologicalNetwork:
        """Create network where pathways are nodes connected by gene overlap.
        
        Constructs a BiologicalNetwork where each pathway becomes a node,
        and edges connect pathways that share genes. Edge weights represent
        Jaccard similarity. Useful for pathway clustering and hierarchical
        analysis.

        Args:
            min_overlap: Minimum number of shared genes required for edge.
                Higher values create sparser networks.
            min_jaccard: Minimum Jaccard similarity (0-1) required for edge.
                Higher values require more similar pathway composition.

        Returns:
            BiologicalNetwork object with pathways as nodes. Nodes retain
            pathway metadata as attributes. Edge weights are Jaccard indices.
            
        Examples:
            >>> pn = PathwayNetwork()
            >>> pn.add_pathway("path1", ["GENE1", "GENE2", "GENE3"])
            >>> pn.add_pathway("path2", ["GENE2", "GENE3", "GENE4"])
            >>> network = pn.create_pathway_network(min_overlap=1, min_jaccard=0.3)
            >>> network.num_edges()
            1  # path1 and path2 connected
        """
        pathway_ids = list(self.pathways.keys())
        network = create_network(pathway_ids)

        # Add pathway metadata as node attributes
        for pathway_id in pathway_ids:
            metadata = self.pathway_metadata.get(pathway_id, {})
            network.add_node(pathway_id, **metadata)

        # Add edges based on pathway overlap
        for i, p1 in enumerate(pathway_ids):
            for j, p2 in enumerate(pathway_ids):
                if i >= j:  # Avoid duplicates and self-loops
                    continue

                overlap_genes, jaccard = self.pathway_overlap(p1, p2)

                if len(overlap_genes) >= min_overlap and jaccard >= min_jaccard:
                    network.add_edge(p1, p2, weight=jaccard)

        return network


def load_pathway_database(pathway_file: str, format: str = "gmt") -> PathwayNetwork:
    """Load pathway database from file.
    
    Supports multiple file formats for pathway databases including GMT
    (Gene Matrix Transposed) format commonly used for pathway databases
    like MSigDB, and tabular formats.
    
    Args:
        pathway_file: Path to pathway database file
        format: File format:
            - "gmt": GMT format (tab-separated: pathway_id, description, genes...)
            - "csv": CSV format (columns: pathway_id, gene_id, optionally pathway_name)
            - "tsv": TSV format (same as CSV but tab-separated)
            
    Returns:
        PathwayNetwork object with loaded pathways
        
    Raises:
        ValueError: If format is unsupported
        
    Examples:
        >>> # Load GMT format (MSigDB, KEGG, etc.)
        >>> pathway_db = load_pathway_database("kegg_pathways.gmt", format="gmt")
        >>> pathway_db.get_pathway_genes("KEGG_PATH:00010")
        {'GENE1', 'GENE2', ...}
        
        >>> # Load CSV format
        >>> pathway_db = load_pathway_database("custom_pathways.csv", format="csv")
        >>> pathway_db.num_pathways()
        50
    """
    pathway_net = PathwayNetwork()

    if format.lower() == "gmt":
        with open(pathway_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    pathway_id = parts[0]
                    description = parts[1]
                    genes = parts[2:]

                    pathway_net.add_pathway(pathway_id, genes, metadata={"description": description})

    elif format.lower() in ["csv", "tsv"]:
        sep = "," if format.lower() == "csv" else "\t"
        df = pd.read_csv(pathway_file, sep=sep)

        # Assume columns: pathway_id, gene_id (and optionally pathway_name)
        for pathway_id, group in df.groupby("pathway_id"):
            genes = group["gene_id"].tolist()
            metadata = {}

            if "pathway_name" in df.columns:
                metadata["name"] = group["pathway_name"].iloc[0]

            pathway_net.add_pathway(pathway_id, genes, metadata)

    else:
        raise ValueError(f"Unsupported format: {format}")

    return pathway_net


def pathway_enrichment(
    gene_list: List[str], pathway_network: PathwayNetwork, background_genes: Optional[List[str]] = None
) -> List[Dict[str, Any]]:
    """Perform pathway enrichment analysis for a gene list.
    
    Identifies pathways that are overrepresented in a query gene set compared
    to a background set. Uses hypergeometric-like statistics to assess enrichment.
    
    Args:
        gene_list: List of gene identifiers of interest (query set)
        pathway_network: PathwayNetwork object containing pathways to test
        background_genes: Optional background gene set. If None, uses all
            genes present in the pathway network
            
    Returns:
        List of dictionaries, one per enriched pathway, sorted by fold enrichment.
        Each dictionary contains:
        - pathway_id: Pathway identifier
        - pathway_size: Number of genes in pathway (in background)
        - overlap_size: Number of query genes in pathway
        - query_size: Number of query genes in background
        - background_size: Total background genes
        - fold_enrichment: Enrichment ratio
        - p_value: Statistical significance (approximate)
        - overlapping_genes: List of genes in both query and pathway
        - metadata: Pathway metadata dictionary
        
    Examples:
        >>> query_genes = ["GENE1", "GENE2", "GENE5"]
        >>> enrichment = pathway_enrichment(query_genes, pathway_network)
        >>> enrichment[0]["pathway_id"]
        'path:00010'
        >>> enrichment[0]["fold_enrichment"]
        3.5
        
    Note:
        This uses simplified statistics. For rigorous analysis, consider
        using Fisher's exact test or hypergeometric test implementations.
    """
    gene_set = set(gene_list)

    if background_genes is None:
        # Use all genes in pathway network as background
        background_genes = set()
        for genes in pathway_network.pathways.values():
            background_genes.update(genes)
        background_genes = list(background_genes)

    background_set = set(background_genes)

    enrichment_results = []

    for pathway_id, pathway_genes in pathway_network.pathways.items():
        # Intersection with background
        pathway_bg = pathway_genes.intersection(background_set)

        if len(pathway_bg) == 0:
            continue

        # Count overlaps
        overlap = gene_set.intersection(pathway_bg)
        k = len(overlap)  # genes in both query and pathway
        n = len(gene_set.intersection(background_set))  # query genes in background
        K = len(pathway_bg)  # pathway genes in background
        N = len(background_set)  # total background genes

        if k == 0:
            continue

        # Simple enrichment score (more sophisticated stats could be added)
        expected = (n * K) / N if N > 0 else 0
        fold_enrichment = (k / expected) if expected > 0 else float("inf")

        # Simple p-value approximation (hypergeometric would be better)
        p_value = max(0.001, 1.0 / (1.0 + fold_enrichment))

        result = {
            "pathway_id": pathway_id,
            "pathway_size": len(pathway_bg),
            "overlap_size": k,
            "query_size": n,
            "background_size": N,
            "fold_enrichment": fold_enrichment,
            "p_value": p_value,
            "overlapping_genes": list(overlap),
            "metadata": pathway_network.pathway_metadata.get(pathway_id, {}),
        }

        enrichment_results.append(result)

    # Sort by fold enrichment
    enrichment_results.sort(key=lambda x: x["fold_enrichment"], reverse=True)

    return enrichment_results


def network_enrichment_analysis(
    gene_network: BiologicalNetwork, pathway_network: PathwayNetwork, method: str = "overlap"
) -> Dict[str, Any]:
    """Analyze pathway enrichment for network modules.
    
    Identifies pathways enriched in a gene/protein network, either by
    simple overlap with all network genes or by focusing on central
    hub genes.
    
    Args:
        gene_network: Biological network of genes/proteins
        pathway_network: PathwayNetwork object containing pathways to test
        method: Analysis method:
            - "overlap": Enrichment based on all genes in network
            - "centrality": Enrichment based on top 20% most central genes
              (by degree centrality)
              
    Returns:
        Dictionary containing:
        - method: Analysis method used
        - network_size: Number of genes in network (for "overlap")
        - central_genes: List of central gene IDs (for "centrality")
        - enrichment_results: List of enriched pathways with statistics
        
    Raises:
        ValueError: If method is unknown
        
    Examples:
        >>> from metainformant.networks import create_network
        >>> network = create_network(["GENE1", "GENE2", "GENE3"])
        >>> enrichment = network_enrichment_analysis(network, pathway_network)
        >>> enrichment["method"]
        'overlap'
        >>> len(enrichment["enrichment_results"]) > 0
        True
    """
    if method == "overlap":
        # Simple overlap-based enrichment
        network_genes = list(gene_network.nodes)
        enrichment = pathway_enrichment(network_genes, pathway_network)

        return {"method": method, "network_size": len(network_genes), "enrichment_results": enrichment}

    elif method == "centrality":
        # Enrichment based on central genes in network
        from .graph import centrality_measures

        centralities = centrality_measures(gene_network)

        # Get top 20% most central genes by degree centrality
        degree_cents = centralities["degree"]
        sorted_genes = sorted(degree_cents.items(), key=lambda x: x[1], reverse=True)
        top_genes = [gene for gene, _ in sorted_genes[: len(sorted_genes) // 5]]

        enrichment = pathway_enrichment(top_genes, pathway_network)

        return {"method": method, "central_genes": top_genes, "enrichment_results": enrichment}

    else:
        raise ValueError(f"Unknown method: {method}")
