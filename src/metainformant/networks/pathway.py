"""Pathway network analysis for biological systems."""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from .graph import BiologicalNetwork, create_network


class PathwayNetwork:
    """Biological pathway network representation."""

    def __init__(self, name: str = "pathway_network"):
        """Initialize pathway network.

        Args:
            name: Name of the pathway network
        """
        self.name = name
        self.pathways: Dict[str, Set[str]] = {}  # pathway_id -> gene_set
        self.gene_pathways: Dict[str, Set[str]] = defaultdict(set)  # gene -> pathway_ids
        self.pathway_metadata: Dict[str, Dict[str, Any]] = {}

    def add_pathway(self, pathway_id: str, genes: List[str], metadata: Optional[Dict[str, Any]] = None) -> None:
        """Add a pathway with its associated genes.

        Args:
            pathway_id: Unique pathway identifier
            genes: List of genes in the pathway
            metadata: Optional pathway metadata (name, description, etc.)
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
        """Get genes in a specific pathway."""
        return self.pathways.get(pathway_id, set())

    def get_gene_pathways(self, gene: str) -> Set[str]:
        """Get pathways containing a specific gene."""
        return self.gene_pathways.get(gene, set())

    def pathway_overlap(self, pathway1: str, pathway2: str) -> Tuple[Set[str], float]:
        """Calculate overlap between two pathways.

        Args:
            pathway1: First pathway ID
            pathway2: Second pathway ID

        Returns:
            Tuple of (overlapping_genes, jaccard_index)
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
        """Create network of pathways based on gene overlap.

        Args:
            min_overlap: Minimum number of overlapping genes
            min_jaccard: Minimum Jaccard index for edge creation

        Returns:
            Network where nodes are pathways and edges represent overlap
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

    Args:
        pathway_file: Path to pathway file
        format: File format ("gmt", "csv", "tsv")

    Returns:
        PathwayNetwork object
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
    """Perform pathway enrichment analysis.

    Args:
        gene_list: List of genes of interest
        pathway_network: PathwayNetwork object
        background_genes: Background gene set (optional)

    Returns:
        List of enrichment results
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
    """Analyze enrichment of network modules in pathways.

    Args:
        gene_network: Network of genes/proteins
        pathway_network: PathwayNetwork object
        method: Analysis method ("overlap", "centrality")

    Returns:
        Dictionary of enrichment results
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
