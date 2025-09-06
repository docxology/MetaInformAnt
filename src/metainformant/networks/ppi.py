"""Protein-protein interaction network analysis."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from .graph import BiologicalNetwork, create_network


class ProteinNetwork:
    """Protein-protein interaction network representation."""

    def __init__(self, name: str = "protein_network"):
        """Initialize protein network.

        Args:
            name: Name of the protein network
        """
        self.name = name
        self.interactions: List[Tuple[str, str, Dict[str, Any]]] = []
        self.proteins: Set[str] = set()
        self.protein_metadata: Dict[str, Dict[str, Any]] = {}

    def add_interaction(
        self,
        protein1: str,
        protein2: str,
        confidence: float = 1.0,
        evidence_types: Optional[List[str]] = None,
        **metadata,
    ) -> None:
        """Add protein-protein interaction.

        Args:
            protein1: First protein identifier
            protein2: Second protein identifier
            confidence: Interaction confidence score (0-1)
            evidence_types: Types of experimental evidence
            **metadata: Additional interaction metadata
        """
        if evidence_types is None:
            evidence_types = []

        interaction_data = {"confidence": confidence, "evidence_types": evidence_types, **metadata}

        self.interactions.append((protein1, protein2, interaction_data))
        self.proteins.add(protein1)
        self.proteins.add(protein2)

    def add_protein_metadata(self, protein_id: str, **metadata) -> None:
        """Add metadata for a protein.

        Args:
            protein_id: Protein identifier
            **metadata: Protein metadata (function, localization, etc.)
        """
        self.protein_metadata[protein_id] = metadata

    def get_interactions(self, protein: str, min_confidence: float = 0.0) -> List[Tuple[str, Dict[str, Any]]]:
        """Get interactions for a specific protein.

        Args:
            protein: Protein identifier
            min_confidence: Minimum confidence threshold

        Returns:
            List of (partner_protein, interaction_data) tuples
        """
        interactions = []

        for p1, p2, data in self.interactions:
            if data.get("confidence", 1.0) < min_confidence:
                continue

            if p1 == protein:
                interactions.append((p2, data))
            elif p2 == protein:
                interactions.append((p1, data))

        return interactions

    def create_network(
        self, min_confidence: float = 0.4, evidence_filter: Optional[List[str]] = None
    ) -> BiologicalNetwork:
        """Create BiologicalNetwork from protein interactions.

        Args:
            min_confidence: Minimum confidence for including interactions
            evidence_filter: Only include interactions with these evidence types

        Returns:
            BiologicalNetwork object
        """
        # Create network with all proteins
        network = create_network(list(self.proteins))

        # Add protein metadata as node attributes
        for protein_id, metadata in self.protein_metadata.items():
            if protein_id in self.proteins:
                network.add_node(protein_id, **metadata)

        # Add interactions as edges
        for p1, p2, data in self.interactions:
            confidence = data.get("confidence", 1.0)
            evidence_types = data.get("evidence_types", [])

            # Apply filters
            if confidence < min_confidence:
                continue

            if evidence_filter and not any(ev in evidence_types for ev in evidence_filter):
                continue

            # Add edge with confidence as weight
            network.add_edge(p1, p2, weight=confidence)

        return network

    def network_statistics(self) -> Dict[str, Any]:
        """Calculate basic network statistics."""
        network = self.create_network()
        from .graph import centrality_measures, network_metrics

        metrics = network_metrics(network)
        centralities = centrality_measures(network)

        # Hub proteins (top 10% by degree centrality)
        degree_cents = centralities["degree"]
        sorted_proteins = sorted(degree_cents.items(), key=lambda x: x[1], reverse=True)
        n_hubs = max(1, len(sorted_proteins) // 10)
        hub_proteins = [protein for protein, _ in sorted_proteins[:n_hubs]]

        return {
            "total_proteins": len(self.proteins),
            "total_interactions": len(self.interactions),
            "network_metrics": metrics,
            "hub_proteins": hub_proteins,
            "avg_degree_centrality": np.mean(list(degree_cents.values())),
            "avg_betweenness_centrality": np.mean(list(centralities["betweenness"].values())),
        }


def load_string_interactions(
    string_file: str, score_threshold: int = 400, limit_organisms: Optional[List[str]] = None
) -> ProteinNetwork:
    """Load protein interactions from STRING database format.

    Args:
        string_file: Path to STRING interactions file
        score_threshold: Minimum combined score (0-1000)
        limit_organisms: List of organism taxonomy IDs to include

    Returns:
        ProteinNetwork object
    """
    ppi_network = ProteinNetwork(name="STRING_PPI")

    try:
        # Read STRING format: protein1 protein2 combined_score
        df = pd.read_csv(string_file, sep="\t", header=0)

        # Expected columns: protein1, protein2, combined_score
        required_cols = ["protein1", "protein2", "combined_score"]
        if not all(col in df.columns for col in required_cols):
            warnings.warn(f"Expected columns {required_cols}, found {list(df.columns)}")
            # Try alternative column names
            if "node1" in df.columns and "node2" in df.columns:
                df = df.rename(columns={"node1": "protein1", "node2": "protein2"})
            if "score" in df.columns and "combined_score" not in df.columns:
                df = df.rename(columns={"score": "combined_score"})

        # Filter by score threshold
        df_filtered = df[df["combined_score"] >= score_threshold]

        # Filter by organism if specified
        if limit_organisms:
            # Assume protein IDs contain organism info (e.g., "9606.ENSP00000...")
            organism_filter = df_filtered["protein1"].str.split(".").str[0].isin(limit_organisms)
            df_filtered = df_filtered[organism_filter]

        # Add interactions
        for _, row in df_filtered.iterrows():
            p1 = str(row["protein1"])
            p2 = str(row["protein2"])
            score = float(row["combined_score"])

            # Convert STRING score (0-1000) to confidence (0-1)
            confidence = score / 1000.0

            # Determine evidence types based on available columns
            evidence_types = []
            for col in ["experimental", "database", "textmining", "coexpression"]:
                if col in row and pd.notna(row[col]) and row[col] > 0:
                    evidence_types.append(col)

            ppi_network.add_interaction(
                p1, p2, confidence=confidence, evidence_types=evidence_types, combined_score=score
            )

    except Exception as e:
        warnings.warn(f"Error loading STRING file: {e}")

    return ppi_network


def predict_interactions(
    protein_features: np.ndarray, protein_ids: List[str], method: str = "correlation", threshold: float = 0.7
) -> ProteinNetwork:
    """Predict protein interactions from feature data.

    Args:
        protein_features: Feature matrix (proteins x features)
        protein_ids: List of protein identifiers
        method: Prediction method ("correlation", "coexpression")
        threshold: Minimum correlation/coexpression for interaction

    Returns:
        ProteinNetwork with predicted interactions
    """
    if len(protein_ids) != protein_features.shape[0]:
        raise ValueError("Number of protein IDs must match feature matrix rows")

    ppi_network = ProteinNetwork(name=f"Predicted_PPI_{method}")

    if method == "correlation":
        # Compute pairwise correlations
        corr_matrix = np.corrcoef(protein_features)

        for i in range(len(protein_ids)):
            for j in range(i + 1, len(protein_ids)):
                correlation = corr_matrix[i, j]

                if abs(correlation) >= threshold:
                    confidence = abs(correlation)
                    ppi_network.add_interaction(
                        protein_ids[i],
                        protein_ids[j],
                        confidence=confidence,
                        evidence_types=["predicted_correlation"],
                        correlation=correlation,
                    )

    elif method == "coexpression":
        # Simple coexpression based on feature similarity
        for i in range(len(protein_ids)):
            for j in range(i + 1, len(protein_ids)):
                # Cosine similarity
                vec1 = protein_features[i, :]
                vec2 = protein_features[j, :]

                dot_product = np.dot(vec1, vec2)
                norms = np.linalg.norm(vec1) * np.linalg.norm(vec2)

                if norms > 0:
                    similarity = dot_product / norms

                    if similarity >= threshold:
                        ppi_network.add_interaction(
                            protein_ids[i],
                            protein_ids[j],
                            confidence=similarity,
                            evidence_types=["predicted_coexpression"],
                            cosine_similarity=similarity,
                        )

    else:
        raise ValueError(f"Unknown prediction method: {method}")

    return ppi_network


def functional_enrichment_ppi(
    ppi_network: ProteinNetwork, functional_annotations: Dict[str, List[str]], min_confidence: float = 0.5
) -> Dict[str, Any]:
    """Analyze functional enrichment in PPI network.

    Args:
        ppi_network: Protein interaction network
        functional_annotations: Protein -> function annotations mapping
        min_confidence: Minimum interaction confidence

    Returns:
        Functional enrichment analysis results
    """
    # Create high-confidence network
    network = ppi_network.create_network(min_confidence=min_confidence)

    # Get network components/modules
    from .community import detect_communities

    communities = detect_communities(network, method="louvain")

    # Group proteins by community
    community_proteins = {}
    for protein, comm_id in communities.items():
        if comm_id not in community_proteins:
            community_proteins[comm_id] = []
        community_proteins[comm_id].append(protein)

    # Analyze functional enrichment in each community
    enrichment_results = {}

    for comm_id, proteins in community_proteins.items():
        if len(proteins) < 3:  # Skip small communities
            continue

        # Collect functions for proteins in this community
        community_functions = []
        for protein in proteins:
            if protein in functional_annotations:
                community_functions.extend(functional_annotations[protein])

        # Count function frequencies
        function_counts = {}
        for func in community_functions:
            function_counts[func] = function_counts.get(func, 0) + 1

        # Calculate enrichment (simplified)
        total_proteins = len(proteins)
        enriched_functions = {}

        for func, count in function_counts.items():
            frequency = count / total_proteins
            if frequency > 0.3:  # At least 30% of proteins have this function
                enriched_functions[func] = {
                    "count": count,
                    "frequency": frequency,
                    "proteins": [
                        p for p in proteins if p in functional_annotations and func in functional_annotations[p]
                    ],
                }

        if enriched_functions:
            enrichment_results[comm_id] = {"community_size": total_proteins, "enriched_functions": enriched_functions}

    return {
        "total_communities": len(community_proteins),
        "enriched_communities": len(enrichment_results),
        "enrichment_results": enrichment_results,
    }
