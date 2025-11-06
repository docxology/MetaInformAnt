"""Species interaction network analysis for ecological communities."""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def build_interaction_network(
    interaction_matrix: pd.DataFrame,
    interaction_type: str = "cooccurrence",
    threshold: float = 0.5,
) -> Dict[str, Any]:
    """Build species interaction network from interaction matrix.
    
    Creates a network representation of species interactions (competition,
    mutualism, predation, etc.) from pairwise interaction strengths.
    
    Args:
        interaction_matrix: DataFrame with species as rows and columns,
            values represent interaction strengths
        interaction_type: Type of interactions ('cooccurrence', 'competition', 'mutualism')
        threshold: Minimum interaction strength to include edge
        
    Returns:
        Dictionary with network representation:
        - 'nodes': List of species identifiers
        - 'edges': List of (species1, species2, weight) tuples
        - 'adjacency_matrix': Adjacency matrix representation
        - 'n_nodes': Number of species
        - 'n_edges': Number of interactions
    """
    if interaction_matrix.empty:
        return {
            "nodes": [],
            "edges": [],
            "adjacency_matrix": np.array([]),
            "n_nodes": 0,
            "n_edges": 0,
        }
    
    # Ensure symmetric matrix
    if not interaction_matrix.index.equals(interaction_matrix.columns):
        # Align indices
        common_species = set(interaction_matrix.index) & set(interaction_matrix.columns)
        interaction_matrix = interaction_matrix.loc[list(common_species), list(common_species)]
    
    # Filter by threshold
    adjacency = (interaction_matrix.abs() >= threshold).astype(int)
    
    # Remove self-loops
    np.fill_diagonal(adjacency.values, 0)
    
    # Extract edges
    edges = []
    species_list = list(interaction_matrix.index)
    
    for i, species1 in enumerate(species_list):
        for j, species2 in enumerate(species_list):
            if i < j and adjacency.iloc[i, j] > 0:
                weight = interaction_matrix.iloc[i, j]
                edges.append((species1, species2, float(weight)))
    
    return {
        "nodes": species_list,
        "edges": edges,
        "adjacency_matrix": adjacency.values,
        "n_nodes": len(species_list),
        "n_edges": len(edges),
    }


def calculate_interaction_strength(
    abundance1: List[float],
    abundance2: List[float],
    method: str = "correlation",
) -> float:
    """Calculate interaction strength between two species.
    
    Measures the strength of interaction (positive or negative) between
    two species based on their abundance patterns across sites.
    
    Args:
        abundance1: Abundance values for species 1 across sites
        abundance2: Abundance values for species 2 across sites
        method: Method for calculating strength ('correlation', 'mutual_information')
        
    Returns:
        Interaction strength value (positive = positive interaction, negative = negative interaction)
    """
    if len(abundance1) != len(abundance2):
        raise ValueError("Abundance vectors must have same length")
    
    if method == "correlation":
        try:
            correlation = np.corrcoef(abundance1, abundance2)[0, 1]
            return float(correlation) if not np.isnan(correlation) else 0.0
        except Exception:
            return 0.0
    elif method == "mutual_information":
        try:
            from metainformant.information import mutual_information
            mi = mutual_information(abundance1, abundance2)
            return float(mi)
        except ImportError:
            logger.warning("Information theory module not available, using correlation")
            return calculate_interaction_strength(abundance1, abundance2, method="correlation")
    else:
        raise ValueError(f"Unknown method: {method}")


def identify_keystone_species(
    interaction_network: Dict[str, Any],
    method: str = "betweenness",
) -> List[Tuple[str, float]]:
    """Identify keystone species in interaction network.
    
    Keystone species have disproportionate impact on community structure.
    Identified using network centrality measures.
    
    Args:
        interaction_network: Network dictionary from build_interaction_network
        method: Centrality measure ('betweenness', 'degree', 'closeness')
        
    Returns:
        List of (species_id, centrality_score) tuples, sorted by score (descending)
    """
    if interaction_network["n_nodes"] == 0:
        return []
    
    try:
        import networkx as nx
    except ImportError:
        logger.warning("networkx not available, using simple degree centrality")
        method = "degree"
    
    if method == "degree":
        # Simple degree centrality
        adjacency = interaction_network["adjacency_matrix"]
        degrees = np.sum(adjacency, axis=1)
        
        species_scores = [
            (species, float(degree))
            for species, degree in zip(interaction_network["nodes"], degrees)
        ]
    else:
        # Use networkx for advanced measures
        try:
            G = nx.from_numpy_array(interaction_network["adjacency_matrix"])
            nx.relabel_nodes(G, {i: interaction_network["nodes"][i] for i in range(len(interaction_network["nodes"]))}, copy=False)
            
            if method == "betweenness":
                centrality = nx.betweenness_centrality(G)
            elif method == "closeness":
                centrality = nx.closeness_centrality(G)
            else:
                centrality = nx.degree_centrality(G)
            
            species_scores = [(species, score) for species, score in centrality.items()]
        except Exception as e:
            logger.warning(f"Network analysis failed: {e}, using degree")
            return identify_keystone_species(interaction_network, method="degree")
    
    # Sort by score (descending)
    species_scores.sort(key=lambda x: x[1], reverse=True)
    
    return species_scores

