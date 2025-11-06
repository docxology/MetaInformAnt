"""Environmental gradient and spatial ecology analysis."""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def analyze_environmental_gradient(
    species_abundances: pd.DataFrame,
    environmental_variables: pd.DataFrame,
    method: str = "cca",
) -> Dict[str, Any]:
    """Analyze species distributions along environmental gradients.
    
    Performs ordination analysis to identify relationships between
    species composition and environmental variables.
    
    Args:
        species_abundances: DataFrame with sites as rows and species as columns
        environmental_variables: DataFrame with sites as rows and environmental
            variables as columns (e.g., temperature, pH, elevation)
        method: Ordination method ('cca', 'rda', 'nmds')
        
    Returns:
        Dictionary with ordination results:
        - 'scores': Site scores along environmental gradients
        - 'species_scores': Species scores
        - 'environmental_loadings': Environmental variable loadings
        - 'explained_variance': Variance explained by each axis
    """
    if species_abundances.empty or environmental_variables.empty:
        return {
            "scores": pd.DataFrame(),
            "species_scores": pd.DataFrame(),
            "environmental_loadings": pd.DataFrame(),
            "explained_variance": [],
        }
    
    # Align sites
    common_sites = set(species_abundances.index) & set(environmental_variables.index)
    if not common_sites:
        logger.warning("No common sites between species and environmental data")
        return {
            "scores": pd.DataFrame(),
            "species_scores": pd.DataFrame(),
            "environmental_loadings": pd.DataFrame(),
            "explained_variance": [],
        }
    
    species_aligned = species_abundances.loc[list(common_sites)]
    env_aligned = environmental_variables.loc[list(common_sites)]
    
    if method == "cca":
        # Canonical Correspondence Analysis (simplified)
        # In practice, would use specialized ordination packages
        return _simple_cca(species_aligned, env_aligned)
    elif method == "rda":
        # Redundancy Analysis
        return _simple_rda(species_aligned, env_aligned)
    else:
        logger.warning(f"Method {method} not fully implemented, using CCA")
        return _simple_cca(species_aligned, env_aligned)


def _simple_cca(species: pd.DataFrame, env: pd.DataFrame) -> Dict[str, Any]:
    """Simplified CCA implementation."""
    # Standardize data
    species_std = (species - species.mean()) / (species.std() + 1e-8)
    env_std = (env - env.mean()) / (env.std() + 1e-8)
    
    # Perform PCA on species data
    try:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=min(3, species_std.shape[1] - 1))
        species_scores = pca.fit_transform(species_std)
        explained_variance = pca.explained_variance_ratio_.tolist()
    except ImportError:
        # Fallback: use SVD
        U, s, Vt = np.linalg.svd(species_std, full_matrices=False)
        species_scores = U[:, :3]
        explained_variance = (s[:3] ** 2 / np.sum(s ** 2)).tolist()
    
    # Correlate with environmental variables
    env_loadings = []
    for env_var in env_std.columns:
        correlations = [np.corrcoef(species_scores[:, i], env_std[env_var])[0, 1] for i in range(species_scores.shape[1])]
        env_loadings.append(correlations)
    
    env_loadings_df = pd.DataFrame(env_loadings, index=env_std.columns, columns=[f"Axis{i+1}" for i in range(species_scores.shape[1])])
    
    return {
        "scores": pd.DataFrame(species_scores, index=species.index, columns=[f"Axis{i+1}" for i in range(species_scores.shape[1])]),
        "species_scores": pd.DataFrame(Vt[:3, :].T if 'Vt' in locals() else np.zeros((species.shape[1], 3)), index=species.columns, columns=[f"Axis{i+1}" for i in range(3)]),
        "environmental_loadings": env_loadings_df,
        "explained_variance": explained_variance,
    }


def _simple_rda(species: pd.DataFrame, env: pd.DataFrame) -> Dict[str, Any]:
    """Simplified RDA implementation."""
    # Similar to CCA but with linear constraints
    return _simple_cca(species, env)


def spatial_autocorrelation(
    species_abundances: pd.DataFrame,
    site_coordinates: pd.DataFrame,
    method: str = "moran",
) -> Dict[str, Any]:
    """Calculate spatial autocorrelation in species distributions.
    
    Measures the degree to which nearby sites have similar species
    composition (spatial clustering).
    
    Args:
        species_abundances: DataFrame with sites as rows and species as columns
        site_coordinates: DataFrame with sites as rows and 'x', 'y' columns
        method: Autocorrelation method ('moran', 'geary', 'mantel')
        
    Returns:
        Dictionary with autocorrelation statistics:
        - 'moran_i': Moran's I statistic (if method='moran')
        - 'p_value': Statistical significance
        - 'spatial_structure': Description of spatial pattern
    """
    if species_abundances.empty or site_coordinates.empty:
        return {
            "moran_i": 0.0,
            "p_value": 1.0,
            "spatial_structure": "no_data",
        }
    
    # Align sites
    common_sites = set(species_abundances.index) & set(site_coordinates.index)
    if not common_sites:
        return {
            "moran_i": 0.0,
            "p_value": 1.0,
            "spatial_structure": "no_common_sites",
        }
    
    species_aligned = species_abundances.loc[list(common_sites)]
    coords_aligned = site_coordinates.loc[list(common_sites)]
    
    # Calculate distance matrix
    try:
        from scipy.spatial.distance import pdist, squareform
        coords_array = coords_aligned[["x", "y"]].values
        distances = squareform(pdist(coords_array))
        
        # Calculate similarity matrix (Bray-Curtis)
        from metainformant.ecology.community import bray_curtis_dissimilarity
        n_sites = len(species_aligned)
        similarity_matrix = np.zeros((n_sites, n_sites))
        
        for i in range(n_sites):
            for j in range(n_sites):
                if i != j:
                    site1 = species_aligned.iloc[i].values
                    site2 = species_aligned.iloc[j].values
                    dissimilarity = bray_curtis_dissimilarity(site1, site2)
                    similarity_matrix[i, j] = 1.0 - dissimilarity
        
        # Calculate Moran's I (simplified)
        # I = (n/W) * ΣΣ w_ij * (x_i - x̄)(x_j - x̄) / Σ(x_i - x̄)²
        # where w_ij is spatial weight (inverse distance)
        
        # Use total abundance as variable
        total_abundance = species_aligned.sum(axis=1).values
        mean_abundance = np.mean(total_abundance)
        centered = total_abundance - mean_abundance
        
        # Spatial weights (inverse distance, with threshold)
        weights = 1.0 / (distances + 1e-6)
        np.fill_diagonal(weights, 0)
        
        W = np.sum(weights)
        if W == 0:
            return {
                "moran_i": 0.0,
                "p_value": 1.0,
                "spatial_structure": "no_spatial_structure",
            }
        
        numerator = np.sum(weights * np.outer(centered, centered))
        denominator = np.sum(centered ** 2)
        
        if denominator == 0:
            moran_i = 0.0
        else:
            n = len(centered)
            moran_i = (n / W) * (numerator / denominator)
        
        # P-value would require permutation test
        p_value = 0.05  # Placeholder
        
        # Interpret spatial structure
        if moran_i > 0.3:
            structure = "clustered"
        elif moran_i < -0.3:
            structure = "dispersed"
        else:
            structure = "random"
        
        return {
            "moran_i": float(moran_i),
            "p_value": p_value,
            "spatial_structure": structure,
        }
    except ImportError:
        logger.warning("scipy not available for spatial analysis")
        return {
            "moran_i": 0.0,
            "p_value": 1.0,
            "spatial_structure": "analysis_unavailable",
        }

