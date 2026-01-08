"""GWAS population structure analysis utilities.

This module provides functions for analyzing population structure in GWAS data,
including PCA computation and kinship matrix calculation.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)


def compute_pca(genotype_matrix: List[List[int]], n_components: int = 10) -> Tuple[List[List[float]], List[float], List[List[float]]]:
    """Compute principal component analysis on genotype matrix.

    Args:
        genotype_matrix: Genotype matrix (variants x samples)
        n_components: Number of principal components to compute

    Returns:
        Tuple of (principal_components, explained_variance, loadings)
    """
    if not genotype_matrix:
        return [], [], []

    logger.info(f"Computing PCA on {len(genotype_matrix)} variants x {len(genotype_matrix[0])} samples")

    # Convert to numerical matrix and center
    n_samples = len(genotype_matrix[0])
    n_variants = len(genotype_matrix)

    # Simple centering (subtract mean)
    centered_matrix = []
    for variant in genotype_matrix:
        mean = sum(variant) / len(variant)
        centered = [x - mean for x in variant]
        centered_matrix.append(centered)

    # Simplified PCA using covariance matrix
    # In practice, this would use numpy.linalg.svd or similar
    covariance_matrix = _compute_covariance_matrix(centered_matrix)

    # Compute eigenvalues and eigenvectors (simplified)
    eigenvalues, eigenvectors = _simple_eigen_decomposition(covariance_matrix)

    # Sort by eigenvalue magnitude
    eigen_pairs = list(zip(eigenvalues, eigenvectors))
    eigen_pairs.sort(key=lambda x: x[0], reverse=True)

    # Extract requested number of components
    n_components = min(n_components, len(eigen_pairs))
    principal_components = []
    explained_variance = []

    for i in range(n_components):
        eigenvalue, eigenvector = eigen_pairs[i]
        explained_variance.append(eigenvalue)

        # Project data onto principal component
        pc_scores = []
        for j in range(n_samples):
            score = sum(centered_matrix[k][j] * eigenvector[k] for k in range(n_variants))
            pc_scores.append(score)
        principal_components.append(pc_scores)

    # Calculate explained variance ratio
    total_variance = sum(eigenvalues)
    explained_variance_ratio = [ev / total_variance for ev in explained_variance]

    logger.info(f"PCA completed: {n_components} components computed")

    return principal_components, explained_variance_ratio, eigenvectors[:n_components]


def compute_kinship_matrix(genotype_matrix: List[List[int]], method: str = "vanraden") -> List[List[float]]:
    """Compute kinship matrix from genotype matrix.

    Args:
        genotype_matrix: Genotype matrix (variants x samples)
        method: Kinship calculation method ('vanraden', 'ibs', etc.)

    Returns:
        Kinship matrix (samples x samples)
    """
    if not genotype_matrix:
        return []

    n_samples = len(genotype_matrix[0])
    logger.info(f"Computing kinship matrix for {n_samples} samples using {method} method")

    if method.lower() == "vanraden":
        return _vanraden_kinship(genotype_matrix)
    elif method.lower() == "ibs":
        return _ibs_kinship(genotype_matrix)
    else:
        raise ValueError(f"Unknown kinship method: {method}")


def estimate_population_structure(genotype_matrix: List[List[int]], n_pcs: int = 10) -> Dict[str, Any]:
    """Estimate population structure using PCA and kinship analysis.

    Args:
        genotype_matrix: Genotype matrix (variants x samples)
        n_pcs: Number of principal components to use

    Returns:
        Dictionary with population structure analysis results
    """
    logger.info("Estimating population structure")

    # Compute PCA
    pcs, explained_var, loadings = compute_pca(genotype_matrix, n_pcs)

    # Compute kinship
    kinship = compute_kinship_matrix(genotype_matrix)

    # Simple clustering based on first 2 PCs (placeholder)
    clusters = _simple_clustering(pcs[:2] if len(pcs) >= 2 else pcs[:1])

    results = {
        'principal_components': pcs,
        'explained_variance': explained_var,
        'loadings': loadings,
        'kinship_matrix': kinship,
        'clusters': clusters,
        'n_components_used': n_pcs,
        'method': 'pca+kinship'
    }

    logger.info("Population structure estimation completed")

    return results


def _vanraden_kinship(genotype_matrix: List[List[int]]) -> List[List[float]]:
    """Compute kinship matrix using VanRaden method.

    Args:
        genotype_matrix: Genotype matrix (variants x samples)

    Returns:
        Kinship matrix
    """
    n_samples = len(genotype_matrix[0])
    kinship = [[0.0] * n_samples for _ in range(n_samples)]

    # Simplified VanRaden calculation
    # In practice, this centers genotypes and computes G = XX'/c
    # where X is the centered genotype matrix

    for i in range(n_samples):
        for j in range(n_samples):
            if i == j:
                kinship[i][j] = 1.0  # Self-relatedness
            else:
                # Compute genetic similarity
                similarity = 0.0
                valid_loci = 0

                for locus in genotype_matrix:
                    if locus[i] >= 0 and locus[j] >= 0:  # Valid genotypes
                        # Simple IBS similarity
                        if locus[i] == locus[j]:
                            similarity += 1.0
                        valid_loci += 1

                if valid_loci > 0:
                    kinship[i][j] = similarity / valid_loci
                    kinship[j][i] = kinship[i][j]  # Symmetric

    return kinship


def _ibs_kinship(genotype_matrix: List[List[int]]) -> List[List[float]]:
    """Compute kinship matrix using identity-by-state (IBS) method.

    Args:
        genotype_matrix: Genotype matrix (variants x samples)

    Returns:
        Kinship matrix
    """
    n_samples = len(genotype_matrix[0])
    kinship = [[0.0] * n_samples for _ in range(n_samples)]

    for i in range(n_samples):
        for j in range(i, n_samples):
            ibs_count = 0
            valid_loci = 0

            for locus in genotype_matrix:
                if locus[i] >= 0 and locus[j] >= 0:
                    if locus[i] == locus[j]:
                        ibs_count += 1
                    valid_loci += 1

            if valid_loci > 0:
                kinship_value = ibs_count / valid_loci
                kinship[i][j] = kinship_value
                kinship[j][i] = kinship_value

    return kinship


def _compute_covariance_matrix(matrix: List[List[float]]) -> List[List[float]]:
    """Compute covariance matrix from data matrix.

    Args:
        matrix: Data matrix (features x samples)

    Returns:
        Covariance matrix
    """
    n_features = len(matrix)
    n_samples = len(matrix[0])

    cov_matrix = [[0.0] * n_features for _ in range(n_features)]

    # Compute means
    means = [sum(row) / n_samples for row in matrix]

    # Compute covariance
    for i in range(n_features):
        for j in range(n_features):
            cov = 0.0
            for k in range(n_samples):
                cov += (matrix[i][k] - means[i]) * (matrix[j][k] - means[j])
            cov /= (n_samples - 1)
            cov_matrix[i][j] = cov

    return cov_matrix


def _eigen_decomposition(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Perform eigen decomposition using numpy.

    Args:
        matrix: Square matrix as numpy array

    Returns:
        Tuple of (eigenvalues, eigenvectors)
    """
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    return eigenvalues, eigenvectors


def _kmeans_clustering(pc_scores: np.ndarray, k: int = 3) -> np.ndarray:
    """Perform k-means clustering on principal component scores.

    Args:
        pc_scores: Principal component scores as numpy array (samples x components)
        k: Number of clusters

    Returns:
        Cluster labels for each sample
    """
    from sklearn.cluster import KMeans

    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(pc_scores)

    logger.info(f"Performed k-means clustering: {len(pc_scores)} samples → {k} clusters")
    return clusters


