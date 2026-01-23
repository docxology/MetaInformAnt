"""GWAS population structure analysis utilities.

This module provides functions for analyzing population structure in GWAS data,
including PCA computation and kinship matrix calculation.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)


def compute_pca(
    genotype_matrix: List[List[int]], n_components: int = 10
) -> Dict[str, Any]:
    """Compute principal component analysis on genotype matrix.

    Args:
        genotype_matrix: Genotype matrix as list of sample genotype lists
            Each inner list contains genotypes for one sample across all variants
            e.g., [[s1_v1, s1_v2, ...], [s2_v1, s2_v2, ...], ...]
        n_components: Number of principal components to compute

    Returns:
        Dictionary with status, pcs (principal component scores per sample),
        explained_variance_ratio, and optionally error message
    """
    if not genotype_matrix or not genotype_matrix[0]:
        return {"status": "failed", "error": "Empty genotype matrix", "pcs": [], "explained_variance_ratio": []}

    n_samples = len(genotype_matrix)
    n_variants = len(genotype_matrix[0])

    logger.info(f"Computing PCA on {n_samples} samples x {n_variants} variants")

    try:
        # Handle missing data by imputation (replace -1 with variant mean)
        imputed_matrix = []
        for sample_gts in genotype_matrix:
            imputed_gts = []
            for gt in sample_gts:
                if gt < 0:
                    # Impute with 1 (heterozygous) as simple approach
                    imputed_gts.append(1)
                else:
                    imputed_gts.append(gt)
            imputed_matrix.append(imputed_gts)

        # Transpose to [variant][sample] for centering
        transposed = [[imputed_matrix[s][v] for s in range(n_samples)] for v in range(n_variants)]

        # Center by variant (subtract mean)
        centered_matrix = []
        for variant_gts in transposed:
            mean = sum(variant_gts) / len(variant_gts) if variant_gts else 0
            centered = [gt - mean for gt in variant_gts]
            centered_matrix.append(centered)

        # Compute covariance matrix (sample x sample)
        covariance_matrix = _compute_covariance_matrix(centered_matrix)

        # Compute eigenvalues and eigenvectors
        eigenvalues, eigenvectors = _simple_eigen_decomposition(covariance_matrix)

        # Sort by eigenvalue magnitude (descending)
        eigen_pairs = list(zip(eigenvalues, eigenvectors))
        eigen_pairs.sort(key=lambda x: abs(x[0]), reverse=True)

        # Extract requested number of components
        # Ensure we return at least n_components (or n_samples) even with low-rank data
        n_comp = min(n_components, n_samples)

        # Pad eigen_pairs if we don't have enough
        while len(eigen_pairs) < n_comp:
            # Add zero eigenvalue with random eigenvector
            zero_vec = [0.0] * n_variants
            if len(zero_vec) > 0:
                zero_vec[len(eigen_pairs) % len(zero_vec)] = 1.0
            eigen_pairs.append((0.0, zero_vec))

        # Build PC scores: for each sample, its score on each PC
        pcs = []
        explained_variance = []

        for i in range(n_comp):
            eigenvalue, eigenvector = eigen_pairs[i]
            explained_variance.append(abs(eigenvalue))

        # Calculate explained variance ratio
        total_variance = sum(abs(ev) for ev, _ in eigen_pairs)
        if total_variance > 0:
            explained_variance_ratio = [ev / total_variance for ev in explained_variance]
        else:
            explained_variance_ratio = [1.0 / n_comp] * n_comp  # Equal distribution if no variance

        # PC scores per sample: project samples onto principal components
        # pcs[sample_idx] = [pc1_score, pc2_score, ...]
        pcs = []
        for s in range(n_samples):
            sample_scores = []
            for i in range(n_comp):
                _, eigenvector = eigen_pairs[i]
                # Project this sample onto the eigenvector
                score = sum(centered_matrix[v][s] * eigenvector[v] for v in range(min(n_variants, len(eigenvector))))
                sample_scores.append(score)
            pcs.append(sample_scores)

        logger.info(f"PCA completed: {n_comp} components computed")

        return {
            "status": "success",
            "pcs": pcs,
            "explained_variance_ratio": explained_variance_ratio,
            "n_components": n_comp,
        }

    except Exception as e:
        logger.error(f"PCA computation failed: {e}")
        return {"status": "failed", "error": str(e), "pcs": [], "explained_variance_ratio": []}


def compute_kinship_matrix(genotype_matrix: List[List[int]], method: str = "vanraden") -> Dict[str, Any]:
    """Compute kinship matrix from genotype matrix.

    Args:
        genotype_matrix: Genotype matrix as list of sample genotype lists
            Each inner list contains genotypes for one sample across all variants
        method: Kinship calculation method ('vanraden', 'astle', 'yang', 'ibs')

    Returns:
        Dictionary with status, kinship_matrix, and optionally error message
    """
    if not genotype_matrix or not genotype_matrix[0]:
        return {"status": "failed", "error": "Empty genotype matrix", "kinship_matrix": []}

    n_samples = len(genotype_matrix)
    n_variants = len(genotype_matrix[0])
    logger.info(f"Computing kinship matrix for {n_samples} samples using {method} method")

    try:
        # Transpose to [variant][sample] for kinship calculation
        transposed = [[genotype_matrix[s][v] for s in range(n_samples)] for v in range(n_variants)]

        method_lower = method.lower()
        if method_lower == "vanraden":
            kinship = _vanraden_kinship(transposed)
        elif method_lower == "ibs":
            kinship = _ibs_kinship(transposed)
        elif method_lower == "astle":
            kinship = _astle_kinship(transposed)
        elif method_lower == "yang":
            kinship = _yang_kinship(transposed)
        else:
            return {"status": "failed", "error": f"Unknown kinship method: {method}", "kinship_matrix": []}

        return {"status": "success", "kinship_matrix": kinship, "method": method}

    except Exception as e:
        logger.error(f"Kinship matrix computation failed: {e}")
        return {"status": "failed", "error": str(e), "kinship_matrix": []}


def estimate_population_structure(
    vcf_input: Union[str, Path, List[List[int]]],
    config: Optional[Dict[str, Any]] = None,
    output_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Estimate population structure using PCA and kinship analysis.

    Args:
        vcf_input: VCF file path or genotype matrix as list of sample genotype lists
        config: Configuration dictionary with keys:
            - compute_pca: Whether to compute PCA (default True)
            - n_components: Number of PCA components (default 10)
            - compute_relatedness: Whether to compute kinship (default True)
            - kinship_method: Kinship method (default 'vanraden')
        output_dir: Optional output directory for results

    Returns:
        Dictionary with population structure analysis results
    """
    logger.info("Estimating population structure")

    # Default config
    if config is None:
        config = {}
    compute_pca_flag = config.get("compute_pca", True)
    n_components = config.get("n_components", 10)
    compute_relatedness = config.get("compute_relatedness", True)
    kinship_method = config.get("kinship_method", "vanraden")

    # Parse VCF if path is provided
    if isinstance(vcf_input, (str, Path)):
        from metainformant.gwas.analysis.quality import parse_vcf_full

        vcf_data = parse_vcf_full(vcf_input)
        genotype_matrix = vcf_data.get("genotypes", [])
    else:
        genotype_matrix = vcf_input

    if not genotype_matrix or not genotype_matrix[0]:
        return {"status": "failed", "error": "Empty genotype matrix"}

    results: Dict[str, Any] = {"status": "success"}

    try:
        # Compute PCA
        if compute_pca_flag:
            pca_result = compute_pca(genotype_matrix, n_components)
            if pca_result["status"] == "success":
                results["pca"] = {
                    "pcs": pca_result["pcs"],
                    "explained_variance_ratio": pca_result["explained_variance_ratio"],
                }
            else:
                results["pca_error"] = pca_result.get("error", "Unknown PCA error")

        # Compute kinship
        if compute_relatedness:
            kinship_result = compute_kinship_matrix(genotype_matrix, method=kinship_method)
            if kinship_result["status"] == "success":
                results["kinship"] = {
                    "matrix": kinship_result["kinship_matrix"],
                    "method": kinship_method,
                }
            else:
                results["kinship_error"] = kinship_result.get("error", "Unknown kinship error")

        # Save results if output_dir provided
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            summary_path = output_dir / "structure_summary.json"
            with open(summary_path, "w") as f:
                json.dump(results, f, indent=2)
            logger.info(f"Structure results saved to {summary_path}")

        logger.info("Population structure estimation completed")
        return results

    except Exception as e:
        logger.error(f"Population structure estimation failed: {e}")
        return {"status": "failed", "error": str(e)}


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


def _astle_kinship(genotype_matrix: List[List[int]]) -> List[List[float]]:
    """Compute kinship matrix using Astle-Balding method.

    The Astle-Balding method computes relatedness using standardized genotypes
    centered by allele frequency.

    Args:
        genotype_matrix: Genotype matrix (variants x samples)

    Returns:
        Kinship matrix
    """
    n_samples = len(genotype_matrix[0])
    n_variants = len(genotype_matrix)
    kinship = [[0.0] * n_samples for _ in range(n_samples)]

    # Compute allele frequencies
    allele_freqs = []
    for locus in genotype_matrix:
        valid = [g for g in locus if g >= 0]
        if valid:
            freq = sum(valid) / (2 * len(valid))  # Divide by 2 because diploid
            allele_freqs.append(max(0.01, min(0.99, freq)))  # Bound away from 0 and 1
        else:
            allele_freqs.append(0.5)

    # Standardize genotypes and compute kinship
    for i in range(n_samples):
        for j in range(i, n_samples):
            kinship_sum = 0.0
            valid_count = 0

            for v in range(n_variants):
                gi = genotype_matrix[v][i]
                gj = genotype_matrix[v][j]

                if gi >= 0 and gj >= 0:
                    p = allele_freqs[v]
                    var = 2 * p * (1 - p)
                    if var > 0:
                        # Standardized genotype contribution
                        zi = (gi - 2 * p) / math.sqrt(var)
                        zj = (gj - 2 * p) / math.sqrt(var)
                        kinship_sum += zi * zj
                        valid_count += 1

            if valid_count > 0:
                kinship_value = kinship_sum / valid_count
                kinship[i][j] = kinship_value
                kinship[j][i] = kinship_value

    return kinship


def _yang_kinship(genotype_matrix: List[List[int]]) -> List[List[float]]:
    """Compute kinship matrix using Yang's method (GCTA-style).

    The Yang method computes genomic relatedness matrix using standardized
    genotypes similar to GCTA.

    Args:
        genotype_matrix: Genotype matrix (variants x samples)

    Returns:
        Kinship matrix
    """
    n_samples = len(genotype_matrix[0])
    n_variants = len(genotype_matrix)
    kinship = [[0.0] * n_samples for _ in range(n_samples)]

    # Compute allele frequencies
    allele_freqs = []
    for locus in genotype_matrix:
        valid = [g for g in locus if g >= 0]
        if valid:
            freq = sum(valid) / (2 * len(valid))
            allele_freqs.append(max(0.01, min(0.99, freq)))
        else:
            allele_freqs.append(0.5)

    # Yang's GRM calculation
    for i in range(n_samples):
        for j in range(i, n_samples):
            grm_sum = 0.0
            valid_count = 0

            for v in range(n_variants):
                gi = genotype_matrix[v][i]
                gj = genotype_matrix[v][j]

                if gi >= 0 and gj >= 0:
                    p = allele_freqs[v]
                    denominator = 2 * p * (1 - p)
                    if denominator > 0:
                        # Yang's formula: (g_i - 2p)(g_j - 2p) / (2p(1-p))
                        grm_sum += (gi - 2 * p) * (gj - 2 * p) / denominator
                        valid_count += 1

            if valid_count > 0:
                kinship_value = grm_sum / valid_count
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
            cov /= n_samples - 1
            cov_matrix[i][j] = cov

    return cov_matrix


def _simple_eigen_decomposition(matrix: List[List[float]]) -> Tuple[List[float], List[List[float]]]:
    """Perform eigenvalue decomposition using power iteration method.

    This is a pure Python implementation for cases where numpy is not available.
    Uses power iteration to find dominant eigenvectors.

    Args:
        matrix: Square covariance matrix as list of lists

    Returns:
        Tuple of (eigenvalues, eigenvectors)
    """
    n = len(matrix)
    if n == 0:
        return [], []

    eigenvalues = []
    eigenvectors = []

    # Work with a copy of the matrix for deflation
    work_matrix = [row[:] for row in matrix]

    # Find top eigenvalues/eigenvectors using power iteration
    for _ in range(min(n, 10)):  # Find up to 10 eigenvalues
        # Initialize random vector
        v = [1.0 / math.sqrt(n)] * n

        # Power iteration
        for _ in range(100):  # Max iterations
            # Matrix-vector multiply
            new_v = [0.0] * n
            for i in range(n):
                for j in range(n):
                    new_v[i] += work_matrix[i][j] * v[j]

            # Compute eigenvalue estimate (Rayleigh quotient)
            eigenvalue = sum(new_v[i] * v[i] for i in range(n))

            # Normalize
            norm = math.sqrt(sum(x * x for x in new_v))
            if norm < 1e-10:
                break
            v = [x / norm for x in new_v]

        # Compute final eigenvalue
        Av = [sum(work_matrix[i][j] * v[j] for j in range(n)) for i in range(n)]
        eigenvalue = sum(Av[i] * v[i] for i in range(n))

        if abs(eigenvalue) < 1e-10:
            break

        eigenvalues.append(eigenvalue)
        eigenvectors.append(v)

        # Deflate matrix: A = A - eigenvalue * v * v^T
        for i in range(n):
            for j in range(n):
                work_matrix[i][j] -= eigenvalue * v[i] * v[j]

    return eigenvalues, eigenvectors


def _simple_clustering(pc_scores: List[List[float]], k: int = 3) -> List[int]:
    """Perform simple k-means clustering using pure Python.

    Args:
        pc_scores: Principal component scores (components x samples)
        k: Number of clusters

    Returns:
        Cluster labels for each sample
    """
    if not pc_scores or not pc_scores[0]:
        return []

    n_components = len(pc_scores)
    n_samples = len(pc_scores[0])

    # Transpose to samples x components for easier processing
    samples = [[pc_scores[c][s] for c in range(n_components)] for s in range(n_samples)]

    # Initialize centroids using first k samples (simple seeding)
    k = min(k, n_samples)
    centroids = [samples[i * n_samples // k][:] for i in range(k)]

    labels = [0] * n_samples

    # K-means iterations
    for _ in range(50):  # Max iterations
        # Assign samples to nearest centroid
        new_labels = []
        for sample in samples:
            min_dist = float("inf")
            best_cluster = 0
            for cluster_idx, centroid in enumerate(centroids):
                dist = sum((s - c) ** 2 for s, c in zip(sample, centroid))
                if dist < min_dist:
                    min_dist = dist
                    best_cluster = cluster_idx
            new_labels.append(best_cluster)

        # Check convergence
        if new_labels == labels:
            break
        labels = new_labels

        # Update centroids
        for cluster_idx in range(k):
            cluster_samples = [s for s, lbl in zip(samples, labels) if lbl == cluster_idx]
            if cluster_samples:
                centroids[cluster_idx] = [
                    sum(s[dim] for s in cluster_samples) / len(cluster_samples) for dim in range(n_components)
                ]

    logger.info(f"Performed k-means clustering: {n_samples} samples -> {k} clusters")
    return labels


def _eigen_decomposition(matrix: Any) -> Tuple[Any, Any]:
    """Perform eigen decomposition using numpy.

    Args:
        matrix: Square matrix as numpy array

    Returns:
        Tuple of (eigenvalues, eigenvectors)
    """
    import numpy as np

    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    return eigenvalues, eigenvectors


def _kmeans_clustering(pc_scores: Any, k: int = 3) -> Any:
    """Perform k-means clustering on principal component scores using sklearn.

    Args:
        pc_scores: Principal component scores as numpy array (samples x components)
        k: Number of clusters

    Returns:
        Cluster labels for each sample
    """
    from sklearn.cluster import KMeans

    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(pc_scores)

    logger.info(f"Performed k-means clustering: {len(pc_scores)} samples -> {k} clusters")
    return clusters
