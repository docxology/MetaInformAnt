"""GWAS population structure analysis utilities.

This module provides functions for analyzing population structure in GWAS data,
including PCA computation and kinship matrix calculation.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Union, cast

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def _is_empty_genotype_matrix(genotype_matrix: Any) -> bool:
    """Return whether a list- or NumPy-backed genotype matrix has no cells."""
    if genotype_matrix is None:
        return True
    if hasattr(genotype_matrix, "size"):
        return int(genotype_matrix.size) == 0
    try:
        return len(genotype_matrix) == 0 or len(genotype_matrix[0]) == 0
    except TypeError, IndexError:
        return True


def compute_pca(genotype_matrix: List[List[int]], n_components: int = 10) -> Dict[str, Any]:
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
    if _is_empty_genotype_matrix(genotype_matrix):
        return {"status": "failed", "error": "Empty genotype matrix", "pcs": [], "explained_variance_ratio": []}

    try:
        # Convert to numpy array (impute missing -1 with mean)
        # genotype_matrix is n_samples x n_variants
        X = np.array(genotype_matrix, dtype=float)

        # Handle missing data (-1)
        # For simplicity in this optimization, replace < 0 with 0 (ref) or mean
        # A vectorized mean imputation is better
        mask = X < 0
        if np.any(mask):
            # Compute column means excluding -1
            col_means = np.sum(np.where(mask, 0, X), axis=0) / (np.sum(~mask, axis=0) + 1e-10)
            # Find indices where mask is True
            rows, cols = np.where(mask)
            X[rows, cols] = col_means[cols]

        # Center the data (subtract mean of each variant)
        X -= np.mean(X, axis=0)

        # Perform SVD: X = U S V^T
        # PC scores are X V = U S (if using full SVD definition)
        # but typical GWAS PCA projects samples.
        # sklearn PCA projects X onto components.
        # U, s, Vt = np.linalg.svd(X, full_matrices=False)
        # PCs = U * s

        U, s, Vt = np.linalg.svd(X, full_matrices=False)

        # Select top components
        n_comp = min(n_components, len(s))

        # PC scores (samples x components)
        # Scale U by singular values to get projected scores
        pcs_array = U[:, :n_comp] * s[:n_comp]
        pcs = pcs_array.tolist()

        # Explained variance ratio
        n_samples_svd = X.shape[0]
        eigenvalues = s**2 / (n_samples_svd - 1)
        total_variance = np.sum(eigenvalues)
        explained_variance_ratio = (eigenvalues[:n_comp] / total_variance).tolist()

        logger.info(f"PCA completed (SVD): {n_comp} components computed")

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
    if _is_empty_genotype_matrix(genotype_matrix):
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

        # Count missing data (genotype value -1)
        total_missing = sum(1 for sample_gts in genotype_matrix for g in sample_gts if g < 0)
        total_genotypes = n_samples * n_variants
        missing_rate = total_missing / total_genotypes if total_genotypes > 0 else 0.0

        return {
            "status": "success",
            "kinship_matrix": kinship,
            "method": method,
            "num_samples": n_samples,
            "missing_data_stats": {
                "total_missing": total_missing,
                "missing_rate": missing_rate,
            },
        }

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

    if _is_empty_genotype_matrix(genotype_matrix):
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
    """Compute kinship matrix using VanRaden (2008) Method 1.

    G = ZZ' / (2 * sum(p_i * (1-p_i)))
    where Z_ij = X_ij - 2*p_j and p_j is the allele frequency at locus j.

    Reference: VanRaden PM (2008) J Dairy Sci 91:4414-4423.

    Args:
        genotype_matrix: Genotype matrix (variants x samples)

    Returns:
        Kinship matrix
    """
    n_variants = len(genotype_matrix)
    n_samples = len(genotype_matrix[0])

    # Convert to numpy array
    X = np.array(genotype_matrix, dtype=float)

    # Handle missing data (< 0)
    # Compute allele frequencies on valid data
    # Creates a mask of valid entries
    mask = X >= 0

    # We need to compute p for each variant based on valid entries
    # p = sum(g) / (2 * count)
    # Avoid division by zero
    counts = 2.0 * np.sum(mask, axis=1)
    sums = np.sum(np.where(mask, X, 0), axis=1)

    # Initialize p with 0.5 (or fit) for variants with no data
    p = np.zeros(n_variants)
    valid_idx = counts > 0
    p[valid_idx] = sums[valid_idx] / counts[valid_idx]

    # Filter monomorphic
    valid_p_mask = (p > 0.0) & (p < 1.0)

    if not np.any(valid_p_mask):
        return np.eye(n_samples).tolist()

    # Filter genotype matrix to valid variants
    X = X[valid_p_mask, :]
    p = p[valid_p_mask]

    # Update n_variants
    n_variants = X.shape[0]

    # Impute missing values with 2*p (mean)
    # Z = X - 2p
    # If X is missing, Z = 0 (mean centered)
    # So we can replace missing X with 2p, then subtract 2p -> 0

    # Reshape p for broadcasting (n_variants, 1)
    p_col = p.reshape(-1, 1)

    # Replace missing with mean (2*p)
    # We need to iterate or use careful masking.
    # X is (variants, samples)
    # missing locations:
    missing_mask = ~mask[valid_p_mask, :]

    # Create Z matrix
    Z = X.copy()

    # Fill missing with 2*p
    # Create a matrix of means
    means = 2.0 * p_col
    # Broadcast standard numpy way if possible, or use where
    # Z[missing_mask] approach requires matching shapes which is tricky with boolean index
    # Easier: Z = np.where(missing, means, X)
    Z = np.where(missing_mask, means, Z)

    # Subtract mean (2*p)
    Z = Z - means

    # Compute scaling: 2 * sum(p * (1-p))
    scaling = float(2.0 * np.sum(p * (1.0 - p)))

    if scaling <= 0:
        return np.eye(n_samples).tolist()

    # Compute Kinship: K = Z.T @ Z / scaling
    # Z is (variants, samples). Z.T is (samples, variants)
    # Result is (samples, samples)
    K = np.matmul(Z.T, Z) / scaling

    return cast(List[List[float]], K.tolist())


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
    """Compute kinship matrix using the Astle-Balding standardization.

    Note: with the standardization ``(g - 2p) / sqrt(2p(1-p))`` and
    per-locus-pair normalization by the valid-locus count, this is numerically
    identical to the Yang GRM implemented in ``_yang_kinship``; both estimators
    are named in the literature, so both method names are offered. (A weighted
    Astle-Balding variant would normalize by the sum of ``2p(1-p)`` across loci
    instead.)

    Args:
        genotype_matrix: Genotype matrix (variants x samples)

    Returns:
        Kinship matrix
    """
    return _yang_kinship(genotype_matrix)


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
