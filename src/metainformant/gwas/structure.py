"""Population structure analysis (PCA and kinship)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from ..core.io import dump_json, ensure_directory, write_delimited
from ..core.logging import get_logger
from .quality import parse_vcf_full

logger = get_logger(__name__)

# Try importing scipy for advanced matrix operations
try:
    from scipy.linalg import eigh

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    eigh = None


def compute_pca(
    genotype_matrix: list[list[int]],
    n_components: int = 10,
    *,
    missing_imputation: str = "mean",
    max_missing_per_variant: float = 0.5,
) -> dict[str, Any]:
    """Compute Principal Component Analysis on genotype matrix.

    Args:
        genotype_matrix: Genotype matrix (samples x variants), encoded as 0/1/2/-1
        n_components: Number of principal components to compute
        missing_imputation: Method for handling missing data. Options:
            - "mean": Mean imputation per variant (default)
            - "median": Median imputation per variant
            - "mode": Mode imputation per variant (most common genotype)
            - "skip": Skip variants with missing data (set to 0)
        max_missing_per_variant: Maximum fraction of missing data per variant.
            Variants with higher missingness will be excluded. Default 0.5 (50%).

    Returns:
        Dictionary with:
        - pcs: Principal components (samples x n_components)
        - explained_variance: Explained variance per component
        - explained_variance_ratio: Explained variance ratio per component
        - missing_data_stats: Statistics about missing data handling
    """
    logger.info(f"compute_pca: Computing {n_components} principal components")

    if not genotype_matrix or not genotype_matrix[0]:
        return {
            "status": "failed",
            "error": "Empty genotype matrix",
        }

    # Convert to numpy array
    try:
        X = np.array(genotype_matrix, dtype=float)
    except Exception as exc:
        return {
            "status": "failed",
            "error": f"Failed to convert to numpy array: {exc}",
        }

    num_samples, num_variants = X.shape
    logger.info(f"compute_pca: Matrix shape: {num_samples} samples x {num_variants} variants")

    # Track missing data statistics
    total_missing = np.sum(X == -1)
    total_cells = num_samples * num_variants
    overall_missing_rate = total_missing / total_cells if total_cells > 0 else 0.0
    
    # Count missing per variant and filter high-missingness variants
    missing_per_variant = np.sum(X == -1, axis=0)
    missing_rate_per_variant = missing_per_variant / num_samples
    high_missing_variants = np.sum(missing_rate_per_variant > max_missing_per_variant)
    
    # Initialize exclude_mask (all False if no high-missing variants)
    if high_missing_variants > 0:
        logger.warning(
            f"compute_pca: Excluding {high_missing_variants} variants with >{max_missing_per_variant:.1%} missing data"
        )
        # Set high-missing variants to 0 (will be excluded in analysis)
        exclude_mask = missing_rate_per_variant > max_missing_per_variant
        X[:, exclude_mask] = 0.0
    else:
        exclude_mask = np.zeros(num_variants, dtype=bool)
    
    # Handle missing data: imputation per variant
    imputed_count = 0
    for var_idx in range(num_variants):
        var_data = X[:, var_idx]
        missing_mask = var_data == -1
        
        if missing_mask.any() and not exclude_mask[var_idx]:
            imputed_count += np.sum(missing_mask)
            
            if missing_imputation == "mean":
                impute_val = np.mean(var_data[~missing_mask]) if np.any(~missing_mask) else 0.0
            elif missing_imputation == "median":
                impute_val = np.median(var_data[~missing_mask]) if np.any(~missing_mask) else 0.0
            elif missing_imputation == "mode":
                # Mode: most common non-missing genotype
                non_missing = var_data[~missing_mask]
                if len(non_missing) > 0:
                    values, counts = np.unique(non_missing, return_counts=True)
                    impute_val = values[np.argmax(counts)]
                else:
                    impute_val = 0.0
            elif missing_imputation == "skip":
                impute_val = 0.0  # Set to 0 for skipped variants
            else:
                logger.warning(f"compute_pca: Unknown imputation method '{missing_imputation}', using mean")
                impute_val = np.mean(var_data[~missing_mask]) if np.any(~missing_mask) else 0.0
            
            X[missing_mask, var_idx] = impute_val
    
    # Log missing data statistics
    if total_missing > 0:
        logger.info(
            f"compute_pca: Missing data: {total_missing}/{total_cells} ({overall_missing_rate:.2%}) "
            f"cells, imputed {imputed_count} values using {missing_imputation} method"
        )
    
    missing_stats = {
        "total_missing": int(total_missing),
        "total_cells": int(total_cells),
        "overall_missing_rate": float(overall_missing_rate),
        "variants_excluded": int(high_missing_variants),
        "values_imputed": int(imputed_count),
        "imputation_method": missing_imputation,
    }

    # Center the data (subtract mean per variant)
    X_centered = X - np.mean(X, axis=0, keepdims=True)

    # Compute covariance matrix (samples x samples)
    # Note: For GWAS, we usually compute variants x variants covariance and get PCs of samples
    # But with large variant counts, it's more efficient to compute sample x sample covariance
    try:
        # Standardize (optional but often done)
        X_std = X_centered / (np.std(X_centered, axis=0, keepdims=True) + 1e-10)

        # Compute sample covariance matrix (samples x samples)
        # Cov = (1/n) * X^T * X where X is (samples x variants), so Cov is (samples x samples)
        cov_matrix = np.dot(X_std, X_std.T) / (num_variants - 1)

        # Compute eigenvalues and eigenvectors
        if SCIPY_AVAILABLE and eigh:
            eigenvalues, eigenvectors = eigh(cov_matrix)
        else:
            # Fallback to numpy
            logger.warning(
                "compute_pca: scipy not available, using numpy.linalg.eigh. "
                "Install scipy for better performance and numerical stability."
            )
            eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)

        # Sort by eigenvalue (descending)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        # Select top n_components
        n_components_actual = min(n_components, len(eigenvalues))
        pcs = eigenvectors[:, :n_components_actual]
        explained_variance = eigenvalues[:n_components_actual]
        explained_variance_ratio = explained_variance / np.sum(eigenvalues)

        logger.info(
            f"compute_pca: Computed {n_components_actual} components, "
            f"explained variance: {np.sum(explained_variance_ratio) * 100:.2f}%"
        )

        return {
            "status": "success",
            "pcs": pcs.tolist(),
            "explained_variance": explained_variance.tolist(),
            "explained_variance_ratio": explained_variance_ratio.tolist(),
            "n_components": n_components_actual,
            "missing_data_stats": missing_stats,
        }

    except Exception as exc:
        logger.error(f"compute_pca: Error computing PCA: {exc}")
        return {
            "status": "failed",
            "error": str(exc),
        }


def compute_kinship_matrix(
    genotype_matrix: list[list[int]],
    method: str = "vanraden",
    *,
    missing_imputation: str = "mean",
    max_missing_per_variant: float = 0.5,
) -> dict[str, Any]:
    """Compute kinship matrix from genotype data.

    Supported methods:
    - vanraden: VanRaden method (default) - K = (M * M.T) / (2 * sum(p * (1-p)))
      where M is centered genotypes
    - astle: Astle-Balding method
    - yang: Yang et al. method

    Args:
        genotype_matrix: Genotype matrix (samples x variants), encoded as 0/1/2/-1
        method: Kinship computation method
        missing_imputation: Method for handling missing data. Options:
            - "mean": Mean imputation per variant (default)
            - "median": Median imputation per variant
            - "mode": Mode imputation per variant (most common genotype)
            - "skip": Skip variants with missing data (set to 0)
        max_missing_per_variant: Maximum fraction of missing data per variant.
            Variants with higher missingness will be excluded. Default 0.5 (50%).

    Returns:
        Dictionary with:
        - kinship_matrix: Kinship matrix (samples x samples)
        - method: Method used
        - missing_data_stats: Statistics about missing data handling
    """
    logger.info(f"compute_kinship_matrix: Computing kinship using {method} method")

    if not genotype_matrix or not genotype_matrix[0]:
        return {
            "status": "failed",
            "error": "Empty genotype matrix",
        }

    try:
        X = np.array(genotype_matrix, dtype=float)
        num_samples, num_variants = X.shape

        # Track missing data statistics
        total_missing = np.sum(X == -1)
        total_cells = num_samples * num_variants
        overall_missing_rate = total_missing / total_cells if total_cells > 0 else 0.0
        
        # Count missing per variant and filter high-missingness variants
        missing_per_variant = np.sum(X == -1, axis=0)
        missing_rate_per_variant = missing_per_variant / num_samples
        high_missing_variants = np.sum(missing_rate_per_variant > max_missing_per_variant)
        
        if high_missing_variants > 0:
            logger.warning(
                f"compute_kinship_matrix: Excluding {high_missing_variants} variants with >{max_missing_per_variant:.1%} missing data"
            )
            exclude_mask = missing_rate_per_variant > max_missing_per_variant
            X[:, exclude_mask] = 0.0
        else:
            exclude_mask = np.zeros(num_variants, dtype=bool)

        # Handle missing data: imputation per variant
        imputed_count = 0
        for var_idx in range(num_variants):
            var_data = X[:, var_idx]
            missing_mask = var_data == -1
            
            if missing_mask.any() and not exclude_mask[var_idx]:
                imputed_count += np.sum(missing_mask)
                
                if missing_imputation == "mean":
                    impute_val = np.mean(var_data[~missing_mask]) if np.any(~missing_mask) else 0.0
                elif missing_imputation == "median":
                    impute_val = np.median(var_data[~missing_mask]) if np.any(~missing_mask) else 0.0
                elif missing_imputation == "mode":
                    # Mode: most common non-missing genotype
                    non_missing = var_data[~missing_mask]
                    if len(non_missing) > 0:
                        values, counts = np.unique(non_missing, return_counts=True)
                        impute_val = values[np.argmax(counts)]
                    else:
                        impute_val = 0.0
                elif missing_imputation == "skip":
                    impute_val = 0.0
                else:
                    logger.warning(f"compute_kinship_matrix: Unknown imputation method '{missing_imputation}', using mean")
                    impute_val = np.mean(var_data[~missing_mask]) if np.any(~missing_mask) else 0.0
                
                X[missing_mask, var_idx] = impute_val
        
        # Log missing data statistics
        if total_missing > 0:
            logger.info(
                f"compute_kinship_matrix: Missing data: {total_missing}/{total_cells} ({overall_missing_rate:.2%}) "
                f"cells, imputed {imputed_count} values using {missing_imputation} method"
            )
        
        missing_stats = {
            "total_missing": int(total_missing),
            "total_cells": int(total_cells),
            "overall_missing_rate": float(overall_missing_rate),
            "variants_excluded": int(high_missing_variants),
            "values_imputed": int(imputed_count),
            "imputation_method": missing_imputation,
        }

        if method == "vanraden":
            # VanRaden method: K = (M * M.T) / (2 * sum(p * (1-p)))
            # M is centered genotypes (mean per variant = 0)

            # Compute allele frequencies per variant
            p = np.mean(X, axis=0) / 2.0  # Average genotype / 2 = allele frequency

            # Center genotypes (subtract 2*p per variant)
            M = X - 2.0 * p[np.newaxis, :]

            # Compute denominator: 2 * sum(p * (1-p))
            denominator = 2.0 * np.sum(p * (1.0 - p))
            if denominator == 0:
                denominator = 1.0  # Avoid division by zero

            # Compute kinship: K = (M * M.T) / denominator
            kinship = np.dot(M, M.T) / denominator

        elif method == "astle" or method == "astle-balding":
            # Astle-Balding method (similar to VanRaden but different normalization)
            p = np.mean(X, axis=0) / 2.0
            M = X - 2.0 * p[np.newaxis, :]
            denominator = np.sum(2.0 * p * (1.0 - p))
            if denominator == 0:
                denominator = 1.0
            kinship = np.dot(M, M.T) / denominator

        elif method == "yang":
            # Yang et al. method
            # Normalize by variant
            X_normalized = (X - np.mean(X, axis=0, keepdims=True)) / (
                np.std(X, axis=0, keepdims=True) + 1e-10
            )
            kinship = np.dot(X_normalized, X_normalized.T) / num_variants
        else:
            return {
                "status": "failed",
                "error": f"Unsupported kinship method: {method}",
            }

        logger.info(f"compute_kinship_matrix: Computed {num_samples}x{num_samples} kinship matrix")

        return {
            "status": "success",
            "kinship_matrix": kinship.tolist(),
            "method": method,
            "num_samples": num_samples,
            "missing_data_stats": missing_stats,
        }

    except Exception as exc:
        logger.error(f"compute_kinship_matrix: Error: {exc}")
        return {
            "status": "failed",
            "error": str(exc),
        }


def estimate_population_structure(
    vcf_path: str | Path,
    config: dict[str, Any],
    output_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Estimate population structure from VCF file.

    Args:
        vcf_path: Path to VCF file
        config: Structure configuration with compute_pca, n_components, etc.
        output_dir: Optional directory to write PCA and kinship results

    Returns:
        Dictionary with PCA and kinship results
    """
    logger.info(f"estimate_population_structure: Analyzing structure for {vcf_path}")

    compute_pca_flag = config.get("compute_pca", True)
    n_components = config.get("n_components", 10)
    compute_relatedness = config.get("compute_relatedness", True)
    kinship_method = config.get("kinship_method", "vanraden")

    # Parse VCF
    vcf_data = parse_vcf_full(vcf_path)
    samples = vcf_data["samples"]
    genotypes = vcf_data["genotypes"]

    if not genotypes:
        return {
            "status": "failed",
            "error": "No genotype data found",
        }

    results: dict[str, Any] = {
        "vcf_path": str(Path(vcf_path)),
        "num_samples": len(samples),
        "samples": samples,
    }

    # Compute PCA
    if compute_pca_flag:
        pca_result = compute_pca(genotypes, n_components=n_components)
        if pca_result.get("status") == "success":
            results["pca"] = pca_result
            logger.info("estimate_population_structure: PCA computation successful")
        else:
            logger.warning(f"estimate_population_structure: PCA failed: {pca_result.get('error')}")

    # Compute kinship
    if compute_relatedness:
        kinship_result = compute_kinship_matrix(genotypes, method=kinship_method)
        if kinship_result.get("status") == "success":
            results["kinship"] = kinship_result
            logger.info("estimate_population_structure: Kinship computation successful")
        else:
            logger.warning(
                f"estimate_population_structure: Kinship failed: {kinship_result.get('error')}"
            )

    # Write results if output directory provided
    if output_dir:
        out_dir = ensure_directory(output_dir)

        # Write PCA results as TSV
        if "pca" in results:
            pca_data = results["pca"]
            pcs = np.array(pca_data["pcs"])
            pca_table = [["sample_id"] + [f"PC{i+1}" for i in range(pcs.shape[1])]]
            for sample_idx, sample_id in enumerate(samples):
                pca_table.append([sample_id] + pcs[sample_idx, :].tolist())

            pca_path = out_dir / "pca_components.tsv"
            from ..core.io import write_tsv
            write_tsv(pca_table, pca_path)
            logger.info(f"estimate_population_structure: Wrote PCA to {pca_path}")

        # Write kinship matrix as TSV
        if "kinship" in results:
            kinship_path = out_dir / "kinship_matrix.tsv"
            kinship_matrix = np.array(results["kinship"]["kinship_matrix"])

            # Write header with sample IDs
            kinship_table = [[""] + samples]
            for sample_idx, sample_id in enumerate(samples):
                kinship_table.append([sample_id] + kinship_matrix[sample_idx, :].tolist())

            from ..core.io import write_tsv
            write_tsv(kinship_table, kinship_path)
            logger.info(f"estimate_population_structure: Wrote kinship matrix to {kinship_path}")

        # Write summary JSON
        summary_path = out_dir / "structure_summary.json"
        dump_json(results, summary_path, indent=2)
        logger.info(f"estimate_population_structure: Wrote summary to {summary_path}")

    results["status"] = "success"
    return results

