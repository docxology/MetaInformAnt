"""Information-theoretic integration across biological data types.

This module provides functions for integrating information-theoretic analysis
across different biological data modalities (DNA, RNA, single-cell, multi-omics).
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Import information theory functions
try:
    from .syntactic import shannon_entropy, mutual_information

    HAS_SYNTACTIC = True
except ImportError:
    HAS_SYNTACTIC = False
    logger.warning("Syntactic information functions not available")

try:
    from .analysis import information_profile

    HAS_ANALYSIS = True
except ImportError:
    HAS_ANALYSIS = False
    logger.warning("Analysis functions not available")


def dna_integration(sequences: List[str], k_values: Optional[List[int]] = None) -> Dict[str, Any]:
    """Integrate information-theoretic analysis for DNA sequences.

    Args:
        sequences: List of DNA sequences
        k_values: k-mer sizes to analyze (default: [1, 2, 3])

    Returns:
        Dictionary with integrated DNA information analysis
    """
    if not sequences:
        raise ValueError("No sequences provided")

    if k_values is None:
        k_values = [1, 2, 3]

    results = {"n_sequences": len(sequences), "k_values": k_values, "sequence_info": {}, "integrated_metrics": {}}

    # Individual sequence analysis
    sequence_entropies = []
    for i, seq in enumerate(sequences):
        if HAS_ANALYSIS:
            try:
                seq_info = information_profile([seq], k=max(k_values))
                results["sequence_info"][f"sequence_{i}"] = seq_info
                if "entropy" in seq_info:
                    sequence_entropies.append(seq_info["entropy"])
            except Exception as e:
                logger.warning(f"Failed to analyze sequence {i}: {e}")
        else:
            # Basic entropy calculation
            if HAS_SYNTACTIC:
                try:
                    # Simple nucleotide frequency entropy
                    nucleotides = ["A", "T", "C", "G"]
                    counts = [seq.upper().count(nuc) for nuc in nucleotides]
                    total = sum(counts)
                    if total > 0:
                        probs = [c / total for c in counts]
                        entropy = shannon_entropy(probs)
                        sequence_entropies.append(entropy)
                        results["sequence_info"][f"sequence_{i}"] = {"entropy": entropy}
                except Exception as e:
                    logger.warning(f"Failed to analyze sequence {i}: {e}")

    # Integrated metrics
    if sequence_entropies:
        results["integrated_metrics"] = {
            "mean_entropy": float(np.mean(sequence_entropies)),
            "std_entropy": float(np.std(sequence_entropies)),
            "min_entropy": float(np.min(sequence_entropies)),
            "max_entropy": float(np.max(sequence_entropies)),
        }

    return results


def rna_integration(expression_matrix: Any, normalize: bool = True) -> Dict[str, Any]:
    """Integrate information-theoretic analysis for RNA expression data.

    Args:
        expression_matrix: Expression matrix (pandas DataFrame or numpy array)
        normalize: Whether to normalize expression values

    Returns:
        Dictionary with RNA information analysis
    """
    try:
        import pandas as pd
    except ImportError:
        logger.warning("pandas not available for RNA integration")
        return {"error": "pandas required"}

    # Convert to DataFrame if needed
    if isinstance(expression_matrix, np.ndarray):
        df = pd.DataFrame(expression_matrix)
    elif hasattr(expression_matrix, "values"):
        df = pd.DataFrame(expression_matrix.values)
    else:
        df = pd.DataFrame(expression_matrix)

    results = {"n_genes": df.shape[0], "n_samples": df.shape[1], "integrated_metrics": {}}

    # Expression entropy analysis
    gene_entropies = []
    sample_entropies = []

    for gene_idx in range(df.shape[0]):
        gene_expr = df.iloc[gene_idx, :].values
        if normalize and gene_expr.sum() > 0:
            gene_expr = gene_expr / gene_expr.sum()

        if HAS_SYNTACTIC and len(gene_expr) > 1:
            try:
                # Discretize expression for entropy calculation
                bins = np.linspace(np.min(gene_expr), np.max(gene_expr), 10)
                digitized = np.digitize(gene_expr, bins) - 1
                unique_vals, counts = np.unique(digitized, return_counts=True)
                probs = counts / len(counts)
                entropy = shannon_entropy(probs)
                gene_entropies.append(entropy)
            except Exception as e:
                logger.debug(f"Failed gene entropy calculation: {e}")

    for sample_idx in range(df.shape[1]):
        sample_expr = df.iloc[:, sample_idx].values
        if normalize and sample_expr.sum() > 0:
            sample_expr = sample_expr / sample_expr.sum()

        if HAS_SYNTACTIC and len(sample_expr) > 1:
            try:
                # Discretize expression for entropy calculation
                bins = np.linspace(np.min(sample_expr), np.max(sample_expr), 10)
                digitized = np.digitize(sample_expr, bins) - 1
                unique_vals, counts = np.unique(digitized, return_counts=True)
                probs = counts / len(counts)
                entropy = shannon_entropy(probs)
                sample_entropies.append(entropy)
            except Exception as e:
                logger.debug(f"Failed sample entropy calculation: {e}")

    # Integrated metrics
    results["integrated_metrics"] = {
        "gene_entropy_mean": float(np.mean(gene_entropies)) if gene_entropies else None,
        "gene_entropy_std": float(np.std(gene_entropies)) if gene_entropies else None,
        "sample_entropy_mean": float(np.mean(sample_entropies)) if sample_entropies else None,
        "sample_entropy_std": float(np.std(sample_entropies)) if sample_entropies else None,
    }

    return results


def singlecell_integration(adata: Any) -> Dict[str, Any]:
    """Integrate information-theoretic analysis for single-cell data.

    Args:
        adata: AnnData object with single-cell data

    Returns:
        Dictionary with single-cell information analysis
    """
    try:
        import scanpy as sc

        HAS_SCANPY = True
    except ImportError:
        HAS_SCANPY = False
        logger.warning("scanpy not available for single-cell integration")

    results = {"n_cells": getattr(adata, "n_obs", 0), "n_genes": getattr(adata, "n_vars", 0), "integrated_metrics": {}}

    # Basic expression entropy
    if hasattr(adata, "X") and adata.X is not None:
        try:
            # Convert to dense if sparse
            expr_matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X

            # Gene expression entropy
            gene_entropies = []
            for gene_idx in range(min(expr_matrix.shape[1], 1000)):  # Limit for performance
                gene_expr = expr_matrix[:, gene_idx]
                if gene_expr.sum() > 0:
                    # Normalize and discretize
                    gene_expr = gene_expr / gene_expr.sum()
                    bins = np.linspace(np.min(gene_expr), np.max(gene_expr), 10)
                    digitized = np.digitize(gene_expr, bins) - 1
                    unique_vals, counts = np.unique(digitized, return_counts=True)
                    if len(counts) > 1:
                        probs = counts / len(counts)
                        if HAS_SYNTACTIC:
                            entropy = shannon_entropy(probs)
                            gene_entropies.append(entropy)

            results["integrated_metrics"] = {
                "gene_expression_entropy_mean": float(np.mean(gene_entropies)) if gene_entropies else None,
                "gene_expression_entropy_std": float(np.std(gene_entropies)) if gene_entropies else None,
            }

        except Exception as e:
            logger.warning(f"Failed single-cell expression analysis: {e}")

    return results


def multiomics_integration(omics_data: Dict[str, Any]) -> Dict[str, Any]:
    """Integrate information across multiple omics platforms.

    Args:
        omics_data: Dictionary with different omics data types

    Returns:
        Dictionary with multi-omics information integration
    """
    results = {"omics_types": list(omics_data.keys()), "integrated_metrics": {}, "cross_platform_analysis": {}}

    # Analyze each omics type
    omics_results = {}
    for omics_type, data in omics_data.items():
        try:
            if omics_type.lower() in ["dna", "sequences"]:
                omics_results[omics_type] = dna_integration(data)
            elif omics_type.lower() in ["rna", "expression"]:
                omics_results[omics_type] = rna_integration(data)
            elif omics_type.lower() in ["singlecell", "sc"]:
                omics_results[omics_type] = singlecell_integration(data)
            else:
                logger.warning(f"Unknown omics type: {omics_type}")
                continue
        except Exception as e:
            logger.warning(f"Failed to analyze {omics_type}: {e}")

    results["omics_results"] = omics_results

    # Cross-platform integration
    if len(omics_results) > 1:
        # Simple integration: average entropies across platforms
        entropy_measures = []
        for omics_type, result in omics_results.items():
            metrics = result.get("integrated_metrics", {})
            for key, value in metrics.items():
                if "entropy" in key and value is not None:
                    entropy_measures.append(value)

        if entropy_measures:
            results["cross_platform_analysis"] = {
                "mean_entropy_across_platforms": float(np.mean(entropy_measures)),
                "std_entropy_across_platforms": float(np.std(entropy_measures)),
                "n_entropy_measures": len(entropy_measures),
            }

    return results


def ml_integration(features: Any, labels: Any, method: str = "mutual_info") -> Dict[str, Any]:
    """Integrate information theory with machine learning feature analysis.

    Args:
        features: Feature matrix
        labels: Target labels
        method: Integration method ('mutual_info', 'entropy', 'correlation')

    Returns:
        Dictionary with ML-information integration results
    """
    try:
        import pandas as pd
        from sklearn.feature_selection import mutual_info_classif, mutual_info_regression
    except ImportError:
        logger.warning("scikit-learn not available for ML integration")
        return {"error": "scikit-learn required"}

    # Convert to numpy arrays
    if hasattr(features, "values"):
        X = features.values
    else:
        X = np.array(features)

    if hasattr(labels, "values"):
        y = labels.values
    else:
        y = np.array(labels)

    results = {"n_features": X.shape[1], "n_samples": X.shape[0], "method": method, "feature_scores": {}}

    try:
        if method == "mutual_info":
            # Mutual information between features and target
            if len(np.unique(y)) <= 10:  # Classification
                scores = mutual_info_classif(X, y)
            else:  # Regression
                scores = mutual_info_regression(X, y)

            # Store top features
            feature_indices = np.argsort(scores)[::-1]
            results["feature_scores"] = {f"feature_{i}": float(scores[i]) for i in feature_indices[:10]}  # Top 10

        elif method == "entropy":
            # Feature entropy analysis
            feature_entropies = []
            for col in range(X.shape[1]):
                feature_vals = X[:, col]
                # Discretize for entropy calculation
                bins = np.linspace(np.min(feature_vals), np.max(feature_vals), 10)
                digitized = np.digitize(feature_vals, bins) - 1
                unique_vals, counts = np.unique(digitized, return_counts=True)
                if len(counts) > 1 and HAS_SYNTACTIC:
                    probs = counts / len(counts)
                    entropy = shannon_entropy(probs)
                    feature_entropies.append(entropy)

            if feature_entropies:
                results["feature_scores"] = {
                    "mean_feature_entropy": float(np.mean(feature_entropies)),
                    "std_feature_entropy": float(np.std(feature_entropies)),
                    "max_feature_entropy": float(np.max(feature_entropies)),
                    "min_feature_entropy": float(np.min(feature_entropies)),
                }

        elif method == "correlation":
            # Correlation with target
            correlations = []
            for col in range(X.shape[1]):
                feature_vals = X[:, col]
                try:
                    corr = np.corrcoef(feature_vals, y)[0, 1]
                    if not np.isnan(corr):
                        correlations.append(abs(corr))
                except (ValueError, TypeError, FloatingPointError):
                    continue

            if correlations:
                results["feature_scores"] = {
                    "mean_correlation": float(np.mean(correlations)),
                    "max_correlation": float(np.max(correlations)),
                    "n_significant_features": sum(1 for c in correlations if c > 0.1),
                }

    except Exception as e:
        logger.warning(f"Failed ML integration: {e}")
        results["error"] = str(e)

    return results
