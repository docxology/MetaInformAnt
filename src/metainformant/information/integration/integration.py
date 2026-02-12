"""Information-theoretic integration across biological data types.

This module provides functions for integrating information-theoretic analysis
across different biological data modalities (DNA, RNA, single-cell, multi-omics,
and machine learning feature analysis).
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List, Optional, Union

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Import information theory functions
try:
    from metainformant.information.metrics.core.syntactic import (
        mutual_information,
        shannon_entropy,
    )

    HAS_SYNTACTIC = True
except ImportError:
    HAS_SYNTACTIC = False
    logger.warning("Syntactic information functions not available")

try:
    from metainformant.information.metrics.analysis.analysis import information_profile

    HAS_ANALYSIS = True
except ImportError:
    HAS_ANALYSIS = False
    logger.warning("Analysis functions not available")


def _discretize_and_entropy(values: np.ndarray, n_bins: int = 10) -> float:
    """Discretize continuous values and compute Shannon entropy."""
    if len(values) < 2:
        return 0.0
    val_range = float(np.max(values) - np.min(values))
    if val_range == 0:
        return 0.0
    bins = np.linspace(np.min(values), np.max(values), n_bins)
    digitized = np.digitize(values, bins) - 1
    _, counts = np.unique(digitized, return_counts=True)
    if len(counts) <= 1:
        return 0.0
    total = counts.sum()
    probs = counts / total
    if HAS_SYNTACTIC:
        return shannon_entropy(probs)
    # Fallback: compute manually
    probs = probs[probs > 0]
    return float(-np.sum(probs * np.log2(probs)))


def dna_integration(
    sequences: Union[List[str], Dict[str, str]],
    *,
    k: int = 1,
    k_values: Optional[List[int]] = None,
    analysis_type: str = "entropy",
) -> Dict[str, Any]:
    """Integrate information-theoretic analysis for DNA sequences.

    Args:
        sequences: List of DNA sequences or dict mapping names to sequences
        k: k-mer size for analysis (convenience for single k)
        k_values: k-mer sizes to analyze (overrides k if provided)
        analysis_type: Type of analysis ('entropy' or 'profile')

    Returns:
        Dictionary with integrated DNA information analysis

    Raises:
        ValueError: If no sequences provided
    """
    # Accept dict or list
    if isinstance(sequences, dict):
        seq_names = list(sequences.keys())
        seq_list = list(sequences.values())
    else:
        seq_list = list(sequences)
        seq_names = [f"sequence_{i}" for i in range(len(seq_list))]

    if not seq_list:
        raise ValueError("No sequences provided")

    if k_values is None:
        k_values = [k]

    results: Dict[str, Any] = {
        "n_sequences": len(seq_list),
        "k_values": k_values,
        "analysis_type": analysis_type,
    }

    if analysis_type == "entropy":
        entropy_analysis: Dict[str, Any] = {}
        sequence_entropies: List[float] = []

        for name, seq in zip(seq_names, seq_list):
            nucleotides = ["A", "T", "C", "G"]
            counts = [seq.upper().count(nuc) for nuc in nucleotides]
            total = sum(counts)
            if total > 0 and HAS_SYNTACTIC:
                probs = [c / total for c in counts if c > 0]
                prob_sum = sum(probs)
                probs = [p / prob_sum for p in probs]
                entropy = shannon_entropy(probs)
            else:
                entropy = 0.0
            entropy_analysis[name] = {"entropy": entropy, "length": len(seq)}
            sequence_entropies.append(entropy)

        results["entropy_analysis"] = entropy_analysis

        if sequence_entropies:
            results["integrated_metrics"] = {
                "mean_entropy": float(np.mean(sequence_entropies)),
                "std_entropy": float(np.std(sequence_entropies)),
                "min_entropy": float(np.min(sequence_entropies)),
                "max_entropy": float(np.max(sequence_entropies)),
            }
        else:
            results["integrated_metrics"] = {}

    elif analysis_type == "profile":
        if HAS_ANALYSIS:
            try:
                profile = information_profile(seq_list, k=max(k_values))
                results["profile"] = profile
            except Exception as e:
                logger.warning(f"Failed profile analysis: {e}")
                results["profile"] = {}
        else:
            results["profile"] = {}

    return results


def rna_integration(
    expression_matrix: Any,
    *,
    method: str = "entropy",
    normalize: bool = True,
) -> Dict[str, Any]:
    """Integrate information-theoretic analysis for RNA expression data.

    Args:
        expression_matrix: Expression matrix (numpy array or DataFrame).
            Rows = samples, columns = genes.
        method: Analysis method ('entropy' or 'mutual_information')
        normalize: Whether to normalize expression values

    Returns:
        Dictionary with RNA information analysis
    """
    # Convert to numpy array
    if hasattr(expression_matrix, "values"):
        data = np.array(expression_matrix.values, dtype=float)
    else:
        data = np.array(expression_matrix, dtype=float)

    n_samples, n_genes = data.shape
    results: Dict[str, Any] = {"n_samples": n_samples, "n_genes": n_genes, "method": method}

    if method == "entropy":
        gene_entropies: List[float] = []
        for gene_idx in range(n_genes):
            gene_expr = data[:, gene_idx].copy()
            if normalize and np.sum(np.abs(gene_expr)) > 0:
                gene_expr = gene_expr / np.sum(np.abs(gene_expr))
            entropy = _discretize_and_entropy(gene_expr)
            gene_entropies.append(entropy)

        results["gene_entropies"] = gene_entropies
        results["integrated_metrics"] = {
            "gene_entropy_mean": float(np.mean(gene_entropies)) if gene_entropies else 0.0,
            "gene_entropy_std": float(np.std(gene_entropies)) if gene_entropies else 0.0,
        }

    elif method == "mutual_information":
        # Compute pairwise MI between genes (discretized)
        mi_matrix = np.zeros((n_genes, n_genes))
        discretized = np.zeros_like(data, dtype=int)
        for g in range(n_genes):
            vals = data[:, g]
            val_range = np.max(vals) - np.min(vals)
            if val_range > 0:
                bins = np.linspace(np.min(vals), np.max(vals), 10)
                discretized[:, g] = np.digitize(vals, bins) - 1
            else:
                discretized[:, g] = 0

        for i in range(n_genes):
            for j in range(i, n_genes):
                if i == j:
                    mi_matrix[i][j] = 0.0
                else:
                    if HAS_SYNTACTIC:
                        mi_val = mutual_information(
                            list(discretized[:, i]),
                            list(discretized[:, j]),
                        )
                    else:
                        mi_val = 0.0
                    mi_matrix[i][j] = mi_val
                    mi_matrix[j][i] = mi_val

        results["mi_matrix"] = mi_matrix.tolist()

    return results


def singlecell_integration(
    data: Any,
    *,
    cell_types: Optional[List[str]] = None,
    method: str = "gene_entropy",
) -> Dict[str, Any]:
    """Integrate information-theoretic analysis for single-cell data.

    Args:
        data: Count matrix (numpy array, shape: cells x genes) or AnnData object
        cell_types: Optional list of cell type labels (one per cell)
        method: Analysis method ('gene_entropy' or 'cell_type_entropy')

    Returns:
        Dictionary with single-cell information analysis
    """
    # Extract expression matrix
    if hasattr(data, "X"):
        # AnnData object
        expr_matrix = data.X.toarray() if hasattr(data.X, "toarray") else np.array(data.X)
        n_cells = data.n_obs if hasattr(data, "n_obs") else expr_matrix.shape[0]
        n_genes = data.n_vars if hasattr(data, "n_vars") else expr_matrix.shape[1]
    else:
        expr_matrix = np.array(data, dtype=float)
        n_cells, n_genes = expr_matrix.shape

    results: Dict[str, Any] = {
        "n_cells": n_cells,
        "n_genes": n_genes,
        "method": method,
    }

    if method == "gene_entropy":
        gene_entropies: List[float] = []
        for gene_idx in range(n_genes):
            gene_expr = expr_matrix[:, gene_idx].astype(float)
            if gene_expr.sum() > 0:
                gene_expr_norm = gene_expr / gene_expr.sum()
                entropy = _discretize_and_entropy(gene_expr_norm)
            else:
                entropy = 0.0
            gene_entropies.append(entropy)

        results["gene_entropies"] = gene_entropies

    elif method == "cell_type_entropy":
        if cell_types is None:
            raise ValueError("cell_types required for cell_type_entropy method")

        unique_types = list(set(cell_types))
        results["num_cell_types"] = len(unique_types)

        # Compute entropy of cell type distribution
        type_counts = Counter(cell_types)
        total = sum(type_counts.values())
        probs = [count / total for count in type_counts.values()]
        if HAS_SYNTACTIC:
            ct_entropy = shannon_entropy(probs)
        else:
            probs_arr = np.array([p for p in probs if p > 0])
            ct_entropy = float(-np.sum(probs_arr * np.log2(probs_arr)))

        results["cell_type_entropy"] = ct_entropy

        # Per-type gene entropy
        type_entropies: Dict[str, float] = {}
        for ct in unique_types:
            mask = [i for i, t in enumerate(cell_types) if t == ct]
            subset = expr_matrix[mask, :]
            type_gene_entropies = []
            for gene_idx in range(n_genes):
                gene_expr = subset[:, gene_idx].astype(float)
                if gene_expr.sum() > 0:
                    gene_expr_norm = gene_expr / gene_expr.sum()
                    entropy = _discretize_and_entropy(gene_expr_norm)
                else:
                    entropy = 0.0
                type_gene_entropies.append(entropy)
            type_entropies[ct] = float(np.mean(type_gene_entropies))

        results["per_type_entropy"] = type_entropies

    return results


def multiomics_integration(
    *,
    omics_data: Optional[Dict[str, Any]] = None,
    method: str = "platform_entropy",
    feature_indices: Optional[Dict[str, Any]] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Integrate information across multiple omics platforms.

    Args:
        omics_data: Dictionary mapping omics type names to data arrays
        method: Integration method ('platform_entropy' or 'cross_platform_mi')
        feature_indices: Optional dict mapping platform names to feature indices
        **kwargs: Named platform data (e.g., genomics_data=..., transcriptomics_data=...)

    Returns:
        Dictionary with multi-omics information integration
    """
    # Build platform data from kwargs or omics_data
    # Separate sequence data (lists of strings) from numeric arrays
    platforms: Dict[str, np.ndarray] = {}
    sequence_platforms: Dict[str, List[str]] = {}

    if omics_data is not None:
        for name, data in omics_data.items():
            if isinstance(data, list) and data and isinstance(data[0], str):
                sequence_platforms[name] = data
            else:
                platforms[name] = np.array(data, dtype=float)

    # Extract named platform data from kwargs
    for key, value in kwargs.items():
        if key.endswith("_data"):
            platform_name = key[:-5]  # Remove "_data" suffix
            if isinstance(value, list) and value and isinstance(value[0], str):
                sequence_platforms[platform_name] = value
            else:
                platforms[platform_name] = np.array(value, dtype=float)

    all_names = list(sequence_platforms.keys()) + list(platforms.keys())

    results: Dict[str, Any] = {
        "omics_types": all_names,
        "method": method,
    }

    if method == "platform_entropy":
        # Handle sequence platforms (DNA) via dna_integration
        for platform_name, seqs in sequence_platforms.items():
            try:
                dna_result = dna_integration(seqs)
                metrics = dna_result.get("integrated_metrics", {})
                results[f"{platform_name}_entropy"] = metrics.get("mean_entropy", 0.0)
            except (ValueError, Exception) as e:
                logger.warning(f"Failed {platform_name} integration: {e}")
                results[f"{platform_name}_entropy"] = 0.0

        # Handle numeric platforms
        for platform_name, data in platforms.items():
            n_features = data.shape[1] if data.ndim > 1 else 1
            feature_entropies = []
            for f in range(n_features):
                col = data[:, f] if data.ndim > 1 else data
                entropy = _discretize_and_entropy(col)
                feature_entropies.append(entropy)
            results[f"{platform_name}_entropy"] = float(np.mean(feature_entropies))

    elif method == "cross_platform_mi":
        # Compute pairwise MI between platforms
        platform_names = list(platforms.keys())

        for i in range(len(platform_names)):
            for j in range(i + 1, len(platform_names)):
                name_i = platform_names[i]
                name_j = platform_names[j]
                data_i = platforms[name_i]
                data_j = platforms[name_j]

                # Get feature indices
                idx_i = feature_indices.get(name_i, None) if feature_indices else None
                idx_j = feature_indices.get(name_j, None) if feature_indices else None

                # Normalize indices to lists
                if idx_i is not None:
                    if isinstance(idx_i, (int, np.integer)):
                        idx_i = [int(idx_i)]
                    idx_i = list(idx_i)
                else:
                    idx_i = [0]  # Default to first feature

                if idx_j is not None:
                    if isinstance(idx_j, (int, np.integer)):
                        idx_j = [int(idx_j)]
                    idx_j = list(idx_j)
                else:
                    idx_j = [0]

                # Compute MI matrix between selected features
                mi_matrix = []
                for fi in idx_i:
                    row = []
                    for fj in idx_j:
                        col_i = data_i[:, fi] if data_i.ndim > 1 else data_i
                        col_j = data_j[:, fj] if data_j.ndim > 1 else data_j
                        # Discretize both
                        bins_i = np.linspace(np.min(col_i), np.max(col_i), 10)
                        bins_j = np.linspace(np.min(col_j), np.max(col_j), 10)
                        disc_i = list(np.digitize(col_i, bins_i) - 1)
                        disc_j = list(np.digitize(col_j, bins_j) - 1)
                        if HAS_SYNTACTIC:
                            mi_val = mutual_information(disc_i, disc_j)
                        else:
                            mi_val = 0.0
                        row.append(float(mi_val))
                    mi_matrix.append(row)

                key_prefix = f"{name_i}_{name_j}"
                results[f"{key_prefix}_mi_matrix"] = mi_matrix
                flat = [v for row in mi_matrix for v in row]
                results[f"{key_prefix}_mean_mi"] = float(np.mean(flat))
                results[f"{key_prefix}_max_mi"] = float(np.max(flat))
                results[f"{key_prefix}_min_mi"] = float(np.min(flat))

    return results


def ml_integration(
    features: Any,
    labels: Any,
    method: str = "feature_mi",
) -> Dict[str, Any]:
    """Integrate information theory with machine learning feature analysis.

    Args:
        features: Feature matrix (numpy array or DataFrame)
        labels: Target labels (classification or regression)
        method: Integration method ('feature_mi', 'mutual_info', 'feature_entropy',
                'entropy', 'correlation')

    Returns:
        Dictionary with ML-information integration results
    """
    try:
        from sklearn.feature_selection import mutual_info_classif, mutual_info_regression
    except ImportError:
        logger.warning("scikit-learn not available for ML integration")
        return {"error": "scikit-learn required"}

    # Convert to numpy arrays
    if hasattr(features, "values"):
        X = features.values
    else:
        X = np.array(features, dtype=float)

    if hasattr(labels, "values"):
        y = labels.values
    else:
        y = np.array(labels)

    n_samples, n_features = X.shape
    results: Dict[str, Any] = {"n_features": n_features, "n_samples": n_samples, "method": method}

    try:
        if method in ("feature_mi", "mutual_info"):
            # Mutual information between features and target
            if len(np.unique(y)) <= 10:
                scores = mutual_info_classif(X, y)
            else:
                scores = mutual_info_regression(X, y)

            results["feature_mis"] = scores.tolist()

            # Top features (sorted by MI score)
            top_indices = np.argsort(scores)[::-1][: min(10, n_features)]
            results["top_features"] = [{"index": int(idx), "mi": float(scores[idx])} for idx in top_indices]
            results["feature_scores"] = {f"feature_{i}": float(scores[i]) for i in top_indices}

        elif method in ("feature_entropy", "entropy"):
            feature_entropies = []
            for col in range(n_features):
                feature_vals = X[:, col]
                entropy = _discretize_and_entropy(feature_vals)
                feature_entropies.append(entropy)

            results["feature_entropies"] = feature_entropies
            results["feature_scores"] = {
                "mean_feature_entropy": float(np.mean(feature_entropies)),
                "std_feature_entropy": float(np.std(feature_entropies)),
                "max_feature_entropy": float(np.max(feature_entropies)),
                "min_feature_entropy": float(np.min(feature_entropies)),
            }

        elif method == "correlation":
            correlations = []
            for col in range(n_features):
                feature_vals = X[:, col]
                try:
                    corr = np.corrcoef(feature_vals, y.astype(float))[0, 1]
                    if not np.isnan(corr):
                        correlations.append(abs(float(corr)))
                except (ValueError, TypeError, FloatingPointError):
                    continue

            if correlations:
                results["feature_scores"] = {
                    "mean_correlation": float(np.mean(correlations)),
                    "max_correlation": float(np.max(correlations)),
                    "n_significant_features": sum(1 for c in correlations if c > 0.1),
                }
            else:
                results["feature_scores"] = {}

    except Exception as e:
        logger.warning(f"Failed ML integration: {e}")
        results["error"] = str(e)
        results["feature_scores"] = {}

    return results
