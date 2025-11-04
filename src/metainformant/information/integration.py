"""Integration functions for information theory with other modules.

This module provides wrapper functions and integration patterns for
applying information theory to data from DNA, RNA, single-cell, multi-omics,
and machine learning modules.
"""

from __future__ import annotations

from collections import Counter
from typing import Any

import numpy as np

from metainformant.information.analysis import information_profile
from metainformant.information.syntactic import (
    mutual_information,
    shannon_entropy_from_counts,
)


def _normalize_feature_indices(
    feature_indices: dict[str, int | list[int]] | None,
    platform: str
) -> list[int]:
    """Normalize feature indices to a list.
    
    Args:
        feature_indices: Dictionary mapping platform names to feature indices
        platform: Platform name to extract indices for
        
    Returns:
        List of feature indices (defaults to [0] if not specified)
    """
    if feature_indices is None:
        return [0]
    
    indices = feature_indices.get(platform, [0])
    if isinstance(indices, int):
        return [indices]
    return indices


def _calculate_cross_platform_mi(
    data1: np.ndarray,
    data2: np.ndarray,
    indices1: list[int],
    indices2: list[int],
    platform1: str,
    platform2: str
) -> dict[str, Any]:
    """Calculate mutual information between two platforms.
    
    Args:
        data1: First platform data matrix (samples x features)
        data2: Second platform data matrix (samples x features)
        indices1: Feature indices for first platform
        indices2: Feature indices for second platform
        platform1: Name of first platform
        platform2: Name of second platform
        
    Returns:
        Dictionary with MI matrix and statistics
        
    Raises:
        ValueError: If data shapes don't match or indices are out of range
    """
    if data1.shape[0] != data2.shape[0]:
        raise ValueError(
            f"Data matrices must have same number of samples: "
            f"{data1.shape[0]} vs {data2.shape[0]}"
        )
    
    # Validate feature indices
    if any(idx >= data1.shape[1] for idx in indices1):
        raise ValueError(
            f"{platform1} feature indices {indices1} out of range "
            f"(max: {data1.shape[1] - 1})"
        )
    if any(idx >= data2.shape[1] for idx in indices2):
        raise ValueError(
            f"{platform2} feature indices {indices2} out of range "
            f"(max: {data2.shape[1] - 1})"
        )
    
    # Calculate MI for all feature pairs
    n_samples = data1.shape[0]
    mi_matrix = []
    for idx1 in indices1:
        mi_row = []
        for idx2 in indices2:
            feature1 = [int(x) for x in data1[:n_samples, idx1]]
            feature2 = [int(y) for y in data2[:n_samples, idx2]]
            mi = mutual_information(feature1, feature2)
            mi_row.append(mi)
        mi_matrix.append(mi_row)
    
    # Always return matrix format with consistent structure
    result_key = f"{platform1}_{platform2}"
    return {
        f"{result_key}_mi_matrix": mi_matrix,
        f"{platform1}_feature_indices": indices1,
        f"{platform2}_feature_indices": indices2,
        f"{result_key}_mean_mi": float(np.mean(mi_matrix)),
        f"{result_key}_max_mi": float(np.max(mi_matrix)),
        f"{result_key}_min_mi": float(np.min(mi_matrix)),
    }


def dna_integration(
    sequences: dict[str, str] | list[str],
    k: int = 1,
    analysis_type: str = "entropy"
) -> dict[str, Any]:
    """Information-theoretic analysis of DNA sequences.
    
    Wrapper for analyzing DNA sequences from the DNA module.
    
    Args:
        sequences: Dictionary of sequence ID -> sequence, or list of sequences
        k: K-mer size for analysis
        analysis_type: Type of analysis ("entropy", "profile", "comparison")
        
    Returns:
        Analysis results dictionary
        
    Examples:
        >>> from metainformant.dna import sequences
        >>> dna_seqs = sequences.read_fasta("data/sequences.fasta")
        >>> results = dna_integration(dna_seqs, k=2)
    """
    # Convert to list if dictionary
    if isinstance(sequences, dict):
        seq_list = list(sequences.values())
    else:
        seq_list = sequences
    
    if not seq_list:
        raise ValueError("Sequences cannot be empty")
    
    if k < 1:
        raise ValueError(f"K-mer size k must be >= 1, got {k}")
    
    if analysis_type == "entropy":
        # Calculate entropy for each sequence
        results = []
        for i, seq in enumerate(seq_list):
            kmer_counts = Counter()
            for j in range(len(seq) - k + 1):
                kmer = seq[j : j + k]
                kmer_counts[kmer] += 1
            
            entropy = shannon_entropy_from_counts(kmer_counts)
            results.append({"sequence_index": i, "entropy": entropy})
        
        return {"entropy_analysis": results, "k": k}
    
    elif analysis_type == "profile":
        # Information profile
        profile = information_profile(seq_list, k=k)
        return {"profile": profile, "k": k}
    
    else:
        raise ValueError(f"Unknown analysis type: {analysis_type}")


def rna_integration(
    expression_data: np.ndarray,
    gene_names: list[str] | None = None,
    method: str = "entropy"
) -> dict[str, Any]:
    """Information-theoretic analysis of RNA expression data.
    
    Args:
        expression_data: Expression matrix (samples x genes)
        gene_names: Optional list of gene names
        method: Analysis method ("entropy", "mutual_information")
        
    Returns:
        Analysis results dictionary
        
    Examples:
        >>> import numpy as np
        >>> expression = np.random.randn(100, 50)  # 100 samples, 50 genes
        >>> results = rna_integration(expression, method="entropy")
    """
    if expression_data.size == 0:
        raise ValueError("Expression data cannot be empty")
    
    if len(expression_data.shape) != 2:
        raise ValueError(
            f"Expression data must be 2D array (samples x genes), got shape {expression_data.shape}"
        )
    
    if method == "entropy":
        # Calculate entropy for each gene
        gene_entropies = []
        for i in range(expression_data.shape[1]):
            gene_expr = expression_data[:, i]
            # Discretize for entropy calculation
            counts = Counter(gene_expr)
            entropy = shannon_entropy_from_counts(counts)
            gene_entropies.append({
                "gene_index": i,
                "gene_name": gene_names[i] if gene_names and i < len(gene_names) else f"Gene_{i}",
                "entropy": entropy,
            })
        
        return {
            "method": method,
            "gene_entropies": gene_entropies,
            "mean_entropy": float(np.mean([g["entropy"] for g in gene_entropies])),
        }
    
    elif method == "mutual_information":
        # Calculate pairwise MI between genes
        n_genes = expression_data.shape[1]
        mi_matrix = np.zeros((n_genes, n_genes))
        
        for i in range(n_genes):
            for j in range(i + 1, n_genes):
                # Discretize
                gene_i = [int(x) for x in expression_data[:, i]]
                gene_j = [int(y) for y in expression_data[:, j]]
                mi = mutual_information(gene_i, gene_j)
                mi_matrix[i, j] = mi
                mi_matrix[j, i] = mi
        
        return {
            "method": method,
            "mi_matrix": mi_matrix.tolist(),
            "mean_mi": float(np.mean(mi_matrix[np.triu_indices_from(mi_matrix, k=1)])),
        }
    
    else:
        raise ValueError(f"Unknown method: {method}")


def singlecell_integration(
    count_matrix: np.ndarray,
    cell_types: list[str] | None = None,
    method: str = "cell_type_entropy"
) -> dict[str, Any]:
    """Information-theoretic analysis of single-cell data.
    
    Args:
        count_matrix: Count matrix (cells x genes)
        cell_types: Optional list of cell type labels
        method: Analysis method ("cell_type_entropy", "gene_entropy")
        
    Returns:
        Analysis results dictionary
        
    Raises:
        ValueError: If count_matrix is empty or has invalid shape
    """
    if count_matrix.size == 0:
        raise ValueError("Count matrix cannot be empty")
    
    if len(count_matrix.shape) != 2:
        raise ValueError(
            f"Count matrix must be 2D array (cells x genes), got shape {count_matrix.shape}"
        )
    
    if method == "cell_type_entropy":
        # Calculate entropy of cell type distribution
        if cell_types:
            cell_type_counts = Counter(cell_types)
            entropy = shannon_entropy_from_counts(cell_type_counts)
            return {
                "method": method,
                "cell_type_entropy": entropy,
                "num_cell_types": len(set(cell_types)),
            }
        else:
            return {"method": method, "error": "Cell types not provided"}
    
    elif method == "gene_entropy":
        # Calculate entropy for each gene across cells
        gene_entropies = []
        for i in range(count_matrix.shape[1]):
            gene_counts = count_matrix[:, i]
            counts = Counter(gene_counts)
            entropy = shannon_entropy_from_counts(counts)
            gene_entropies.append({"gene_index": i, "entropy": entropy})
        
        return {
            "method": method,
            "gene_entropies": gene_entropies,
            "mean_entropy": float(np.mean([g["entropy"] for g in gene_entropies])),
        }
    
    else:
        raise ValueError(f"Unknown method: {method}")


def multiomics_integration(
    genomics_data: np.ndarray | None = None,
    transcriptomics_data: np.ndarray | None = None,
    proteomics_data: np.ndarray | None = None,
    method: str = "cross_platform_mi",
    feature_indices: dict[str, int | list[int]] | None = None
) -> dict[str, Any]:
    """Information-theoretic analysis across multiple omics platforms.
    
    Args:
        genomics_data: Genomic data matrix (samples x features)
        transcriptomics_data: Transcriptomic data matrix (samples x features)
        proteomics_data: Proteomic data matrix (samples x features)
        method: Analysis method ("cross_platform_mi", "platform_entropy")
        feature_indices: Optional dictionary mapping platform names to feature indices.
            For example: {"genomics": [0, 1, 2], "transcriptomics": [0, 3, 5]}.
            If None, uses first feature [0]. Single integers are converted to lists.
        
    Returns:
        Analysis results dictionary with consistent matrix format:
        - {platform1}_{platform2}_mi_matrix: MI matrix
        - {platform1}_feature_indices: List of indices used
        - {platform2}_feature_indices: List of indices used
        - {platform1}_{platform2}_mean_mi: Mean MI
        - {platform1}_{platform2}_max_mi: Max MI
        - {platform1}_{platform2}_min_mi: Min MI
        
    Examples:
        >>> import numpy as np
        >>> genomics = np.random.randn(100, 50)
        >>> transcriptomics = np.random.randn(100, 50)
        >>> # Use first feature (default)
        >>> results = multiomics_integration(genomics, transcriptomics, method="cross_platform_mi")
        >>> # Use specific features
        >>> results = multiomics_integration(
        ...     genomics, transcriptomics,
        ...     method="cross_platform_mi",
        ...     feature_indices={"genomics": [0, 1], "transcriptomics": [0, 2]}
        ... )
    """
    results: dict[str, Any] = {"method": method}
    
    if method == "platform_entropy":
        # Calculate entropy for each platform
        if genomics_data is not None:
            flat_genomics = genomics_data.flatten()
            counts = Counter(flat_genomics)
            results["genomics_entropy"] = shannon_entropy_from_counts(counts)
        
        if transcriptomics_data is not None:
            flat_transcriptomics = transcriptomics_data.flatten()
            counts = Counter(flat_transcriptomics)
            results["transcriptomics_entropy"] = shannon_entropy_from_counts(counts)
        
        if proteomics_data is not None:
            flat_proteomics = proteomics_data.flatten()
            counts = Counter(flat_proteomics)
            results["proteomics_entropy"] = shannon_entropy_from_counts(counts)
    
    elif method == "cross_platform_mi":
        # Calculate MI between genomics and transcriptomics
        if (
            genomics_data is not None
            and transcriptomics_data is not None
            and genomics_data.shape[0] == transcriptomics_data.shape[0]
        ):
            gen_indices = _normalize_feature_indices(feature_indices, "genomics")
            trans_indices = _normalize_feature_indices(feature_indices, "transcriptomics")
            mi_results = _calculate_cross_platform_mi(
                genomics_data, transcriptomics_data,
                gen_indices, trans_indices,
                "genomics", "transcriptomics"
            )
            results.update(mi_results)
        
        # Calculate MI between genomics and proteomics
        if (
            genomics_data is not None
            and proteomics_data is not None
            and genomics_data.shape[0] == proteomics_data.shape[0]
        ):
            gen_indices = _normalize_feature_indices(feature_indices, "genomics")
            prot_indices = _normalize_feature_indices(feature_indices, "proteomics")
            mi_results = _calculate_cross_platform_mi(
                genomics_data, proteomics_data,
                gen_indices, prot_indices,
                "genomics", "proteomics"
            )
            results.update(mi_results)
        
        # Calculate MI between transcriptomics and proteomics
        if (
            transcriptomics_data is not None
            and proteomics_data is not None
            and transcriptomics_data.shape[0] == proteomics_data.shape[0]
        ):
            trans_indices = _normalize_feature_indices(feature_indices, "transcriptomics")
            prot_indices = _normalize_feature_indices(feature_indices, "proteomics")
            mi_results = _calculate_cross_platform_mi(
                transcriptomics_data, proteomics_data,
                trans_indices, prot_indices,
                "transcriptomics", "proteomics"
            )
            results.update(mi_results)
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return results


def ml_integration(
    X: np.ndarray,
    y: np.ndarray,
    method: str = "feature_mi"
) -> dict[str, Any]:
    """Information-theoretic feature selection for ML.
    
    Args:
        X: Feature matrix (samples x features)
        y: Target vector
        method: Selection method ("feature_mi", "feature_entropy")
        
    Returns:
        Feature analysis results
        
    Examples:
        >>> import numpy as np
        >>> X = np.random.randn(100, 50)
        >>> y = np.random.randint(0, 2, 100)
        >>> results = ml_integration(X, y, method="feature_mi")
        
    Raises:
        ValueError: If X or y are empty, shapes don't match, or method is invalid
    """
    if X.size == 0:
        raise ValueError("Feature matrix X cannot be empty")
    
    if len(X.shape) != 2:
        raise ValueError(f"Feature matrix X must be 2D array (samples x features), got shape {X.shape}")
    
    if len(y) == 0:
        raise ValueError("Target vector y cannot be empty")
    
    if len(y) != X.shape[0]:
        raise ValueError(
            f"Target vector y length ({len(y)}) must match number of samples in X ({X.shape[0]})"
        )
    
    valid_methods = ["feature_mi", "feature_entropy"]
    if method not in valid_methods:
        raise ValueError(f"Method must be one of {valid_methods}, got {method}")
    
    if method == "feature_mi":
        # Calculate MI between each feature and target
        feature_mis = []
        for i in range(X.shape[1]):
            feature = X[:, i]
            # Discretize
            feature_discrete = [int(x) for x in feature]
            y_discrete = [int(yl) for yl in y]
            mi = mutual_information(feature_discrete, y_discrete)
            feature_mis.append({"feature_index": i, "mutual_information": mi})
        
        # Sort by MI
        feature_mis.sort(key=lambda x: x["mutual_information"], reverse=True)
        
        return {
            "method": method,
            "feature_mis": feature_mis,
            "top_features": [f["feature_index"] for f in feature_mis[:10]],
        }
    
    elif method == "feature_entropy":
        # Calculate entropy for each feature
        feature_entropies = []
        for i in range(X.shape[1]):
            feature = X[:, i]
            counts = Counter(feature)
            entropy = shannon_entropy_from_counts(counts)
            feature_entropies.append({"feature_index": i, "entropy": entropy})
        
        return {
            "method": method,
            "feature_entropies": feature_entropies,
            "mean_entropy": float(np.mean([f["entropy"] for f in feature_entropies])),
        }
    
    else:
        raise ValueError(f"Unknown method: {method}")

