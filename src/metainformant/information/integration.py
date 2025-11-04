"""Integration functions for information theory with other modules.

This module provides wrapper functions and integration patterns for
applying information theory to data from DNA, RNA, single-cell, multi-omics,
and machine learning modules.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from metainformant.information.analysis import information_profile
from metainformant.information.syntactic import (
    mutual_information,
    shannon_entropy_from_counts,
)


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
    
    if analysis_type == "entropy":
        # Calculate entropy for each sequence
        results = []
        for i, seq in enumerate(seq_list):
            from collections import Counter
            
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
    if method == "entropy":
        # Calculate entropy for each gene
        gene_entropies = []
        for i in range(expression_data.shape[1]):
            gene_expr = expression_data[:, i]
            # Discretize for entropy calculation
            from collections import Counter
            
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
    """
    if method == "cell_type_entropy":
        # Calculate entropy of cell type distribution
        if cell_types:
            from collections import Counter
            
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
            from collections import Counter
            
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
    method: str = "cross_platform_mi"
) -> dict[str, Any]:
    """Information-theoretic analysis across multiple omics platforms.
    
    Args:
        genomics_data: Genomic data matrix
        transcriptomics_data: Transcriptomic data matrix
        proteomics_data: Proteomic data matrix
        method: Analysis method ("cross_platform_mi", "platform_entropy")
        
    Returns:
        Analysis results dictionary
    """
    results: dict[str, Any] = {"method": method}
    
    if method == "platform_entropy":
        # Calculate entropy for each platform
        if genomics_data is not None:
            from collections import Counter
            
            # Flatten and discretize
            flat_genomics = genomics_data.flatten()
            counts = Counter(flat_genomics)
            results["genomics_entropy"] = shannon_entropy_from_counts(counts)
        
        if transcriptomics_data is not None:
            from collections import Counter
            
            flat_transcriptomics = transcriptomics_data.flatten()
            counts = Counter(flat_transcriptomics)
            results["transcriptomics_entropy"] = shannon_entropy_from_counts(counts)
        
        if proteomics_data is not None:
            from collections import Counter
            
            flat_proteomics = proteomics_data.flatten()
            counts = Counter(flat_proteomics)
            results["proteomics_entropy"] = shannon_entropy_from_counts(counts)
    
    elif method == "cross_platform_mi":
        # Calculate MI between platforms (if aligned)
        if (
            genomics_data is not None
            and transcriptomics_data is not None
            and genomics_data.shape[0] == transcriptomics_data.shape[0]
        ):
            # Sample-level MI
            n_samples = min(genomics_data.shape[0], transcriptomics_data.shape[0])
            # Use first feature from each platform
            gen_feature = [int(x) for x in genomics_data[:n_samples, 0]]
            trans_feature = [int(y) for y in transcriptomics_data[:n_samples, 0]]
            mi = mutual_information(gen_feature, trans_feature)
            results["genomics_transcriptomics_mi"] = mi
    
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
    """
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
            from collections import Counter
            
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

