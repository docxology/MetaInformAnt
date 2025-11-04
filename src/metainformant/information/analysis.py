"""High-level information-theoretic analysis functions.

This module provides analysis workflows and convenience functions for
applying information theory to biological data.
"""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from metainformant.information.syntactic import (
    conditional_entropy,
    joint_entropy,
    kl_divergence,
    mutual_information,
    shannon_entropy,
    shannon_entropy_from_counts,
    transfer_entropy,
)


def information_profile(
    sequences: list[str],
    k: int = 1
) -> dict[str, Any]:
    """Calculate comprehensive information profile for a set of sequences.
    
    Computes multiple information-theoretic measures including entropy,
    k-mer frequencies, and sequence complexity.
    
    Args:
        sequences: List of sequences (DNA, RNA, or protein)
        k: K-mer size for analysis
        
    Returns:
        Dictionary containing:
        - entropy: Shannon entropy of k-mer frequencies
        - kmer_frequencies: Dictionary of k-mer counts
        - unique_kmers: Number of unique k-mers
        - sequence_complexity: Normalized entropy (entropy / max_entropy)
        
    Examples:
        >>> seqs = ["ATCG", "ATCG", "AAAA"]
        >>> profile = information_profile(seqs, k=1)
        >>> profile["entropy"] > 0
        True
    """
    if not sequences:
        return {
            "entropy": 0.0,
            "kmer_frequencies": {},
            "unique_kmers": 0,
            "sequence_complexity": 0.0,
        }
    
    # Count k-mers across all sequences
    kmer_counts: dict[str, int] = Counter()
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            kmer_counts[kmer] += 1
    
    # Calculate entropy
    entropy = shannon_entropy_from_counts(kmer_counts)
    
    # Calculate maximum possible entropy (uniform distribution)
    unique_kmers = len(kmer_counts)
    max_entropy = math.log2(unique_kmers) if unique_kmers > 0 else 0.0
    
    # Sequence complexity (normalized entropy)
    complexity = entropy / max_entropy if max_entropy > 0 else 0.0
    
    return {
        "entropy": entropy,
        "kmer_frequencies": dict(kmer_counts),
        "unique_kmers": unique_kmers,
        "sequence_complexity": complexity,
    }


def information_signature(
    data: np.ndarray | list[list[float]],
    method: str = "entropy"
) -> dict[str, Any]:
    """Calculate information signature for multivariate data.
    
    Computes information-theoretic measures for feature vectors or
    time series data.
    
    Args:
        data: 2D array where rows are samples and columns are features
        method: Analysis method ("entropy", "mutual_information", "total_correlation")
        
    Returns:
        Dictionary containing information signature metrics
    """
    if isinstance(data, list):
        data = np.array(data)
    
    if data.size == 0:
        return {"signature": {}, "method": method}
    
    n_samples, n_features = data.shape
    
    if method == "entropy":
        # Calculate entropy for each feature
        entropies = []
        for i in range(n_features):
            feature_data = data[:, i]
            # Discretize for entropy calculation
            counts = Counter(feature_data)
            entropy = shannon_entropy_from_counts(counts)
            entropies.append(entropy)
        
        return {
            "signature": {
                "feature_entropies": entropies,
                "mean_entropy": np.mean(entropies),
                "std_entropy": np.std(entropies),
            },
            "method": method,
        }
    
    elif method == "mutual_information":
        # Calculate pairwise mutual information
        mi_matrix = np.zeros((n_features, n_features))
        for i in range(n_features):
            for j in range(i + 1, n_features):
                # Discretize for MI calculation
                x_discrete = [int(x) for x in data[:, i]]
                y_discrete = [int(y) for y in data[:, j]]
                mi = mutual_information(x_discrete, y_discrete)
                mi_matrix[i, j] = mi
                mi_matrix[j, i] = mi
        
        return {
            "signature": {
                "mi_matrix": mi_matrix.tolist(),
                "mean_mi": np.mean(mi_matrix[np.triu_indices_from(mi_matrix, k=1)]),
            },
            "method": method,
        }
    
    else:
        raise ValueError(f"Unknown method: {method}")


def analyze_sequence_information(
    sequence: str,
    k_values: list[int] | None = None
) -> dict[str, Any]:
    """Analyze information content of a single sequence.
    
    Args:
        sequence: Sequence to analyze
        k_values: List of k-mer sizes to analyze (default: [1, 2, 3])
        
    Returns:
        Dictionary with information analysis results
    """
    if k_values is None:
        k_values = [1, 2, 3]
    
    results: dict[str, Any] = {
        "sequence_length": len(sequence),
        "kmer_analyses": {},
    }
    
    for k in k_values:
        if k > len(sequence):
            continue
        
        # Count k-mers
        kmer_counts: dict[str, int] = Counter()
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i : i + k]
            kmer_counts[kmer] += 1
        
        # Calculate entropy
        entropy = shannon_entropy_from_counts(kmer_counts)
        
        # Maximum entropy
        unique_kmers = len(kmer_counts)
        max_entropy = math.log2(unique_kmers) if unique_kmers > 0 else 0.0
        
        results["kmer_analyses"][k] = {
            "entropy": entropy,
            "max_entropy": max_entropy,
            "normalized_entropy": entropy / max_entropy if max_entropy > 0 else 0.0,
            "unique_kmers": unique_kmers,
            "total_kmers": sum(kmer_counts.values()),
        }
    
    return results


def compare_sequences_information(
    seq1: str,
    seq2: str,
    k: int = 1
) -> dict[str, Any]:
    """Compare two sequences using information-theoretic measures.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        k: K-mer size for comparison
        
    Returns:
        Dictionary containing:
        - entropy_1, entropy_2: Individual entropies
        - kl_divergence: KL divergence between k-mer distributions
        - joint_entropy: Joint entropy of both sequences
        - mutual_information: Mutual information between sequences
    """
    # Count k-mers for both sequences
    kmer_counts1: dict[str, int] = Counter()
    kmer_counts2: dict[str, int] = Counter()
    
    for i in range(len(seq1) - k + 1):
        kmer = seq1[i : i + k]
        kmer_counts1[kmer] += 1
    
    for i in range(len(seq2) - k + 1):
        kmer = seq2[i : i + k]
        kmer_counts2[kmer] += 1
    
    # Get all unique k-mers
    all_kmers = set(kmer_counts1.keys()) | set(kmer_counts2.keys())
    
    # Create probability distributions
    total1 = sum(kmer_counts1.values())
    total2 = sum(kmer_counts2.values())
    
    p1 = [kmer_counts1.get(kmer, 0) / total1 if total1 > 0 else 0 for kmer in sorted(all_kmers)]
    p2 = [kmer_counts2.get(kmer, 0) / total2 if total2 > 0 else 0 for kmer in sorted(all_kmers)]
    
    # Calculate measures
    entropy1 = shannon_entropy_from_counts(kmer_counts1)
    entropy2 = shannon_entropy_from_counts(kmer_counts2)
    kl = kl_divergence(p1, p2)
    
    # For joint entropy and MI, we need aligned sequences
    # For simplicity, use k-mer co-occurrence
    joint_counts: dict[tuple[str, str], int] = {}
    for i in range(min(len(seq1), len(seq2)) - k + 1):
        kmer1 = seq1[i : i + k]
        kmer2 = seq2[i : i + k]
        key = (kmer1, kmer2)
        joint_counts[key] = joint_counts.get(key, 0) + 1
    
    total_joint = sum(joint_counts.values())
    joint_probs = {k: v / total_joint for k, v in joint_counts.items()} if total_joint > 0 else {}
    h_joint = joint_entropy(joint_probs)
    
    mi = entropy1 + entropy2 - h_joint if h_joint > 0 else 0.0
    
    return {
        "entropy_1": entropy1,
        "entropy_2": entropy2,
        "kl_divergence": kl,
        "joint_entropy": h_joint,
        "mutual_information": max(0.0, mi),
    }

