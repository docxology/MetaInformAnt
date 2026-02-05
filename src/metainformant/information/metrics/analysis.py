"""High-level information analysis functions for biological sequences.

This module provides comprehensive analysis functions that combine multiple
information-theoretic measures to characterize biological sequences and data.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Union

import numpy as np

from metainformant.core import logging
from . import syntactic, semantic, continuous, estimation

logger = logging.get_logger(__name__)


def information_profile(
    sequences: List[str], k: int = 1, normalize: bool = True, method: str = "plugin"
) -> Dict[str, Any]:
    """Calculate information profile for a set of sequences.

    The information profile characterizes the information content
    across different sequence positions using k-mer entropy.

    Args:
        sequences: List of DNA/protein sequences (must be aligned)
        k: k-mer size for analysis
        normalize: Whether to normalize entropy values
        method: Entropy estimation method

    Returns:
        Dictionary with position-wise entropy and profile statistics

    Raises:
        ValueError: If sequences have different lengths or are too short
    """
    if not sequences:
        raise ValueError("Sequence list cannot be empty")

    seq_length = len(sequences[0])
    if not all(len(seq) == seq_length for seq in sequences):
        raise ValueError("All sequences must have the same length (aligned)")

    if seq_length < k:
        raise ValueError(f"Sequence length {seq_length} too short for k={k}")

    n_sequences = len(sequences)
    profile = []

    # Calculate entropy at each position (sliding window)
    for i in range(seq_length - k + 1):
        # Extract k-mer column
        kmers = [seq[i : i + k] for seq in sequences]

        # Count k-mer frequencies
        from collections import Counter

        kmer_counts = Counter(kmers)

        # Estimate entropy
        entropy_val = estimation.entropy_estimator(kmer_counts, method=method, bias_correction=True)

        # Normalize by maximum possible entropy (log2 of alphabet size)
        if normalize:
            # Estimate alphabet size from data
            alphabet_size = len(set("".join(sequences)))
            max_entropy = math.log2(alphabet_size**k)
            if max_entropy > 0:
                entropy_val /= max_entropy

        profile.append(
            {
                "position": i,
                "entropy": entropy_val,
                "k": k,
                "unique_kmers": len(kmer_counts),
                "total_kmers": n_sequences,
            }
        )

    # Calculate profile statistics
    entropy_values = [pos["entropy"] for pos in profile]

    profile_stats = {
        "mean_entropy": float(np.mean(entropy_values)),
        "std_entropy": float(np.std(entropy_values)),
        "min_entropy": float(np.min(entropy_values)),
        "max_entropy": float(np.max(entropy_values)),
        "entropy_range": float(np.ptp(entropy_values)),
        "low_complexity_regions": _identify_low_complexity_regions(profile),
        "high_complexity_regions": _identify_high_complexity_regions(profile),
    }

    result = {
        "entropy": profile_stats["mean_entropy"],
        "positions": [pos["position"] for pos in profile],
        "profile": profile,
        "statistics": profile_stats,
        "parameters": {
            "k": k,
            "normalize": normalize,
            "method": method,
            "n_sequences": n_sequences,
            "sequence_length": seq_length,
        },
    }

    logger.info(
        f"Calculated information profile: {len(profile)} positions, "
        f"mean entropy = {profile_stats['mean_entropy']:.3f}"
    )

    return result


def _identify_low_complexity_regions(
    profile: List[Dict[str, Any]], threshold_percentile: float = 25.0, min_length: int = 3
) -> List[Dict[str, Any]]:
    """Identify regions with low information content."""
    entropy_values = [pos["entropy"] for pos in profile]
    threshold = np.percentile(entropy_values, threshold_percentile)

    regions = []
    current_region = None

    for i, pos_data in enumerate(profile):
        entropy_val = pos_data["entropy"]

        if entropy_val <= threshold:
            if current_region is None:
                current_region = {"start": i, "entropy_values": []}
            current_region["entropy_values"].append(entropy_val)
        else:
            if current_region is not None and len(current_region["entropy_values"]) >= min_length:
                current_region["end"] = i - 1
                current_region["length"] = len(current_region["entropy_values"])
                current_region["mean_entropy"] = np.mean(current_region["entropy_values"])
                regions.append(current_region)
            current_region = None

    # Handle region at end
    if current_region is not None and len(current_region["entropy_values"]) >= min_length:
        current_region["end"] = len(profile) - 1
        current_region["length"] = len(current_region["entropy_values"])
        current_region["mean_entropy"] = np.mean(current_region["entropy_values"])
        regions.append(current_region)

    return regions


def _identify_high_complexity_regions(
    profile: List[Dict[str, Any]], threshold_percentile: float = 75.0, min_length: int = 3
) -> List[Dict[str, Any]]:
    """Identify regions with high information content."""
    entropy_values = [pos["entropy"] for pos in profile]
    threshold = np.percentile(entropy_values, threshold_percentile)

    regions = []
    current_region = None

    for i, pos_data in enumerate(profile):
        entropy_val = pos_data["entropy"]

        if entropy_val >= threshold:
            if current_region is None:
                current_region = {"start": i, "entropy_values": []}
            current_region["entropy_values"].append(entropy_val)
        else:
            if current_region is not None and len(current_region["entropy_values"]) >= min_length:
                current_region["end"] = i - 1
                current_region["length"] = len(current_region["entropy_values"])
                current_region["mean_entropy"] = np.mean(current_region["entropy_values"])
                regions.append(current_region)
            current_region = None

    # Handle region at end
    if current_region is not None and len(current_region["entropy_values"]) >= min_length:
        current_region["end"] = len(profile) - 1
        current_region["length"] = len(current_region["entropy_values"])
        current_region["mean_entropy"] = np.mean(current_region["entropy_values"])
        regions.append(current_region)

    return regions


def information_signature(
    data: Union[np.ndarray, List[List[float]]], method: str = "entropy", **kwargs: Any
) -> Dict[str, Any]:
    """Calculate information signature of multivariate data.

    The information signature characterizes the dependence structure
    and information content of multivariate datasets.

    Args:
        data: 2D array (n_samples x n_features)
        method: Analysis method ('entropy', 'mutual_info', 'copula')
        **kwargs: Additional parameters for analysis

    Returns:
        Dictionary with information signature

    Raises:
        ValueError: If data is not 2D or has insufficient samples
    """
    data = np.asarray(data)

    if data.ndim != 2:
        raise ValueError("Data must be 2D (n_samples x n_features)")

    n_samples, n_features = data.shape

    if n_samples < 10:
        raise ValueError("Need at least 10 samples for analysis")
    if n_features < 2:
        raise ValueError("Need at least 2 features for signature analysis")

    signature = {
        "method": method,
        "n_samples": n_samples,
        "n_features": n_features,
        "features": {},
    }

    if method == "entropy":
        # Individual feature entropies
        feature_entropies = []
        for i in range(n_features):
            entropy_val = continuous.differential_entropy(data[:, i], **kwargs)
            feature_entropies.append(entropy_val)
            signature["features"][f"feature_{i}"] = {
                "entropy": entropy_val,
                "index": i,
            }

        signature["summary"] = {
            "mean_entropy": float(np.mean(feature_entropies)),
            "std_entropy": float(np.std(feature_entropies)),
            "total_entropy": float(np.sum(feature_entropies)),
        }

    elif method == "mutual_info":
        # Pairwise mutual information matrix
        mi_matrix = np.zeros((n_features, n_features))

        for i in range(n_features):
            for j in range(i + 1, n_features):
                mi_val = continuous.mutual_information_continuous(data[:, i], data[:, j], **kwargs)
                mi_matrix[i, j] = mi_val
                mi_matrix[j, i] = mi_val

        signature["mutual_information_matrix"] = mi_matrix.tolist()
        signature["summary"] = {
            "mean_mi": float(np.mean(mi_matrix)),
            "max_mi": float(np.max(mi_matrix)),
            "mi_sparsity": float(np.count_nonzero(mi_matrix) / mi_matrix.size),
        }

    elif method == "copula":
        # Copula entropy for dependence analysis
        copula_ent = continuous.copula_entropy(data, **kwargs)

        signature["copula_entropy"] = copula_ent
        signature["summary"] = {
            "copula_entropy": copula_ent,
            "normalized_entropy": copula_ent / math.log2(n_features) if n_features > 1 else 0,
        }

    else:
        raise ValueError(f"Unknown method: {method}")

    logger.info(f"Calculated information signature using {method}: " f"{n_samples} samples × {n_features} features")

    return signature


def analyze_sequence_information(
    sequence: str, k_values: Optional[List[int]] = None, methods: Optional[List[str]] = None
) -> Dict[str, Any]:
    """Comprehensive information analysis of a biological sequence.

    Args:
        sequence: Biological sequence (DNA, RNA, or protein)
        k_values: List of k-mer sizes to analyze
        methods: List of analysis methods

    Returns:
        Dictionary with comprehensive sequence analysis

    Raises:
        ValueError: If sequence is too short
    """
    if len(sequence) < 10:
        raise ValueError("Sequence too short for analysis")

    if k_values is None:
        k_values = [1, 2, 3]

    if methods is None:
        methods = ["entropy", "complexity", "patterns"]

    analysis_result = {
        "sequence_length": len(sequence),
        "sequence_type": _infer_sequence_type(sequence),
        "k_mer_analysis": {},
        "methods": {},
    }

    # K-mer analysis
    for k in k_values:
        if len(sequence) >= k:
            kmers = [sequence[i : i + k] for i in range(len(sequence) - k + 1)]
            from collections import Counter

            kmer_counts = Counter(kmers)

            analysis_result["k_mer_analysis"][f"k{k}"] = {
                "n_unique": len(kmer_counts),
                "n_total": len(kmers),
                "most_common": kmer_counts.most_common(5),
                "entropy": estimation.entropy_estimator(kmer_counts),
            }

    # Method-specific analyses
    for method in methods:
        if method == "entropy":
            # Overall sequence entropy
            chars = list(sequence.upper())
            char_counts = Counter(chars)
            analysis_result["methods"]["entropy"] = {
                "shannon_entropy": estimation.entropy_estimator(char_counts),
                "alphabet_size": len(char_counts),
                "gc_content": _calculate_gc_content(sequence) if _is_nucleotide(sequence) else None,
            }

        elif method == "complexity":
            # Sequence complexity measures
            analysis_result["methods"]["complexity"] = {
                "linguistic_complexity": _linguistic_complexity(sequence),
                "compression_ratio": _compression_complexity(sequence),
            }

        elif method == "patterns":
            # Pattern analysis
            analysis_result["methods"]["patterns"] = {
                "repeats": _find_repeats(sequence),
                "palindromes": _find_palindromes(sequence),
                "low_complexity_regions": _find_low_complexity_regions(sequence),
            }

    logger.info(f"Completed sequence information analysis: {len(sequence)} bp/aa")
    return analysis_result


def _infer_sequence_type(sequence: str) -> str:
    """Infer biological sequence type."""
    sequence = sequence.upper()
    chars = set(sequence)

    # Nucleotide characters
    nucleotides = set("ATCGU")
    # Amino acid characters (simplified)
    amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

    # If ALL characters are valid nucleotides, classify as DNA/RNA
    if chars <= nucleotides:
        if "U" in chars:
            return "RNA"
        else:
            return "DNA"

    nuc_overlap = len(chars & nucleotides)
    aa_overlap = len(chars & amino_acids)

    if nuc_overlap > aa_overlap:
        if "U" in chars:
            return "RNA"
        else:
            return "DNA"
    elif aa_overlap > nuc_overlap:
        return "protein"
    else:
        return "unknown"


def _calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of nucleotide sequence."""
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    total = len(sequence) - sequence.count("N")  # Exclude ambiguous bases

    return gc_count / total if total > 0 else 0.0


def _is_nucleotide(sequence: str) -> bool:
    """Check if sequence appears to be nucleotide."""
    seq_type = _infer_sequence_type(sequence)
    return seq_type in ["DNA", "RNA"]


def _linguistic_complexity(sequence: str) -> float:
    """Calculate linguistic complexity using LZW compression."""
    # Simple LZW-like complexity measure
    dictionary = {}
    next_code = 256
    compressed_length = 0

    w = ""
    for c in sequence:
        wc = w + c
        if wc in dictionary:
            w = wc
        else:
            compressed_length += 1
            dictionary[wc] = next_code
            next_code += 1
            w = c

    if w:
        compressed_length += 1

    # Complexity = compressed_length / original_length
    return compressed_length / len(sequence) if len(sequence) > 0 else 0.0


def _compression_complexity(sequence: str) -> float:
    """Calculate compression-based complexity."""
    # Use run-length encoding as simple compression
    compressed = []
    current_char = sequence[0]
    count = 1

    for char in sequence[1:]:
        if char == current_char:
            count += 1
        else:
            compressed.extend([current_char, str(count)])
            current_char = char
            count = 1

    compressed.extend([current_char, str(count)])
    compressed_str = "".join(compressed)

    # Compression ratio
    return len(compressed_str) / len(sequence)


def _find_repeats(sequence: str, min_length: int = 3) -> List[Dict[str, Any]]:
    """Find repeated subsequences."""
    repeats = []
    seq_len = len(sequence)

    for length in range(min_length, seq_len // 2 + 1):
        for i in range(seq_len - length + 1):
            pattern = sequence[i : i + length]

            # Count occurrences
            count = sequence.count(pattern)
            if count > 1:
                # Check if this repeat is already recorded
                already_recorded = any(r["pattern"] == pattern and r["start"] == i for r in repeats)

                if not already_recorded:
                    repeats.append(
                        {
                            "pattern": pattern,
                            "length": length,
                            "count": count,
                            "start": i,
                        }
                    )

    return repeats


def _find_palindromes(sequence: str, min_length: int = 4) -> List[Dict[str, Any]]:
    """Find palindromic sequences."""
    palindromes = []
    seq_len = len(sequence)

    for i in range(seq_len - min_length + 1):
        for j in range(min_length, seq_len - i + 1):
            substring = sequence[i : i + j]
            if _is_palindrome(substring):
                palindromes.append(
                    {
                        "sequence": substring,
                        "start": i,
                        "length": j,
                    }
                )

    return palindromes


def _is_palindrome(sequence: str) -> bool:
    """Check if sequence is a palindrome."""
    # For DNA/RNA, consider complement
    if _is_nucleotide(sequence):
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "U": "A", "N": "N"}
        complement_seq = "".join(complement.get(c, c) for c in sequence[::-1].upper())
        return sequence.upper() == complement_seq
    else:
        # For proteins or other, simple palindrome check
        return sequence.upper() == sequence[::-1].upper()


def _find_low_complexity_regions(sequence: str, window_size: int = 20, threshold: float = 0.7) -> List[Dict[str, Any]]:
    """Find low complexity regions in sequence."""
    regions = []
    seq_len = len(sequence)

    for i in range(0, seq_len - window_size + 1, window_size // 2):
        window = sequence[i : i + window_size]

        # Calculate complexity (entropy)
        from collections import Counter

        char_counts = Counter(window)
        entropy = estimation.entropy_estimator(char_counts)

        # Maximum possible entropy
        alphabet_size = len(set(window))
        max_entropy = math.log2(alphabet_size) if alphabet_size > 0 else 1

        complexity_ratio = entropy / max_entropy if max_entropy > 0 else 0

        if complexity_ratio < threshold:
            regions.append(
                {
                    "start": i,
                    "end": min(i + window_size, seq_len),
                    "complexity": complexity_ratio,
                    "sequence": window,
                }
            )

    return regions


def compare_sequences_information(seq1: str, seq2: str, k: int = 1, method: str = "mutual_info") -> Dict[str, Any]:
    """Compare information content between two sequences.

    Args:
        seq1: First sequence
        seq2: Second sequence
        k: k-mer size for comparison
        method: Comparison method

    Returns:
        Dictionary with comparison results

    Raises:
        ValueError: If sequences have different lengths
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length for comparison")

    comparison = {
        "sequence_lengths": len(seq1),
        "k": k,
        "method": method,
    }

    if method == "mutual_info":
        # Treat sequences as aligned and calculate positional mutual information
        positional_mi = []

        for i in range(len(seq1) - k + 1):
            chars1 = [seq1[j] for j in range(i, min(i + k, len(seq1)))]
            chars2 = [seq2[j] for j in range(i, min(i + k, len(seq2)))]

            if len(chars1) == len(chars2) == k:
                mi_val = syntactic.mutual_information(chars1, chars2)
                positional_mi.append(mi_val)

        comparison["positional_mutual_info"] = positional_mi
        comparison["mean_mi"] = float(np.mean(positional_mi)) if positional_mi else 0.0
        comparison["max_mi"] = float(np.max(positional_mi)) if positional_mi else 0.0

    elif method == "entropy":
        # Compare entropy profiles
        profile1 = information_profile([seq1], k=k)
        profile2 = information_profile([seq2], k=k)

        comparison["entropy_profile_1"] = profile1
        comparison["entropy_profile_2"] = profile2
        comparison["entropy_difference"] = [
            p1["entropy"] - p2["entropy"] for p1, p2 in zip(profile1["profile"], profile2["profile"])
        ]

    elif method == "kld":
        # Kullback-Leibler divergence between position-wise distributions
        kld_values = []

        for i in range(len(seq1) - k + 1):
            window1 = seq1[i : i + k]
            window2 = seq2[i : i + k]

            # Convert to count distributions
            from collections import Counter

            counts1 = Counter(window1)
            counts2 = Counter(window2)

            # Estimate KL divergence
            kld_val = estimation.kl_divergence_estimator(list(counts1.elements()), list(counts2.elements()))
            kld_values.append(kld_val)

        comparison["positional_kld"] = kld_values
        comparison["mean_kld"] = float(np.mean(kld_values)) if kld_values else 0.0

    logger.info(f"Compared sequences using {method}: {len(seq1)} bp/aa each")
    return comparison
