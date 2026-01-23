"""Evolutionary distance calculations for DNA sequences.

This module provides implementations of various evolutionary distance measures
including Jukes-Cantor, Kimura 2-parameter, p-distance, and distance matrix
construction for phylogenetic analysis.
"""

from __future__ import annotations

import math
from typing import Dict

import numpy as np
import pandas as pd

from metainformant.core import logging

logger = logging.get_logger(__name__)


def jukes_cantor_distance(seq1: str, seq2: str) -> float:
    """Calculate Jukes-Cantor evolutionary distance between two sequences.

    The Jukes-Cantor model assumes equal mutation rates between all nucleotides
    and no transition/transversion bias.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence

    Returns:
        Evolutionary distance (substitutions per site)

    Raises:
        ValueError: If sequences have different lengths or contain invalid nucleotides

    Example:
        >>> dist = jukes_cantor_distance("ATCG", "ATCG")
        >>> dist
        0.0
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have equal length")

    if not seq1 or not seq2:
        return 0.0

    # Calculate proportion of differences (p-distance)
    p = p_distance(seq1, seq2)

    # Jukes-Cantor correction
    if p >= 0.75:  # Maximum possible under Jukes-Cantor model
        return float("inf")  # Infinite distance

    distance = -0.75 * math.log(1 - (4 / 3) * p)
    return distance


def kimura_distance(seq1: str, seq2: str) -> float:
    """Calculate Kimura 2-parameter evolutionary distance.

    The Kimura 2-parameter model accounts for different rates of transitions
    and transversions, providing more accurate distances for sequences with
    transition bias.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence

    Returns:
        Evolutionary distance (substitutions per site)

    Raises:
        ValueError: If sequences have different lengths or contain invalid nucleotides

    Example:
        >>> dist = kimura_distance("ATCG", "ATCG")
        >>> dist
        0.0
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have equal length")

    if not seq1 or not seq2:
        return 0.0

    # Count different types of substitutions
    transitions = 0  # A<->G, C<->T
    transversions = 0  # A<->C, A<->T, G<->C, G<->T
    total_sites = 0

    transition_pairs = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    for a, b in zip(seq1.upper(), seq2.upper()):
        if a not in "ATCG" or b not in "ATCG":
            raise ValueError(f"Invalid nucleotide: {a} or {b}")

        if a != b:
            if (a, b) in transition_pairs:
                transitions += 1
            else:
                transversions += 1

        total_sites += 1

    if total_sites == 0:
        return 0.0

    # Proportions
    P = transitions / total_sites  # transition frequency
    Q = transversions / total_sites  # transversion frequency

    # Kimura 2-parameter distance
    if P + Q >= 1.0:  # Maximum possible
        return float("inf")

    # Avoid division by zero and log of negative values
    if P >= 0.5 or Q >= 0.5:
        return float("inf")

    distance = -0.5 * math.log((1 - 2 * P - Q) * math.sqrt(1 - 2 * Q))
    return distance


def p_distance(seq1: str, seq2: str) -> float:
    """Calculate uncorrected p-distance (proportion of differences).

    This is the simplest distance measure, counting the proportion of
    nucleotide sites that differ between two sequences.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence

    Returns:
        Proportion of differing sites (0.0 to 1.0)

    Raises:
        ValueError: If sequences contain invalid nucleotides

    Example:
        >>> p_distance("ATCG", "ATCG")
        0.0
        >>> p_distance("ATCG", "GCTA")
        1.0
    """
    # Handle different lengths by comparing minimum overlapping region
    min_len = min(len(seq1), len(seq2))
    seq1 = seq1[:min_len]
    seq2 = seq2[:min_len]

    if min_len == 0:
        return 0.0

    differences = 0
    total_sites = 0

    for a, b in zip(seq1.upper(), seq2.upper()):
        if a not in "ATCG" or b not in "ATCG":
            raise ValueError(f"Invalid nucleotide: {a} or {b}")

        if a != b:
            differences += 1
        total_sites += 1

    return differences / total_sites if total_sites > 0 else 0.0


def hamming_distance(seq1: str, seq2: str) -> int:
    """Calculate Hamming distance between two sequences.

    The Hamming distance is the number of positions where the sequences differ.
    Sequences must be of equal length.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence

    Returns:
        Number of differing positions

    Raises:
        ValueError: If sequences have different lengths or contain invalid nucleotides

    Example:
        >>> hamming_distance("ATCG", "ATCG")
        0
        >>> hamming_distance("ATCG", "GCTA")
        4
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have equal length")

    if not seq1 or not seq2:
        return 0

    distance = 0
    for a, b in zip(seq1.upper(), seq2.upper()):
        if a not in "ATCG" or b not in "ATCG":
            raise ValueError(f"Invalid nucleotide: {a} or {b}")
        if a != b:
            distance += 1

    return distance


def distance_matrix(sequences: Dict[str, str], method: str = "jukes_cantor") -> pd.DataFrame:
    """Calculate pairwise distance matrix for a set of sequences.

    Args:
        sequences: Dictionary mapping sequence IDs to DNA sequences
        method: Distance method to use ("jukes_cantor", "kimura", "p_distance", "hamming")

    Returns:
        Distance matrix as pandas DataFrame

    Raises:
        ValueError: If invalid method specified or sequences have inconsistent lengths

    Example:
        >>> seqs = {"seq1": "ATCG", "seq2": "ATCG", "seq3": "GCTA"}
        >>> matrix = distance_matrix(seqs, method="p_distance")
        >>> matrix.shape
        (3, 3)
    """
    if not sequences:
        raise ValueError("Must provide at least one sequence")

    seq_ids = list(sequences.keys())

    # Check that all sequences have the same length
    seq_lengths = {seq_id: len(seq) for seq_id, seq in sequences.items()}
    if len(set(seq_lengths.values())) > 1:
        raise ValueError("All sequences must have the same length")

    # Select distance function
    if method == "jukes_cantor":
        distance_func = jukes_cantor_distance
    elif method == "kimura":
        distance_func = kimura_distance
    elif method == "p_distance":
        distance_func = p_distance
    elif method == "hamming":
        distance_func = hamming_distance
    else:
        raise ValueError(f"Unknown distance method: {method}")

    # Calculate pairwise distances
    n = len(seq_ids)
    distance_matrix_array = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            seq1 = sequences[seq_ids[i]]
            seq2 = sequences[seq_ids[j]]

            try:
                distance = distance_func(seq1, seq2)
                # Handle infinite distances
                if distance == float("inf"):
                    distance = np.nan
            except (ValueError, ZeroDivisionError):
                distance = np.nan

            distance_matrix_array[i, j] = distance
            distance_matrix_array[j, i] = distance  # Symmetric matrix

    # Create pandas DataFrame
    df = pd.DataFrame(distance_matrix_array, index=seq_ids, columns=seq_ids)

    return df


def jc69_distance(seq1: str, seq2: str) -> float:
    """Jukes-Cantor 1969 distance (alias for jukes_cantor_distance).

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence

    Returns:
        Jukes-Cantor distance
    """
    return jukes_cantor_distance(seq1, seq2)


def kimura_2p_distance(seq1: str, seq2: str) -> float:
    """Kimura 2-parameter distance (alias for kimura_distance).

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence

    Returns:
        Kimura 2-parameter distance
    """
    return kimura_distance(seq1, seq2)


def kmer_distance(seq1: str, seq2: str, k: int = 3) -> float:
    """Calculate k-mer based distance between sequences.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence
        k: K-mer size

    Returns:
        K-mer distance (cosine distance between k-mer frequency vectors)
    """
    from collections import Counter
    import math

    def get_kmer_counts(seq: str, k: int) -> Counter:
        """Get k-mer counts for a sequence."""
        return Counter(seq[i : i + k] for i in range(len(seq) - k + 1))

    if len(seq1) < k or len(seq2) < k:
        return 1.0  # Maximum distance for very short sequences

    counts1 = get_kmer_counts(seq1.upper(), k)
    counts2 = get_kmer_counts(seq2.upper(), k)

    # Get all unique k-mers
    all_kmers = set(counts1.keys()) | set(counts2.keys())

    # Calculate cosine similarity
    dot_product = sum(counts1.get(kmer, 0) * counts2.get(kmer, 0) for kmer in all_kmers)
    norm1 = math.sqrt(sum(count**2 for count in counts1.values()))
    norm2 = math.sqrt(sum(count**2 for count in counts2.values()))

    if norm1 == 0 or norm2 == 0:
        return 1.0

    cosine_similarity = dot_product / (norm1 * norm2)

    # Return distance (1 - similarity)
    return 1.0 - cosine_similarity


def tajima_nei_distance(seq1: str, seq2: str, kappa: float = 2.0) -> float:
    """Alias for tamura_nei_distance (Tajima-Nei distance model).

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence
        kappa: Transition/transversion rate ratio

    Returns:
        Tajima-Nei distance
    """
    return tamura_nei_distance(seq1, seq2, kappa)


def tamura_nei_distance(seq1: str, seq2: str, kappa: float = 2.0) -> float:
    """Calculate Tamura-Nei evolutionary distance.

    The Tamura-Nei model extends the Kimura 2-parameter model by allowing
    different rates for different types of transversions.

    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence
        kappa: Transition/transversion rate ratio (default: 2.0)

    Returns:
        Evolutionary distance (substitutions per site)

    Raises:
        ValueError: If sequences have different lengths or contain invalid nucleotides
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have equal length")

    if not seq1 or not seq2:
        return 0.0

    # Count different types of substitutions
    # Transitions: A<->G, C<->T
    # Transversions: A<->C, A<->T, G<->C, G<->T
    transitions = 0
    transversions = 0
    total_sites = 0

    transition_pairs = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    # Count GC content for transversion rate calculation
    gc_count = 0

    for a, b in zip(seq1.upper(), seq2.upper()):
        if a not in "ATCG" or b not in "ATCG":
            raise ValueError(f"Invalid nucleotide: {a} or {b}")

        if a == "G" or a == "C":
            gc_count += 1

        if a != b:
            if (a, b) in transition_pairs:
                transitions += 1
            else:
                transversions += 1

        total_sites += 1

    if total_sites == 0:
        return 0.0

    # Proportions
    P = transitions / total_sites  # transition frequency
    Q = transversions / total_sites  # transversion frequency

    # GC content proportion
    gc_prop = gc_count / total_sites

    # Tamura-Nei distance calculation
    # This is a simplified version; full implementation would require
    # more complex parameter estimation

    if P + Q >= 1.0:
        return float("inf")

    # Simplified Tamura-Nei (using kappa for transition bias)
    try:
        term1 = -kappa * math.log(1 - P / kappa - Q)
        term2 = -(1 / kappa) * math.log(1 - Q)

        if term1 < 0 or term2 < 0:
            return float("inf")

        distance = term1 + term2
        return distance
    except (ValueError, ZeroDivisionError):
        return float("inf")


def kmer_distance_matrix(sequences: Dict[str, str], k: int = 3) -> pd.DataFrame:
    """Calculate k-mer distance matrix for a set of sequences.

    Args:
        sequences: Dictionary mapping sequence IDs to DNA sequences
        k: k-mer size

    Returns:
        Distance matrix as pandas DataFrame

    Example:
        >>> seqs = {"seq1": "ATCG", "seq2": "ATCG", "seq3": "GCTA"}
        >>> matrix = kmer_distance_matrix(seqs, k=2)
        >>> matrix.shape
        (3, 3)
    """
    if not sequences:
        raise ValueError("Must provide at least one sequence")

    seq_ids = list(sequences.keys())

    # Calculate pairwise k-mer distances
    n = len(seq_ids)
    distance_matrix_array = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            seq1 = sequences[seq_ids[i]]
            seq2 = sequences[seq_ids[j]]

            try:
                distance = kmer_distance(seq1, seq2, k)
                # Handle infinite distances
                if distance == float("inf"):
                    distance = np.nan
            except (ValueError, ZeroDivisionError):
                distance = np.nan

            distance_matrix_array[i, j] = distance
            distance_matrix_array[j, i] = distance  # Symmetric matrix

    # Create pandas DataFrame
    df = pd.DataFrame(distance_matrix_array, index=seq_ids, columns=seq_ids)

    return df


def sequence_identity_matrix(sequences: Dict[str, str]) -> pd.DataFrame:
    """Calculate sequence identity matrix (1 - p_distance).

    Args:
        sequences: Dictionary mapping sequence IDs to DNA sequences

    Returns:
        Identity matrix as pandas DataFrame (values from 0.0 to 1.0)

    Example:
        >>> seqs = {"seq1": "ATCG", "seq2": "ATCG", "seq3": "GCTA"}
        >>> matrix = sequence_identity_matrix(seqs)
        >>> matrix.shape
        (3, 3)
    """
    if not sequences:
        raise ValueError("Must provide at least one sequence")

    seq_ids = list(sequences.keys())

    # Calculate pairwise identities
    n = len(seq_ids)
    identity_matrix_array = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            seq1 = sequences[seq_ids[i]]
            seq2 = sequences[seq_ids[j]]

            try:
                distance = p_distance(seq1, seq2)
                identity = 1.0 - distance
            except (ValueError, ZeroDivisionError):
                identity = np.nan

            identity_matrix_array[i, j] = identity
            identity_matrix_array[j, i] = identity  # Symmetric matrix

    # Diagonal should be 1.0 (identity with self)
    for i in range(n):
        identity_matrix_array[i, i] = 1.0

    # Create pandas DataFrame
    df = pd.DataFrame(identity_matrix_array, index=seq_ids, columns=seq_ids)

    return df
