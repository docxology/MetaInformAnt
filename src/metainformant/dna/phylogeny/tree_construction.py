"""Phylogenetic tree construction algorithms.

This module provides functions for building phylogenetic trees from DNA
sequences using neighbor-joining, UPGMA, and k-mer distance methods,
along with the distance matrix computation helpers they depend on.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any, Dict, List

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Type alias for tree representation
Tree = Dict[str, Any]


def neighbor_joining_tree(id_to_seq: Dict[str, str]) -> Tree:
    """Construct a neighbor-joining phylogenetic tree from sequences.

    Args:
        id_to_seq: Dictionary mapping sequence IDs to DNA sequences

    Returns:
        Tree represented as nested dictionary structure
    """
    if len(id_to_seq) < 2:
        raise ValueError("Need at least 2 sequences for tree construction")

    # Calculate distance matrix
    taxa = list(id_to_seq.keys())
    distance_matrix = _calculate_distance_matrix(id_to_seq)

    # Initialize tree with taxa as leaves
    tree = {taxon: None for taxon in taxa}
    active_taxa = set(taxa)

    while len(active_taxa) > 2:
        # Find closest pair
        min_dist = float("inf")
        closest_pair = None

        active_list = list(active_taxa)
        for i in range(len(active_list)):
            for j in range(i + 1, len(active_list)):
                taxon1, taxon2 = active_list[i], active_list[j]
                idx1, idx2 = taxa.index(taxon1), taxa.index(taxon2)

                if distance_matrix[idx1][idx2] < min_dist:
                    min_dist = distance_matrix[idx1][idx2]
                    closest_pair = (taxon1, taxon2)

        if not closest_pair:
            break

        taxon1, taxon2 = closest_pair
        idx1, idx2 = taxa.index(taxon1), taxa.index(taxon2)

        # Calculate branch lengths
        r1 = sum(
            distance_matrix[idx1][j] for j in range(len(taxa)) if taxa[j] in active_taxa and j != idx1 and j != idx2
        )
        r2 = sum(
            distance_matrix[idx2][j] for j in range(len(taxa)) if taxa[j] in active_taxa and j != idx1 and j != idx2
        )

        n_active = len(active_taxa)
        branch1 = (min_dist + (r1 - r2) / (n_active - 2)) / 2
        branch2 = min_dist - branch1

        # Create new internal node
        new_node = f"Node_{len(tree)}"

        # Update tree structure
        tree[new_node] = {taxon1: branch1, taxon2: branch2}

        # Remove old taxa and add new node
        active_taxa.remove(taxon1)
        active_taxa.remove(taxon2)
        active_taxa.add(new_node)

        # Update distance matrix
        distance_matrix = _update_distance_matrix(distance_matrix, taxa, idx1, idx2, new_node, active_taxa)

        # Update taxa list
        taxa.remove(taxon1)
        taxa.remove(taxon2)
        taxa.append(new_node)

    # Connect final two nodes
    if len(active_taxa) == 2:
        remaining = list(active_taxa)
        final_node = f"Root_{len(tree)}"
        tree[final_node] = {
            remaining[0]: distance_matrix[taxa.index(remaining[0])][taxa.index(remaining[1])] / 2,
            remaining[1]: distance_matrix[taxa.index(remaining[0])][taxa.index(remaining[1])] / 2,
        }

    return tree


def upgma_tree(id_to_seq: Dict[str, str]) -> Tree:
    """Construct a UPGMA (Unweighted Pair Group Method with Arithmetic Mean) tree.

    Args:
        id_to_seq: Dictionary mapping sequence IDs to DNA sequences

    Returns:
        Tree represented as nested dictionary structure
    """
    if len(id_to_seq) < 2:
        raise ValueError("Need at least 2 sequences for tree construction")

    # Calculate distance matrix
    taxa = list(id_to_seq.keys())
    distance_matrix = _calculate_distance_matrix(id_to_seq)

    # Initialize tree
    tree = {taxon: None for taxon in taxa}
    active_taxa = set(taxa)
    cluster_sizes = {taxon: 1 for taxon in taxa}

    while len(active_taxa) > 1:
        # Find closest pair
        min_dist = float("inf")
        closest_pair = None

        active_list = list(active_taxa)
        for i in range(len(active_list)):
            for j in range(i + 1, len(active_list)):
                taxon1, taxon2 = active_list[i], active_list[j]
                idx1, idx2 = taxa.index(taxon1), taxa.index(taxon2)

                if distance_matrix[idx1][idx2] < min_dist:
                    min_dist = distance_matrix[idx1][idx2]
                    closest_pair = (taxon1, taxon2)

        if not closest_pair:
            break

        taxon1, taxon2 = closest_pair
        idx1, idx2 = taxa.index(taxon1), taxa.index(taxon2)

        # Create new cluster
        size1 = cluster_sizes[taxon1]
        size2 = cluster_sizes[taxon2]
        total_size = size1 + size2

        new_node = f"Cluster_{len(tree)}"

        # Calculate branch lengths
        branch1 = min_dist / 2
        branch2 = min_dist / 2

        # Update tree
        tree[new_node] = {taxon1: branch1, taxon2: branch2}

        # Update cluster sizes
        cluster_sizes[new_node] = total_size

        # Remove old taxa and add new node
        active_taxa.remove(taxon1)
        active_taxa.remove(taxon2)
        active_taxa.add(new_node)

        # Update distance matrix for UPGMA
        distance_matrix = _update_distance_matrix_upgma(
            distance_matrix, taxa, idx1, idx2, new_node, active_taxa, size1, size2, total_size
        )

        # Update taxa list
        taxa.remove(taxon1)
        taxa.remove(taxon2)
        taxa.append(new_node)

    return tree


def nj_tree_from_kmer(id_to_seq: Dict[str, str], *, k: int = 3, metric: str = "cosine") -> Tree:
    """Construct neighbor-joining tree using k-mer frequencies.

    Args:
        id_to_seq: Dictionary mapping sequence IDs to DNA sequences
        k: k-mer size
        metric: Distance metric ('cosine', 'euclidean', 'jaccard')

    Returns:
        Neighbor-joining tree based on k-mer distances
    """
    if len(id_to_seq) < 2:
        raise ValueError("Need at least 2 sequences")

    # Calculate k-mer distance matrix
    distance_matrix = _calculate_kmer_distance_matrix(id_to_seq, k, metric)

    # Apply full neighbor-joining algorithm using k-mer distances
    taxa = list(id_to_seq.keys())

    # Initialize tree with taxa as leaves
    tree = {taxon: None for taxon in taxa}
    active_taxa = set(taxa)

    # Work with a copy of distance matrix that we can modify
    working_matrix = [row[:] for row in distance_matrix]
    working_taxa = taxa[:]

    while len(active_taxa) > 2:
        # Find closest pair using NJ criterion
        min_dist = float("inf")
        closest_pair = None

        active_list = list(active_taxa)
        n_active = len(active_list)

        for i in range(len(active_list)):
            for j in range(i + 1, len(active_list)):
                taxon1, taxon2 = active_list[i], active_list[j]
                idx1 = working_taxa.index(taxon1)
                idx2 = working_taxa.index(taxon2)

                # Calculate Q matrix element (NJ criterion)
                r_i = sum(
                    working_matrix[idx1][l]
                    for l in range(len(working_taxa))
                    if working_taxa[l] in active_taxa and l != idx1
                )
                r_j = sum(
                    working_matrix[idx2][l]
                    for l in range(len(working_taxa))
                    if working_taxa[l] in active_taxa and l != idx2
                )

                q_ij = (n_active - 2) * working_matrix[idx1][idx2] - r_i - r_j

                if q_ij < min_dist:
                    min_dist = q_ij
                    closest_pair = (taxon1, taxon2)

        if not closest_pair:
            break

        taxon1, taxon2 = closest_pair
        idx1 = working_taxa.index(taxon1)
        idx2 = working_taxa.index(taxon2)

        # Calculate branch lengths using NJ formula
        r1 = sum(
            working_matrix[idx1][j]
            for j in range(len(working_taxa))
            if working_taxa[j] in active_taxa and j != idx1 and j != idx2
        )
        r2 = sum(
            working_matrix[idx2][j]
            for j in range(len(working_taxa))
            if working_taxa[j] in active_taxa and j != idx1 and j != idx2
        )

        d_ij = working_matrix[idx1][idx2]
        n_remaining = len(active_taxa) - 2

        if n_remaining > 0:
            branch1 = max(0, d_ij / 2 + (r1 - r2) / (2 * n_remaining))
            branch2 = max(0, d_ij - branch1)
        else:
            branch1 = d_ij / 2
            branch2 = d_ij / 2

        # Create new internal node
        new_node = f"KmerNode_{len(tree)}"

        # Update tree structure
        tree[new_node] = {taxon1: branch1, taxon2: branch2}

        # Calculate distances to new node
        new_distances = []
        for k_idx in range(len(working_taxa)):
            if working_taxa[k_idx] == taxon1 or working_taxa[k_idx] == taxon2:
                new_distances.append(0.0)
            else:
                # NJ distance to new node
                d_ik = working_matrix[idx1][k_idx]
                d_jk = working_matrix[idx2][k_idx]
                d_uk = (d_ik + d_jk - d_ij) / 2
                new_distances.append(max(0, d_uk))

        # Update working matrix
        new_row = new_distances + [0.0]
        for row in working_matrix:
            row.append(0.0)
        working_matrix.append(new_row)
        for k_idx in range(len(working_matrix) - 1):
            working_matrix[-1][k_idx] = new_distances[k_idx]
            working_matrix[k_idx][-1] = new_distances[k_idx]

        # Update taxa tracking
        working_taxa.append(new_node)
        active_taxa.remove(taxon1)
        active_taxa.remove(taxon2)
        active_taxa.add(new_node)

    # Connect final two nodes
    if len(active_taxa) == 2:
        remaining = list(active_taxa)
        idx1 = working_taxa.index(remaining[0])
        idx2 = working_taxa.index(remaining[1])
        final_dist = working_matrix[idx1][idx2]

        final_node = f"KmerRoot_{len(tree)}"
        tree[final_node] = {
            remaining[0]: final_dist / 2,
            remaining[1]: final_dist / 2,
        }

    return tree


# ---------------------------------------------------------------------------
# Distance matrix computation helpers
# ---------------------------------------------------------------------------


def _calculate_distance_matrix(id_to_seq: Dict[str, str]) -> List[List[float]]:
    """Calculate pairwise distance matrix from sequences."""
    taxa = list(id_to_seq.keys())
    n = len(taxa)
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            seq1 = id_to_seq[taxa[i]]
            seq2 = id_to_seq[taxa[j]]

            # Use p-distance
            distance = _p_distance(seq1, seq2)
            matrix[i][j] = distance
            matrix[j][i] = distance

    return matrix


def _calculate_kmer_distance_matrix(id_to_seq: Dict[str, str], k: int, metric: str) -> List[List[float]]:
    """Calculate k-mer based distance matrix."""
    taxa = list(id_to_seq.keys())
    n = len(taxa)

    # Calculate k-mer frequency vectors
    kmer_vectors = {}
    for taxon, seq in id_to_seq.items():
        kmer_vectors[taxon] = _kmer_frequencies(seq, k)

    # Calculate distance matrix
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            vec1 = kmer_vectors[taxa[i]]
            vec2 = kmer_vectors[taxa[j]]

            distance = _vector_distance(vec1, vec2, metric)
            matrix[i][j] = distance
            matrix[j][i] = distance

    return matrix


def _p_distance(seq1: str, seq2: str) -> float:
    """Calculate p-distance between two sequences."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be aligned (same length)")

    differences = 0
    valid_sites = 0

    for a, b in zip(seq1.upper(), seq2.upper()):
        if a in "ATCG" and b in "ATCG":
            valid_sites += 1
            if a != b:
                differences += 1

    return differences / valid_sites if valid_sites > 0 else 0.0


def _kmer_frequencies(seq: str, k: int) -> Dict[str, float]:
    """Calculate k-mer frequencies for a sequence."""
    seq_upper = seq.upper()
    kmers = []

    for i in range(len(seq_upper) - k + 1):
        kmer = seq_upper[i : i + k]
        if all(c in "ATCG" for c in kmer):
            kmers.append(kmer)

    counts = Counter(kmers)
    total = sum(counts.values())

    frequencies = {}
    for kmer, count in counts.items():
        frequencies[kmer] = count / total if total > 0 else 0.0

    return frequencies


def _vector_distance(vec1: Dict[str, float], vec2: Dict[str, float], metric: str) -> float:
    """Calculate distance between two frequency vectors."""
    all_keys = set(vec1.keys()) | set(vec2.keys())

    if metric == "cosine":
        # Cosine distance
        dot_product = sum(vec1.get(k, 0) * vec2.get(k, 0) for k in all_keys)
        norm1 = math.sqrt(sum(v**2 for v in vec1.values()))
        norm2 = math.sqrt(sum(v**2 for v in vec2.values()))

        if norm1 == 0 or norm2 == 0:
            return 1.0

        return 1 - (dot_product / (norm1 * norm2))

    elif metric == "euclidean":
        # Euclidean distance
        return math.sqrt(sum((vec1.get(k, 0) - vec2.get(k, 0)) ** 2 for k in all_keys))

    elif metric == "jaccard":
        # Jaccard distance
        intersection = sum(min(vec1.get(k, 0), vec2.get(k, 0)) for k in all_keys)
        union = sum(max(vec1.get(k, 0), vec2.get(k, 0)) for k in all_keys)

        if union == 0:
            return 0.0

        return 1 - (intersection / union)

    else:
        raise ValueError(f"Unknown metric: {metric}")


def _update_distance_matrix(
    matrix: List[List[float]], taxa: List[str], idx1: int, idx2: int, new_node: str, active_taxa: set
) -> List[List[float]]:
    """Update distance matrix after joining two taxa in neighbor-joining.

    Uses the NJ formula: d(u,k) = (d(i,k) + d(j,k) - d(i,j)) / 2
    where u is the new node joining i and j.
    """
    n = len(matrix)

    # Ensure idx1 < idx2 for consistent indexing
    if idx1 > idx2:
        idx1, idx2 = idx2, idx1

    # Create new matrix with one less dimension
    new_n = n - 1
    new_matrix = [[0.0] * new_n for _ in range(new_n)]

    # Map old indices to new indices (excluding idx1 and idx2)
    old_to_new = {}
    new_idx = 0
    for old_idx in range(n):
        if old_idx != idx1 and old_idx != idx2:
            old_to_new[old_idx] = new_idx
            new_idx += 1

    # The last index in new matrix is for the joined node
    joined_idx = new_n - 1

    # Copy distances between remaining taxa
    for i in range(n):
        if i == idx1 or i == idx2:
            continue
        for j in range(i + 1, n):
            if j == idx1 or j == idx2:
                continue
            new_i = old_to_new[i]
            new_j = old_to_new[j]
            new_matrix[new_i][new_j] = matrix[i][j]
            new_matrix[new_j][new_i] = matrix[i][j]

    # Calculate distances from new node to remaining taxa
    d_ij = matrix[idx1][idx2]
    for k in range(n):
        if k == idx1 or k == idx2:
            continue
        new_k = old_to_new[k]
        # NJ distance formula
        d_uk = (matrix[idx1][k] + matrix[idx2][k] - d_ij) / 2.0
        new_matrix[joined_idx][new_k] = d_uk
        new_matrix[new_k][joined_idx] = d_uk

    return new_matrix


def _update_distance_matrix_upgma(
    matrix: List[List[float]],
    taxa: List[str],
    idx1: int,
    idx2: int,
    new_node: str,
    active_taxa: set,
    size1: int,
    size2: int,
    total_size: int,
) -> List[List[float]]:
    """Update distance matrix for UPGMA after joining clusters.

    Uses weighted average: d(u,k) = (size1 * d(i,k) + size2 * d(j,k)) / (size1 + size2)
    """
    n = len(matrix)

    # Ensure idx1 < idx2 for consistent indexing
    if idx1 > idx2:
        idx1, idx2 = idx2, idx1
        size1, size2 = size2, size1

    # Create new matrix with one less dimension
    new_n = n - 1
    new_matrix = [[0.0] * new_n for _ in range(new_n)]

    # Map old indices to new indices (excluding idx1 and idx2)
    old_to_new = {}
    new_idx = 0
    for old_idx in range(n):
        if old_idx != idx1 and old_idx != idx2:
            old_to_new[old_idx] = new_idx
            new_idx += 1

    # The last index in new matrix is for the joined cluster
    joined_idx = new_n - 1

    # Copy distances between remaining taxa
    for i in range(n):
        if i == idx1 or i == idx2:
            continue
        for j in range(i + 1, n):
            if j == idx1 or j == idx2:
                continue
            new_i = old_to_new[i]
            new_j = old_to_new[j]
            new_matrix[new_i][new_j] = matrix[i][j]
            new_matrix[new_j][new_i] = matrix[i][j]

    # Calculate distances from new cluster to remaining taxa (weighted average)
    cluster_size = size1 + size2
    for k in range(n):
        if k == idx1 or k == idx2:
            continue
        new_k = old_to_new[k]
        # UPGMA weighted average distance formula
        d_uk = (size1 * matrix[idx1][k] + size2 * matrix[idx2][k]) / cluster_size
        new_matrix[joined_idx][new_k] = d_uk
        new_matrix[new_k][joined_idx] = d_uk

    return new_matrix
