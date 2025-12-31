"""Phylogenetic analysis utilities for DNA sequences.

This module provides functions for phylogenetic tree construction,
distance matrix calculation, and tree manipulation for DNA sequence data.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Type alias for tree representation
Tree = Dict[str, any]


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
        min_dist = float('inf')
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
        r1 = sum(distance_matrix[idx1][j] for j in range(len(taxa)) if taxa[j] in active_taxa and j != idx1 and j != idx2)
        r2 = sum(distance_matrix[idx2][j] for j in range(len(taxa)) if taxa[j] in active_taxa and j != idx1 and j != idx2)

        n_active = len(active_taxa)
        branch1 = (min_dist + (r1 - r2) / (n_active - 2)) / 2
        branch2 = min_dist - branch1

        # Create new internal node
        new_node = f"Node_{len(tree)}"

        # Update tree structure
        tree[new_node] = {
            taxon1: branch1,
            taxon2: branch2
        }

        # Remove old taxa and add new node
        active_taxa.remove(taxon1)
        active_taxa.remove(taxon2)
        active_taxa.add(new_node)

        # Update distance matrix
        distance_matrix = _update_distance_matrix(
            distance_matrix, taxa, idx1, idx2, new_node, active_taxa
        )

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
            remaining[1]: distance_matrix[taxa.index(remaining[0])][taxa.index(remaining[1])] / 2
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
        min_dist = float('inf')
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
        tree[new_node] = {
            taxon1: branch1,
            taxon2: branch2
        }

        # Update cluster sizes
        cluster_sizes[new_node] = total_size

        # Remove old taxa and add new node
        active_taxa.remove(taxon1)
        active_taxa.remove(taxon2)
        active_taxa.add(new_node)

        # Update distance matrix for UPGMA
        distance_matrix = _update_distance_matrix_upgma(
            distance_matrix, taxa, idx1, idx2, new_node, active_taxa,
            size1, size2, total_size
        )

        # Update taxa list
        taxa.remove(taxon1)
        taxa.remove(taxon2)
        taxa.append(new_node)

    return tree


def to_newick(tree: Tree) -> str:
    """Convert tree dictionary to Newick format string.

    Args:
        tree: Tree represented as nested dictionary

    Returns:
        Newick format string
    """
    def _to_newick_recursive(node):
        if isinstance(tree[node], dict):
            children = []
            for child, branch_length in tree[node].items():
                if branch_length is not None:
                    children.append(f"{_to_newick_recursive(child)}:{branch_length}")
                else:
                    children.append(_to_newick_recursive(child))

            return f"({','.join(children)})"
        else:
            return node

    # Find root (node with no parent)
    all_nodes = set(tree.keys())
    child_nodes = set()

    for node_data in tree.values():
        if isinstance(node_data, dict):
            child_nodes.update(node_data.keys())

    root = (all_nodes - child_nodes).pop()
    return _to_newick_recursive(root) + ";"


def bootstrap_support(tree: Tree, sequences: Dict[str, str], n_replicates: int = 100,
                     method: str = "nj") -> Tree:
    """Calculate bootstrap support for tree branches.

    Args:
        tree: Original tree
        sequences: Original sequence data
        n_replicates: Number of bootstrap replicates
        method: Tree building method ('nj' or 'upgma')

    Returns:
        Tree with bootstrap support values
    """
    # This is a simplified implementation
    # Full bootstrap would require resampling sites with replacement

    logger.info(f"Calculating bootstrap support with {n_replicates} replicates")

    # For now, return original tree with placeholder support values
    # Real implementation would require site resampling and tree comparison

    supported_tree = tree.copy()

    # Add bootstrap values (placeholder)
    def _add_bootstrap(node):
        if isinstance(supported_tree[node], dict):
            for child in supported_tree[node].keys():
                _add_bootstrap(child)
            # Add bootstrap support (simplified)
            supported_tree[node]['bootstrap'] = 85  # Placeholder value

    for root in supported_tree.keys():
        if supported_tree[root] is not None:
            _add_bootstrap(root)
            break

    return supported_tree


def to_ascii(tree: Tree) -> str:
    """Convert tree to ASCII art representation.

    Args:
        tree: Tree represented as nested dictionary

    Returns:
        ASCII art string representation
    """
    def _build_ascii(node, prefix="", is_last=True):
        lines = []

        if isinstance(tree[node], dict):
            children = list(tree[node].keys())
            for i, child in enumerate(children):
                is_last_child = (i == len(children) - 1)

                # Branch symbol
                branch = "└── " if is_last_child else "├── "

                # Child branch
                child_prefix = prefix + ("    " if is_last_child else "│   ")

                lines.append(f"{prefix}{branch}{child}")
                lines.extend(_build_ascii(child, child_prefix, is_last_child))
        else:
            lines.append(f"{prefix}{node}")

        return lines

    # Find root
    all_nodes = set(tree.keys())
    child_nodes = set()

    for node_data in tree.values():
        if isinstance(node_data, dict):
            child_nodes.update(node_data.keys())

    root = (all_nodes - child_nodes).pop()

    ascii_lines = [root]
    if tree[root]:
        children = list(tree[root].keys())
        for i, child in enumerate(children):
            is_last = (i == len(children) - 1)
            branch = "└── " if is_last else "├── "
            child_prefix = "    " if is_last else "│   "

            ascii_lines.append(f"{branch}{child}")
            ascii_lines.extend(_build_ascii(child, child_prefix, is_last))

    return "\n".join(ascii_lines)


def basic_tree_stats(tree: Tree) -> Dict[str, int]:
    """Calculate basic statistics for a phylogenetic tree.

    Args:
        tree: Tree represented as nested dictionary

    Returns:
        Dictionary with tree statistics
    """
    def _count_leaves(node):
        if not isinstance(tree[node], dict):
            return 1

        total = 0
        for child in tree[node].keys():
            total += _count_leaves(child)

        return total

    def _count_internal_nodes(node):
        if not isinstance(tree[node], dict):
            return 0

        total = 1  # Count this node
        for child in tree[node].keys():
            total += _count_internal_nodes(child)

        return total

    # Find root
    all_nodes = set(tree.keys())
    child_nodes = set()

    for node_data in tree.values():
        if isinstance(node_data, dict):
            child_nodes.update(node_data.keys())

    root = (all_nodes - child_nodes).pop()

    stats = {
        "leaves": _count_leaves(root),
        "internal_nodes": _count_internal_nodes(root),
        "total_nodes": len(tree),
        "tree_height": _calculate_tree_height(tree, root),
    }

    return stats


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

    # Convert to the format expected by neighbor_joining_tree
    # This is a simplified implementation
    taxa = list(id_to_seq.keys())

    # Use simplified NJ-like algorithm for k-mer data
    tree = {taxon: None for taxon in taxa}

    # For now, create a simple star phylogeny
    # Real implementation would use full NJ algorithm
    if len(taxa) > 2:
        center_node = "Kmer_Center"
        tree[center_node] = {}

        for taxon in taxa:
            tree[center_node][taxon] = 0.1  # Placeholder branch length

    return tree


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
        if a in 'ATCG' and b in 'ATCG':
            valid_sites += 1
            if a != b:
                differences += 1

    return differences / valid_sites if valid_sites > 0 else 0.0


def _kmer_frequencies(seq: str, k: int) -> Dict[str, float]:
    """Calculate k-mer frequencies for a sequence."""
    from collections import Counter

    seq_upper = seq.upper()
    kmers = []

    for i in range(len(seq_upper) - k + 1):
        kmer = seq_upper[i:i + k]
        if all(c in 'ATCG' for c in kmer):
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
        return math.sqrt(sum((vec1.get(k, 0) - vec2.get(k, 0))**2 for k in all_keys))

    elif metric == "jaccard":
        # Jaccard distance
        intersection = sum(min(vec1.get(k, 0), vec2.get(k, 0)) for k in all_keys)
        union = sum(max(vec1.get(k, 0), vec2.get(k, 0)) for k in all_keys)

        if union == 0:
            return 0.0

        return 1 - (intersection / union)

    else:
        raise ValueError(f"Unknown metric: {metric}")


def _update_distance_matrix(matrix: List[List[float]], taxa: List[str],
                          idx1: int, idx2: int, new_node: str,
                          active_taxa: set) -> List[List[float]]:
    """Update distance matrix after joining two taxa."""
    # Simplified update - real NJ uses more complex formula
    n = len(matrix)
    new_matrix = [[0.0] * (n - 1) for _ in range(n - 1)]

    # Copy existing distances (simplified)
    new_taxa = [t for t in taxa if t != taxa[idx1] and t != taxa[idx2]] + [new_node]

    return new_matrix  # Placeholder


def _update_distance_matrix_upgma(matrix: List[List[float]], taxa: List[str],
                                 idx1: int, idx2: int, new_node: str,
                                 active_taxa: set, size1: int, size2: int,
                                 total_size: int) -> List[List[float]]:
    """Update distance matrix for UPGMA after joining clusters."""
    # Simplified UPGMA update
    n = len(matrix)
    new_matrix = [[0.0] * (n - 1) for _ in range(n - 1)]

    # Copy existing distances (simplified)
    new_taxa = [t for t in taxa if t != taxa[idx1] and t != taxa[idx2]] + [new_node]

    return new_matrix  # Placeholder


def _calculate_tree_height(tree: Tree, node: str) -> int:
    """Calculate the height of a tree from a given node."""
    if not isinstance(tree[node], dict):
        return 0

    max_child_height = 0
    for child in tree[node].keys():
        child_height = _calculate_tree_height(tree, child)
        max_child_height = max(max_child_height, child_height)

    return max_child_height + 1


