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
                if child == "bootstrap":
                    continue
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


def bootstrap_support(tree: Tree, sequences: Dict[str, str], n_replicates: int = 100, method: str = "nj") -> Tree:
    """Calculate bootstrap support for tree branches.

    Performs bootstrap analysis by resampling alignment sites with replacement,
    building trees from resampled data, and counting how often each clade appears.

    Args:
        tree: Original tree
        sequences: Original sequence data
        n_replicates: Number of bootstrap replicates
        method: Tree building method ('nj' or 'upgma')

    Returns:
        Tree with bootstrap support values (0-100)
    """
    logger.info(f"Calculating bootstrap support with {n_replicates} replicates")

    if not sequences:
        return tree.copy()

    # Get alignment length
    seq_list = list(sequences.values())
    if not seq_list:
        return tree.copy()
    alignment_length = len(seq_list[0])

    if alignment_length == 0:
        return tree.copy()

    # Extract clades from original tree
    original_clades = _extract_clades(tree)

    # Run bootstrap replicates
    clade_counts = {frozenset(clade): 0 for clade in original_clades}

    np.random.seed(42)  # Reproducibility

    for rep in range(n_replicates):
        # Resample sites with replacement
        sampled_sites = np.random.randint(0, alignment_length, alignment_length)

        # Build resampled sequences
        resampled_seqs = {}
        for name, seq in sequences.items():
            resampled_seqs[name] = ''.join(seq[i] if i < len(seq) else '-' for i in sampled_sites)

        # Build tree from resampled data
        try:
            # Build tree directly from sequences (functions handle distance matrix internally)
            if method == "nj":
                rep_tree = neighbor_joining_tree(resampled_seqs)
            else:
                rep_tree = upgma_tree(resampled_seqs)

            # Extract clades from replicate tree
            rep_clades = _extract_clades(rep_tree)

            # Count matching clades
            for clade in original_clades:
                clade_set = frozenset(clade)
                if clade_set in [frozenset(c) for c in rep_clades]:
                    clade_counts[clade_set] = clade_counts.get(clade_set, 0) + 1

        except Exception as e:
            logger.debug(f"Bootstrap replicate {rep} failed: {e}")
            continue

    # Add bootstrap values to tree
    supported_tree = tree.copy()

    def _add_bootstrap_values(node, current_clade=None):
        if isinstance(supported_tree.get(node), dict):
            children = [c for c in supported_tree[node].keys() if c != "bootstrap"]

            # Calculate clade for this node
            clade = _get_descendant_leaves(supported_tree, node)
            clade_set = frozenset(clade)

            # Get bootstrap support
            if clade_set in clade_counts:
                support = int(100 * clade_counts[clade_set] / max(1, n_replicates))
            else:
                support = 0

            supported_tree[node]["bootstrap"] = support

            for child in children:
                _add_bootstrap_values(child)

    for root in supported_tree.keys():
        if supported_tree.get(root) is not None:
            _add_bootstrap_values(root)
            break

    return supported_tree


def _extract_clades(tree: Tree) -> List[List[str]]:
    """Extract all clades (sets of descendant leaves) from a tree."""
    clades = []

    def _get_leaves(node):
        if tree.get(node) is None or not isinstance(tree.get(node), dict):
            return [node]
        children = [c for c in tree[node].keys() if c != "bootstrap"]
        leaves = []
        for child in children:
            leaves.extend(_get_leaves(child))
        return leaves

    def _traverse(node):
        if isinstance(tree.get(node), dict):
            children = [c for c in tree[node].keys() if c != "bootstrap"]
            if children:
                leaves = _get_leaves(node)
                if len(leaves) > 1:  # Only internal nodes
                    clades.append(leaves)
            for child in children:
                _traverse(child)

    for root in tree.keys():
        _traverse(root)
        break

    return clades


def _get_descendant_leaves(tree: Tree, node: str) -> List[str]:
    """Get all leaf descendants of a node."""
    if tree.get(node) is None or not isinstance(tree.get(node), dict):
        return [node]
    children = [c for c in tree[node].keys() if c != "bootstrap"]
    leaves = []
    for child in children:
        leaves.extend(_get_descendant_leaves(tree, child))
    return leaves


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
            children = [c for c in tree[node].keys() if c != "bootstrap"]
            for i, child in enumerate(children):
                is_last_child = i == len(children) - 1

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
            is_last = i == len(children) - 1
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
            if child == "bootstrap":
                continue
            total += _count_leaves(child)

        return total

    def _count_internal_nodes(node):
        if not isinstance(tree[node], dict):
            return 0

        total = 1  # Count this node
        for child in tree[node].keys():
            if child == "bootstrap":
                continue
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
                r_i = sum(working_matrix[idx1][l] for l in range(len(working_taxa))
                         if working_taxa[l] in active_taxa and l != idx1)
                r_j = sum(working_matrix[idx2][l] for l in range(len(working_taxa))
                         if working_taxa[l] in active_taxa and l != idx2)

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
        r1 = sum(working_matrix[idx1][j] for j in range(len(working_taxa))
                 if working_taxa[j] in active_taxa and j != idx1 and j != idx2)
        r2 = sum(working_matrix[idx2][j] for j in range(len(working_taxa))
                 if working_taxa[j] in active_taxa and j != idx1 and j != idx2)

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
    from collections import Counter

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


def _calculate_tree_height(tree: Tree, node: str) -> int:
    """Calculate the height of a tree from a given node."""
    if not isinstance(tree[node], dict):
        return 0

    max_child_height = 0
    for child in tree[node].keys():
        child_height = _calculate_tree_height(tree, child)
        max_child_height = max(max_child_height, child_height)

    return max_child_height + 1
