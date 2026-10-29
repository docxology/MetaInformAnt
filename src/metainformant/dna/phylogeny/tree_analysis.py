"""Phylogenetic tree analysis, comparison, and I/O utilities.

This module provides functions for tree serialization (Newick format),
bootstrap support, tree comparison (Robinson-Foulds distance), monophyly
testing, pruning, diameter calculation, and basic statistics.
"""

from __future__ import annotations

import copy
from typing import Any, Dict, List, Tuple

import numpy as np

from metainformant.core.utils import logging

from .tree_construction import Tree, neighbor_joining_tree, upgma_tree

logger = logging.get_logger(__name__)


def to_newick(tree: Tree) -> str:
    """Convert tree dictionary to Newick format string.

    Args:
        tree: Tree represented as nested dictionary

    Returns:
        Newick format string
    """

    def _to_newick_recursive(node: str) -> str:
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

    root = _find_root(tree)
    return _to_newick_recursive(root) + ";"


_BOOTSTRAP_METHODS = ("nj", "upgma")


def _support_key(clade: frozenset[str], all_leaves: frozenset[str], method: str) -> Any:
    """Identify what a replicate must reproduce for the node whose descendant leaves are ``clade``.

    UPGMA trees are rooted, so the key is the clade itself. Neighbor joining returns an unrooted tree drawn
    with an arbitrary root, so the key is the bipartition the node's branch induces (the smaller side, as in
    ``_get_bipartitions``), which does not depend on where a replicate happens to be rooted. The key is None for a
    trivial split (one side empty or a single leaf): every tree over these taxa has it.
    """
    if method == "upgma":
        return clade
    complement = all_leaves - clade
    if len(clade) <= 1 or len(complement) <= 1:
        return None
    if len(clade) != len(complement):
        return clade if len(clade) < len(complement) else complement
    return min(clade, complement, key=lambda side: sorted(side))


def bootstrap_support(tree: Tree, sequences: Dict[str, str], n_replicates: int = 100, method: str = "nj") -> Tree:
    """Calculate bootstrap support for tree branches.

    Performs bootstrap analysis by resampling alignment sites with replacement,
    building trees from resampled data, and counting how often each branch appears.
    The resampling uses a fixed seed, so repeated calls give the same values.

    A branch of a UPGMA tree is counted when a replicate has the same clade (the same set of descendant leaves
    under the same root). A branch of a neighbor-joining tree is counted when a replicate has the same split of
    the leaves into two groups: neighbor joining does not determine the root, so a replicate that places its root
    elsewhere still supports the branch. Both children of the root, which lie on one branch, therefore get the
    same value, and a branch to a single leaf (or the root itself) is always supported (100).

    Args:
        tree: Original tree, built from ``sequences`` by the same method
        sequences: Original sequence data
        n_replicates: Number of bootstrap replicates
        method: Tree building method ('nj' or 'upgma')

    Returns:
        A new tree with bootstrap support values (0-100) on its internal nodes. The input
        ``tree`` is never modified: its internal-node dictionaries are copied, not shared.

    Raises:
        ValueError: If ``method`` is not 'nj' or 'upgma', or the leaves of ``tree`` are not the names in
            ``sequences``.
    """
    if method not in _BOOTSTRAP_METHODS:
        raise ValueError(f"Unknown bootstrap method {method!r}; expected one of {list(_BOOTSTRAP_METHODS)}")

    logger.info(f"Calculating bootstrap support with {n_replicates} replicates")

    if not sequences:
        return copy.deepcopy(tree)

    # Get alignment length
    seq_list = list(sequences.values())
    if not seq_list:
        return copy.deepcopy(tree)
    alignment_length = len(seq_list[0])

    if alignment_length == 0:
        return copy.deepcopy(tree)

    all_leaves = frozenset(_get_all_leaves(tree))
    if all_leaves != set(sequences):
        only_in_tree = sorted(all_leaves - set(sequences))
        only_in_sequences = sorted(set(sequences) - all_leaves)
        raise ValueError(
            f"The tree and the sequences must have the same taxa: only in tree {only_in_tree}, "
            f"only in sequences {only_in_sequences}"
        )

    # What each internal node of the original tree needs a replicate to reproduce
    node_keys = {
        node: _support_key(frozenset(_get_descendant_leaves(tree, node)), all_leaves, method)
        for node, entry in tree.items()
        if isinstance(entry, dict)
    }
    key_counts = {key: 0 for key in node_keys.values() if key is not None}

    rng = np.random.default_rng(42)  # Reproducibility (local generator)

    for rep in range(n_replicates):
        # Resample sites with replacement
        sampled_sites = rng.integers(0, alignment_length, alignment_length)

        # Build resampled sequences
        resampled_seqs = {}
        for name, seq in sequences.items():
            resampled_seqs[name] = "".join(seq[i] if i < len(seq) else "-" for i in sampled_sites)

        # Build tree from resampled data
        try:
            # Build tree directly from sequences (functions handle distance matrix internally)
            if method == "nj":
                rep_tree = neighbor_joining_tree(resampled_seqs)
                present = _get_bipartitions(rep_tree, set(all_leaves))
            else:
                rep_tree = upgma_tree(resampled_seqs)
                present = {frozenset(clade) for clade in _extract_clades(rep_tree)}
        except Exception as e:
            logger.debug(f"Bootstrap replicate {rep} failed: {e}")
            continue

        for key in key_counts:
            if key in present:
                key_counts[key] += 1

    # Add bootstrap values to a deep copy: a shallow dict.copy() would share the internal-node
    # dictionaries with the caller's tree and write the support values into it.
    supported_tree = copy.deepcopy(tree)

    for node, key in node_keys.items():
        count = n_replicates if key is None else key_counts[key]
        supported_tree[node]["bootstrap"] = int(100 * count / max(1, n_replicates))

    return supported_tree


def to_ascii(tree: Tree) -> str:
    """Convert tree to ASCII art representation.

    The root is printed on the first line and every other node once, on its own line, indented under its
    parent; leaves are not repeated. A ``"bootstrap"`` support entry is not a node and is not printed.

    Args:
        tree: Tree represented as nested dictionary

    Returns:
        ASCII art string representation
    """

    def _build_ascii(node: str, prefix: str = "") -> list[str]:
        lines: list[str] = []
        entry = tree[node]
        children = [c for c in entry.keys() if c != "bootstrap"] if isinstance(entry, dict) else []
        for i, child in enumerate(children):
            is_last_child = i == len(children) - 1

            # Branch symbol
            branch = "└── " if is_last_child else "├── "

            # Child branch
            child_prefix = prefix + ("    " if is_last_child else "│   ")

            lines.append(f"{prefix}{branch}{child}")
            lines.extend(_build_ascii(child, child_prefix))

        return lines

    root = _find_root(tree)

    return "\n".join([root, *_build_ascii(root)])


def basic_tree_stats(tree: Tree) -> Dict[str, int]:
    """Calculate basic statistics for a phylogenetic tree.

    Args:
        tree: Tree represented as nested dictionary

    Returns:
        Dictionary with tree statistics
    """

    def _count_leaves(node: str) -> int:
        if not isinstance(tree[node], dict):
            return 1

        total = 0
        for child in tree[node].keys():
            if child == "bootstrap":
                continue
            total += _count_leaves(child)

        return total

    def _count_internal_nodes(node: str) -> int:
        if not isinstance(tree[node], dict):
            return 0

        total = 1  # Count this node
        for child in tree[node].keys():
            if child == "bootstrap":
                continue
            total += _count_internal_nodes(child)

        return total

    root = _find_root(tree)

    stats = {
        "leaves": _count_leaves(root),
        "internal_nodes": _count_internal_nodes(root),
        "total_nodes": len(tree),
        "tree_height": _calculate_tree_height(tree, root),
    }

    return stats


def robinson_foulds_distance(tree1: Tree, tree2: Tree) -> int:
    """Calculate the Robinson-Foulds distance between two phylogenetic trees.

    The Robinson-Foulds (RF) distance is the symmetric difference of the sets
    of bipartitions (splits) induced by the two trees. Each internal edge in an
    unrooted tree defines a bipartition of the leaf set; the RF distance counts
    how many bipartitions appear in one tree but not the other.

    Args:
        tree1: First tree represented as nested dictionary.
        tree2: Second tree represented as nested dictionary.

    Returns:
        Robinson-Foulds distance (non-negative integer).

    Raises:
        ValueError: If the two trees do not share the same leaf set.
    """
    leaves1 = set(_get_all_leaves(tree1))
    leaves2 = set(_get_all_leaves(tree2))

    if leaves1 != leaves2:
        raise ValueError(
            f"Trees must have the same leaf set. " f"tree1 leaves: {sorted(leaves1)}, tree2 leaves: {sorted(leaves2)}"
        )

    # Extract bipartitions as frozensets of the smaller side
    splits1 = _get_bipartitions(tree1, leaves1)
    splits2 = _get_bipartitions(tree2, leaves2)

    # Symmetric difference
    return len(splits1.symmetric_difference(splits2))


def is_monophyletic(tree: Tree, taxa: List[str]) -> bool:
    """Check if a set of taxa forms a monophyletic group in the tree.

    A set of taxa is monophyletic if there exists an internal node in the tree
    whose set of descendant leaves is exactly the given set of taxa.

    Args:
        tree: Tree represented as nested dictionary.
        taxa: List of taxon names to test for monophyly.

    Returns:
        True if the taxa form a monophyletic group, False otherwise.
    """
    if not taxa:
        return True

    target = set(taxa)
    all_leaves = set(_get_all_leaves(tree))

    # If the target is the full leaf set, it is trivially monophyletic (the root)
    if target == all_leaves:
        return True

    # Check if any node's descendant leaf set matches the target exactly
    for node in tree:
        if isinstance(tree.get(node), dict):
            descendant_leaves = set(_get_descendant_leaves(tree, node))
            if descendant_leaves == target:
                return True

    return False


def total_branch_length(tree: Tree) -> float:
    """Calculate the sum of all branch lengths in the tree.

    Iterates over every internal node and sums the branch lengths of all
    edges. Branches with no length (None) are treated as zero.

    Args:
        tree: Tree represented as nested dictionary.

    Returns:
        Total branch length as a float.
    """
    total = 0.0

    for node, children in tree.items():
        if isinstance(children, dict):
            for child, branch_length in children.items():
                if child == "bootstrap":
                    continue
                if branch_length is not None:
                    total += float(branch_length)

    return total


def from_newick(newick_str: str) -> Tree:
    """Parse a Newick format string into the Tree dictionary structure.

    Supports branch lengths (e.g. ``(A:0.1,B:0.2):0.3;``), nested clades,
    and unlabeled internal nodes (which receive auto-generated names).

    Args:
        newick_str: Newick format string (must end with ``;``).

    Returns:
        Tree represented as nested dictionary consistent with this module.

    Raises:
        ValueError: If the Newick string is empty or malformed.
    """
    from metainformant.core.utils.newick import from_newick as _from_newick

    return _from_newick(newick_str)


def prune_tree(tree: Tree, taxa_to_keep: List[str]) -> Tree:
    """Prune a tree to keep only the specified taxa.

    Removes all leaves not in ``taxa_to_keep`` and collapses any resulting
    internal nodes that have only a single child by merging their branch
    lengths.

    Args:
        tree: Tree represented as nested dictionary.
        taxa_to_keep: List of leaf taxon names to retain.

    Returns:
        New pruned tree containing only the specified taxa.

    Raises:
        ValueError: If ``taxa_to_keep`` is empty or contains names not in the tree.
    """
    if not taxa_to_keep:
        raise ValueError("taxa_to_keep must not be empty")

    keep_set = set(taxa_to_keep)
    all_leaves = set(_get_all_leaves(tree))
    missing = keep_set - all_leaves
    if missing:
        raise ValueError(f"Taxa not found in tree: {sorted(missing)}")

    # Find root
    root = _find_root(tree)

    # Recursively build pruned tree
    pruned: Tree = {}

    def _prune(node: str) -> Tuple[str | None, float]:
        """Return (pruned_node_name, branch_length_to_parent) or (None, 0) if pruned away."""
        if tree.get(node) is None or not isinstance(tree.get(node), dict):
            # Leaf node
            if node in keep_set:
                pruned[node] = None
                return node, 0.0
            else:
                return None, 0.0

        # Internal node: recurse into children
        children = [c for c in tree[node].keys() if c != "bootstrap"]
        surviving_children: Dict[str, float] = {}

        for child in children:
            child_bl = tree[node][child] if tree[node][child] is not None else 0.0
            result_name, extra_bl = _prune(child)
            if result_name is not None:
                surviving_children[result_name] = child_bl + extra_bl

        if len(surviving_children) == 0:
            # No surviving descendants
            return None, 0.0
        elif len(surviving_children) == 1:
            # Collapse: pass through the single child, accumulating branch length
            child_name = next(iter(surviving_children))
            child_bl = surviving_children[child_name]
            return child_name, child_bl
        else:
            # Multiple surviving children: keep this internal node
            pruned[node] = surviving_children
            return node, 0.0

    _prune(root)
    return pruned


def tree_diameter(tree: Tree) -> float:
    """Calculate the diameter of the tree (longest path between any two leaves).

    The diameter is the maximum sum of branch lengths along the path connecting
    any pair of leaf nodes in the tree.

    Args:
        tree: Tree represented as nested dictionary.

    Returns:
        Tree diameter as a float.
    """
    root = _find_root(tree)

    # For each node, compute the farthest leaf distance and track the global max path
    max_diameter = [0.0]

    def _farthest(node: str) -> float:
        """Return the distance from node to its farthest descendant leaf."""
        if tree.get(node) is None or not isinstance(tree.get(node), dict):
            return 0.0

        children = [c for c in tree[node].keys() if c != "bootstrap"]
        if not children:
            return 0.0

        child_depths: List[float] = []
        for child in children:
            bl = tree[node][child] if tree[node][child] is not None else 0.0
            depth = _farthest(child) + float(bl)
            child_depths.append(depth)

        child_depths.sort(reverse=True)

        # The diameter through this node is the sum of the two longest arms
        if len(child_depths) >= 2:
            candidate = child_depths[0] + child_depths[1]
            if candidate > max_diameter[0]:
                max_diameter[0] = candidate

        # Also check single-arm (in case root is a leaf edge)
        if child_depths[0] > max_diameter[0] and len(child_depths) == 1:
            max_diameter[0] = child_depths[0]

        return child_depths[0]

    _farthest(root)
    return max_diameter[0]


# ---------------------------------------------------------------------------
# Internal helper functions
# ---------------------------------------------------------------------------


def _find_root(tree: Tree) -> str:
    """Find the root node of a tree (node with no parent).

    Args:
        tree: Tree represented as nested dictionary.

    Returns:
        Name of the root node.
    """
    all_nodes = set(tree.keys())
    child_nodes: set[str] = set()

    for node_data in tree.values():
        if isinstance(node_data, dict):
            child_nodes.update(node_data.keys())

    roots = all_nodes - child_nodes
    return roots.pop()


def _get_all_leaves(tree: Tree) -> List[str]:
    """Get all leaf nodes from a tree.

    Args:
        tree: Tree represented as nested dictionary.

    Returns:
        List of leaf node names.
    """
    leaves: List[str] = []
    for node, children in tree.items():
        if children is None or not isinstance(children, dict):
            leaves.append(node)
    return leaves


def _get_descendant_leaves(tree: Tree, node: str) -> List[str]:
    """Get all leaf descendants of a node."""
    if tree.get(node) is None or not isinstance(tree.get(node), dict):
        return [node]
    children = [c for c in tree[node].keys() if c != "bootstrap"]
    leaves = []
    for child in children:
        leaves.extend(_get_descendant_leaves(tree, child))
    return leaves


def _extract_clades(tree: Tree) -> List[List[str]]:
    """Extract all clades (sets of descendant leaves) from a tree."""
    clades = []

    def _get_leaves(node: str) -> list[str]:
        if tree.get(node) is None or not isinstance(tree.get(node), dict):
            return [node]
        children = [c for c in tree[node].keys() if c != "bootstrap"]
        leaves = []
        for child in children:
            leaves.extend(_get_leaves(child))
        return leaves

    def _traverse(node: str) -> None:
        if isinstance(tree.get(node), dict):
            children = [c for c in tree[node].keys() if c != "bootstrap"]
            if children:
                leaves = _get_leaves(node)
                if len(leaves) > 1:  # Only internal nodes
                    clades.append(leaves)
            for child in children:
                _traverse(child)

    _traverse(_find_root(tree))

    return clades


def _get_bipartitions(tree: Tree, all_leaves: set[str]) -> set[frozenset[str]]:
    """Extract the set of non-trivial bipartitions from a tree.

    Each internal edge induces a split of the leaf set into two groups.
    We represent each split by the frozenset of the smaller side (or
    alphabetically first if equal size) to ensure canonical form.

    Args:
        tree: Tree represented as nested dictionary.
        all_leaves: Complete set of leaf names.

    Returns:
        Set of frozensets, each representing one side of a bipartition.
    """
    splits: set[frozenset[str]] = set()

    for node in tree:
        if isinstance(tree.get(node), dict):
            children = [c for c in tree[node].keys() if c != "bootstrap"]
            for child in children:
                descendant_leaves = frozenset(_get_descendant_leaves(tree, child))
                complement = frozenset(all_leaves - descendant_leaves)

                # Skip trivial splits (single leaf or entire tree)
                if len(descendant_leaves) <= 1 or len(complement) <= 1:
                    continue

                # Canonical form: use the smaller side, or sorted-first if equal
                if len(descendant_leaves) < len(complement):
                    canonical = descendant_leaves
                elif len(complement) < len(descendant_leaves):
                    canonical = complement
                else:
                    canonical = min(descendant_leaves, complement, key=lambda s: sorted(s))

                splits.add(canonical)

    return splits


def _calculate_tree_height(tree: Tree, node: str) -> int:
    """Calculate the height of a tree from a given node."""
    if not isinstance(tree[node], dict):
        return 0

    max_child_height = 0
    for child in tree[node].keys():
        if child == "bootstrap":
            continue  # a support value stored in the node, not a child node
        child_height = _calculate_tree_height(tree, child)
        max_child_height = max(max_child_height, child_height)

    return max_child_height + 1
