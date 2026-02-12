"""Phylogenetic tree analysis, comparison, and I/O utilities.

This module provides functions for tree serialization (Newick format),
bootstrap support, tree comparison (Robinson-Foulds distance), monophyly
testing, pruning, diameter calculation, and basic statistics.
"""

from __future__ import annotations

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
            resampled_seqs[name] = "".join(seq[i] if i < len(seq) else "-" for i in sampled_sites)

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
    if not newick_str or not newick_str.strip():
        raise ValueError("Newick string is empty")

    newick_str = newick_str.strip()
    if newick_str.endswith(";"):
        newick_str = newick_str[:-1]

    if not newick_str:
        raise ValueError("Newick string is empty after removing semicolon")

    tree: Tree = {}
    _internal_counter = [0]

    def _parse(s: str) -> Tuple[str, float | None]:
        """Parse a Newick sub-expression, returning (node_name, branch_length)."""
        s = s.strip()

        if s.startswith("("):
            # Find matching closing parenthesis
            depth = 0
            end_paren = -1
            for i, ch in enumerate(s):
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
                    if depth == 0:
                        end_paren = i
                        break

            if end_paren == -1:
                raise ValueError("Unmatched parenthesis in Newick string")

            # Content inside parentheses
            inner = s[1:end_paren]

            # Remainder after closing paren: optional label and branch length
            remainder = s[end_paren + 1 :]

            # Parse label and branch length from remainder
            node_label, branch_length = _parse_label_length(remainder)

            if not node_label:
                node_label = f"Internal_{_internal_counter[0]}"
                _internal_counter[0] += 1

            # Split inner by commas at depth 0
            children_strs = _split_at_top_level(inner)

            # Parse each child
            children_dict: Dict[str, float] = {}
            for child_str in children_strs:
                child_name, child_bl = _parse(child_str)
                children_dict[child_name] = child_bl if child_bl is not None else 0.0

            tree[node_label] = children_dict
            return node_label, branch_length

        else:
            # Leaf node: "name:length" or just "name"
            label, branch_length = _parse_label_length(s)
            if not label:
                raise ValueError(f"Empty leaf label in Newick string: '{s}'")
            tree[label] = None
            return label, branch_length

    def _parse_label_length(s: str) -> Tuple[str, float | None]:
        """Parse 'label:length' returning (label, length)."""
        s = s.strip()
        if ":" in s:
            parts = s.rsplit(":", 1)
            label = parts[0].strip()
            try:
                bl = float(parts[1].strip())
            except ValueError:
                bl = None
            return label, bl
        return s, None

    def _split_at_top_level(s: str) -> List[str]:
        """Split string by commas not inside parentheses."""
        parts: List[str] = []
        depth = 0
        current: List[str] = []
        for ch in s:
            if ch == "(":
                depth += 1
                current.append(ch)
            elif ch == ")":
                depth -= 1
                current.append(ch)
            elif ch == "," and depth == 0:
                parts.append("".join(current))
                current = []
            else:
                current.append(ch)
        if current:
            parts.append("".join(current))
        return parts

    _parse(newick_str)
    return tree


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
        child_height = _calculate_tree_height(tree, child)
        max_child_height = max(max_child_height, child_height)

    return max_child_height + 1
