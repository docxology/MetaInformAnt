"""Phylogenetic analysis utilities for DNA sequences.

This module re-exports all public symbols from the split submodules
``tree_construction`` and ``tree_analysis`` so that existing import paths
continue to work unchanged.
"""

from __future__ import annotations

# Re-export the Tree type alias
from .tree_construction import Tree

# Re-export construction functions
from .tree_construction import (
    neighbor_joining_tree,
    nj_tree_from_kmer,
    upgma_tree,
)

# Re-export analysis functions
from .tree_analysis import (
    _find_root,
    _get_all_leaves,
    _get_bipartitions,
    basic_tree_stats,
    bootstrap_support,
    from_newick,
    is_monophyletic,
    prune_tree,
    robinson_foulds_distance,
    to_ascii,
    to_newick,
    total_branch_length,
    tree_diameter,
)

__all__ = [
    "Tree",
    # Construction
    "neighbor_joining_tree",
    "upgma_tree",
    "nj_tree_from_kmer",
    # Analysis
    "to_newick",
    "from_newick",
    "bootstrap_support",
    "to_ascii",
    "basic_tree_stats",
    "robinson_foulds_distance",
    "is_monophyletic",
    "total_branch_length",
    "prune_tree",
    "tree_diameter",
]
