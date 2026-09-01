"""Newick string parsing.

Domain-neutral tree parsing shared by any clade or datatype module that needs
phylogenetic structure (``dna.phylogeny`` re-exports it for backwards
compatibility; ``rna.analysis`` consumes it directly here rather than
importing across domains).

The parser accepts a Newick string with branch lengths, nested clades, and
unlabeled internal nodes (which receive auto-generated ``Internal_N`` names),
and returns a ``Tree`` nested-dictionary structure: leaf names map to ``None``
and internal node names map to ``{child_name: branch_length}`` dicts.
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

Tree = Dict[str, Any]


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
