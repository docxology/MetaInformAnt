from __future__ import annotations

from typing import Any

import matplotlib

# Non-interactive for tests/headless
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402

try:  # optional Bio.Phylo ASCII exists already; we also expose matplotlib plot
    from Bio import Phylo
except Exception:  # pragma: no cover - optional runtime
    Phylo = None  # type: ignore


def plot_phylo_tree(tree: Any, *, ax: plt.Axes | None = None) -> plt.Axes:
    """Plot a Biopython Phylo tree to matplotlib Axes.

    Accepts any object compatible with Bio.Phylo.draw().
    """
    if Phylo is None:
        raise RuntimeError("Biopython Phylo not available")
    if ax is None:
        _, ax = plt.subplots()
    Phylo.draw(tree, do_show=False, axes=ax)
    return ax
