"""Match gene names to expression columns (shared by the communication and deconvolution modules).

Ligand, receptor and marker genes are matched to the columns of an expression matrix by their exact names,
and by nothing else: a name that is not among the matrix's gene names is absent, never mapped to some
other column. Pure Python on purpose, so that it imports without numpy.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any


def expression_and_gene_names(expression: Any, gene_names: Sequence[str] | None) -> tuple[Any, list[str]]:
    """Split a matrix-or-dataset argument into the matrix and the names of its columns.

    ``expression`` may be a plain matrix (numpy array, list of lists, scipy sparse; ``gene_names`` is then
    required), a pandas DataFrame (cells x genes; its columns are the names), or a dataset object with
    ``.expression`` and ``.gene_names`` (``SpatialDataset``, ``XeniumDataset``, ``MERFISHDataset``).
    Explicit ``gene_names`` always take precedence over the names the object carries.

    Raises:
        ValueError: If no gene names are available.
    """
    matrix = expression
    carried: Sequence[str] | None = None
    if hasattr(expression, "expression") and hasattr(expression, "gene_names"):
        matrix = expression.expression
        carried = expression.gene_names
    elif hasattr(expression, "columns") and hasattr(expression, "to_numpy"):
        matrix = expression.to_numpy()
        carried = list(expression.columns)

    names = gene_names if gene_names is not None else carried
    if names is None:
        raise ValueError(
            "gene_names is required: genes are matched to expression columns by name, so pass gene_names=... "
            "(one name per column) or a dataset / DataFrame that carries them"
        )
    return matrix, [str(name) for name in names]


def matrix_shape(matrix: Any) -> tuple[int, int]:
    """``(n_rows, n_columns)`` of a matrix-like object, which must be two-dimensional.

    Raises:
        ValueError: If the object is not two-dimensional.
    """
    shape = getattr(matrix, "shape", None)
    if shape is None:
        import numpy as np  # callers have already checked that numpy is available

        shape = np.asarray(matrix).shape
    if len(shape) != 2:
        raise ValueError(f"expression must be a 2-D matrix (spots or cells x genes), got shape {tuple(shape)}")
    return int(shape[0]), int(shape[1])


def expression_columns(matrix: Any, columns: Sequence[int]) -> Any:
    """The given columns of a dense or scipy-sparse matrix as a dense float64 array.

    Sparse input is densified one block of columns at a time, never as a whole.
    """
    import numpy as np  # callers have already checked that numpy is available

    selected = list(columns)
    if hasattr(matrix, "tocsc"):
        return np.asarray(matrix.tocsc()[:, selected].toarray(), dtype=np.float64)
    return np.asarray(matrix, dtype=np.float64)[:, selected]


def columns_of(wanted: Iterable[str], gene_names: Sequence[str], *, what: str = "gene") -> dict[str, int]:
    """Column index of each wanted name that is present; absent names are simply left out.

    Raises:
        ValueError: If a wanted name labels more than one column, because then no column is the right one.
    """
    wanted_set = set(wanted)
    positions: dict[str, list[int]] = {}
    for index, name in enumerate(gene_names):
        if name in wanted_set:
            positions.setdefault(name, []).append(index)

    ambiguous = {name: columns for name, columns in positions.items() if len(columns) > 1}
    if ambiguous:
        example = next(iter(ambiguous.items()))
        raise ValueError(
            f"{len(ambiguous)} {what} name(s) label more than one expression column, so they cannot be matched "
            f"unambiguously (for example {example[0]!r} labels columns {example[1][:4]}); make gene_names "
            "unique first, for example by collapsing the duplicate columns"
        )
    return {name: columns[0] for name, columns in positions.items()}
