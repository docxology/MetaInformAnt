"""Batch spatially variable gene (SVG) detection for spatial transcriptomics.

Ranks genes in a spatial expression matrix by spatial autocorrelation, using
Moran's I on a shared spatial weights matrix. Descriptive statistics only:
the p-value and adjusted q-value come from the classical normal approximation
under the null of no spatial autocorrelation and carry no inferential
cross-species claim.

This complements per-vector :func:`morans_i` with a matrix-level interface
suited to whole-slide gene panels (Visium, MERFISH, Xenium aggregated counts).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

from metainformant.core.utils.logging import get_logger
from metainformant.spatial.analysis.autocorrelation import spatial_weights_matrix

logger = get_logger(__name__)

# Optional dependencies
try:
    import numpy as np
    from numpy.typing import NDArray
except ImportError:
    np = None  # type: ignore[assignment]
    NDArray = None  # type: ignore[assignment,misc]

try:
    from scipy.stats import false_discovery_control
except ImportError:  # pragma: no cover - scipy<1.11 fallback
    false_discovery_control = None  # type: ignore[assignment]


@dataclass
class SVGResult:
    """Batch spatially variable gene detection result.

    Attributes:
        genes: Gene names, aligned row-wise with the statistic arrays,
            sorted by Moran's I descending.
        morans_i: Moran's I per gene.
        z_scores: Standardized Z-score per gene.
        p_values: Two-sided normal-approximation p-value per gene.
        q_values: Benjamini-Hochberg adjusted p-values (None when unavailable).
        n_spots: Number of spatial spots used.
    """

    genes: list[str]
    morans_i: Any  # np.ndarray (n_genes,)
    z_scores: Any  # np.ndarray (n_genes,)
    p_values: Any  # np.ndarray (n_genes,)
    q_values: Optional[Any]
    n_spots: int

    def top_genes(self, n: int = 20) -> list[str]:
        """Return the ``n`` genes with the highest Moran's I, ranked descending."""
        order = np.argsort(-self.morans_i)
        return [self.genes[i] for i in order[: min(n, len(self.genes))]]


def detect_spatially_variable_genes(
    expression: Any,
    coordinates: Any,
    genes: Optional[list[str]] = None,
    *,
    method: str = "knn",
    k: int = 6,
    min_mean_expression: Optional[float] = None,
    min_variance: float = 0.0,
) -> SVGResult:
    """Rank all genes in a spatial expression matrix by Moran's I.

    Args:
        expression: Gene-by-spot or spot-by-gene matrix (auto-oriented to
            genes x spots). Genes failing the mean/variance filters are skipped.
        coordinates: Spot coordinates (n_spots x n_dims).
        genes: Optional gene names. Defaults to ``Gene0..GeneN-1``.
        method: Weight scheme passed to :func:`spatial_weights_matrix`
            ("knn", "distance", or "binary").
        k: Neighbors for the knn weight scheme.
        min_mean_expression: Drop genes with mean expression below this
            (None = no mean filter).
        min_variance: Drop genes with variance below this.

    Returns:
        SVGResult with one entry per retained gene, sorted by Moran's I descending.

    Raises:
        ImportError: If NumPy is unavailable.
        ValueError: If expression orientation/shape disagrees with the
            coordinate count, or the gene-name count is wrong.
    """
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")

    X = np.asarray(expression, dtype=np.float64)
    coords = np.asarray(coordinates, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError(f"expression must be 2-D, got shape {X.shape}")

    n_spots = coords.shape[0]
    if X.shape[0] == n_spots and X.shape[1] != n_spots:
        X = X.T  # orient to genes x spots
    if X.shape[1] != n_spots:
        raise ValueError(f"expression has {X.shape[1]} columns after orientation but {n_spots} coordinates")

    if genes is None:
        genes = [f"Gene{i}" for i in range(X.shape[0])]
    if len(genes) != X.shape[0]:
        raise ValueError(f"got {len(genes)} gene names for {X.shape[0]} rows")

    W = spatial_weights_matrix(coords, method=method, k=k)

    from metainformant.spatial.analysis.autocorrelation import morans_i

    gene_names: list[str] = []
    i_vals: list[float] = []
    z_vals: list[float] = []
    p_vals: list[float] = []

    for idx in range(X.shape[0]):
        row = X[idx]
        if min_mean_expression is not None and row.mean() < min_mean_expression:
            continue
        if row.var() <= min_variance:
            continue
        result = morans_i(row, W)
        gene_names.append(genes[idx])
        i_vals.append(result.I)
        z_vals.append(result.z_score)
        p_vals.append(result.p_value)

    if not gene_names:
        logger.warning("No genes passed the expression/variance filters")
        return SVGResult(
            genes=[],
            morans_i=np.array([], dtype=np.float64),
            z_scores=np.array([], dtype=np.float64),
            p_values=np.array([], dtype=np.float64),
            q_values=None,
            n_spots=n_spots,
        )

    i_arr = np.asarray(i_vals, dtype=np.float64)
    order = np.argsort(-i_arr)

    q_arr: Optional[Any] = None
    if false_discovery_control is not None:
        try:
            q_arr = false_discovery_control(np.asarray(p_vals, dtype=np.float64), method="bh")
        except Exception:  # pragma: no cover - degenerate p arrays
            q_arr = None
    if q_arr is not None:
        q_arr = q_arr[order]

    logger.info("SVG detection: %d genes ranked across %d spots", len(gene_names), n_spots)
    return SVGResult(
        genes=[gene_names[i] for i in order],
        morans_i=i_arr[order],
        z_scores=np.asarray(z_vals, dtype=np.float64)[order],
        p_values=np.asarray(p_vals, dtype=np.float64)[order],
        q_values=q_arr,
        n_spots=n_spots,
    )
