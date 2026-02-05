"""Spatial transcriptomics visualization plots.

Publication-quality plotting functions for spatial data including expression
scatter plots, tissue image overlays, cell type maps, neighborhood graphs,
deconvolution pie charts, and spatial autocorrelation (LISA) maps.

All plot functions save to disk and return the matplotlib Figure object.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies
try:
    import numpy as np
except ImportError:
    np = None  # type: ignore[assignment]

try:
    import matplotlib

    matplotlib.use("Agg")  # Non-interactive backend for file output
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Wedge
    import matplotlib.colors as mcolors
except ImportError:
    plt = None  # type: ignore[assignment]
    PatchCollection = None  # type: ignore[assignment,misc]
    Wedge = None  # type: ignore[assignment,misc]
    mcolors = None  # type: ignore[assignment]

try:
    import seaborn as sns
except ImportError:
    sns = None  # type: ignore[assignment]

try:
    from scipy import sparse as sp_sparse
except ImportError:
    sp_sparse = None  # type: ignore[assignment]


def _ensure_plotting_deps() -> None:
    """Check that matplotlib is available."""
    if plt is None:
        raise ImportError("matplotlib is required for plotting: uv pip install matplotlib")
    if np is None:
        raise ImportError("NumPy is required: uv pip install numpy")


def _save_figure(fig: Any, output_path: str | Path, dpi: int = 150) -> None:
    """Save a matplotlib figure to disk, creating directories as needed."""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out), dpi=dpi, bbox_inches="tight", facecolor="white")
    logger.info(f"Saved plot: {out}")


def plot_spatial_scatter(
    coordinates: Any,
    values: Any,
    output_path: str | Path,
    *,
    cmap: str = "viridis",
    title: str = "Spatial Scatter",
    point_size: float = 10.0,
    alpha: float = 0.8,
    figsize: tuple[float, float] = (8, 8),
    colorbar_label: str = "",
    vmin: float | None = None,
    vmax: float | None = None,
) -> Any:
    """Create a spatial scatter plot colored by continuous or categorical values.

    Args:
        coordinates: Spatial coordinates (n x 2).
        values: Values for coloring (length n). Can be numeric or categorical.
        output_path: Path to save the plot image.
        cmap: Matplotlib colormap name.
        title: Plot title.
        point_size: Scatter point size.
        alpha: Point transparency.
        figsize: Figure dimensions (width, height) in inches.
        colorbar_label: Label for the colorbar.
        vmin: Minimum value for color scale.
        vmax: Maximum value for color scale.

    Returns:
        Matplotlib Figure object.
    """
    _ensure_plotting_deps()

    coords = np.asarray(coordinates, dtype=np.float64)
    vals = np.asarray(values)

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Check if values are categorical (string/object) or numeric
    is_categorical = vals.dtype.kind in ("U", "S", "O")

    if is_categorical:
        unique_vals = sorted(set(vals.tolist()))
        color_map = {}
        cmap_obj = plt.cm.get_cmap("tab20", len(unique_vals))
        for i, v in enumerate(unique_vals):
            color_map[v] = cmap_obj(i)

        for v in unique_vals:
            mask = vals == v
            ax.scatter(
                coords[mask, 1], coords[mask, 0],
                c=[color_map[v]],
                s=point_size,
                alpha=alpha,
                label=str(v),
                edgecolors="none",
            )
        ax.legend(
            bbox_to_anchor=(1.05, 1), loc="upper left",
            frameon=False, fontsize=8, markerscale=2,
        )
    else:
        vals_float = vals.astype(np.float64)
        sc = ax.scatter(
            coords[:, 1], coords[:, 0],
            c=vals_float,
            cmap=cmap,
            s=point_size,
            alpha=alpha,
            vmin=vmin,
            vmax=vmax,
            edgecolors="none",
        )
        cbar = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
        if colorbar_label:
            cbar.set_label(colorbar_label, fontsize=10)

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("X coordinate", fontsize=10)
    ax.set_ylabel("Y coordinate", fontsize=10)
    ax.invert_yaxis()  # Match image convention (origin top-left)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=8)

    fig.tight_layout()
    _save_figure(fig, output_path)
    plt.close(fig)
    return fig


def plot_tissue_overlay(
    coordinates: Any,
    values: Any,
    tissue_image: Any,
    output_path: str | Path,
    *,
    cmap: str = "hot",
    title: str = "Tissue Overlay",
    point_size: float = 20.0,
    alpha: float = 0.6,
    figsize: tuple[float, float] = (10, 10),
    scale_factor: float = 1.0,
) -> Any:
    """Overlay expression values on a tissue H&E image.

    Args:
        coordinates: Spatial coordinates (n x 2) in pixel space.
        values: Expression values or cluster labels (length n).
        tissue_image: Tissue image array (H, W, C).
        output_path: Path to save the plot.
        cmap: Colormap for expression overlay.
        title: Plot title.
        point_size: Scatter point size.
        alpha: Overlay transparency.
        figsize: Figure dimensions.
        scale_factor: Scale factor to convert coordinates to image pixel space.

    Returns:
        Matplotlib Figure object.
    """
    _ensure_plotting_deps()

    coords = np.asarray(coordinates, dtype=np.float64) * scale_factor
    vals = np.asarray(values)
    img = np.asarray(tissue_image)

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Display tissue image
    ax.imshow(img, aspect="equal")

    # Overlay spots
    is_categorical = vals.dtype.kind in ("U", "S", "O")

    if is_categorical:
        unique_vals = sorted(set(vals.tolist()))
        cmap_obj = plt.cm.get_cmap("tab20", len(unique_vals))
        color_map = {v: cmap_obj(i) for i, v in enumerate(unique_vals)}
        for v in unique_vals:
            mask = vals == v
            ax.scatter(
                coords[mask, 1], coords[mask, 0],
                c=[color_map[v]], s=point_size, alpha=alpha,
                label=str(v), edgecolors="none",
            )
        ax.legend(
            bbox_to_anchor=(1.05, 1), loc="upper left",
            frameon=False, fontsize=8, markerscale=2,
        )
    else:
        vals_float = vals.astype(np.float64)
        sc = ax.scatter(
            coords[:, 1], coords[:, 0],
            c=vals_float, cmap=cmap, s=point_size, alpha=alpha,
            edgecolors="none",
        )
        fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.02)

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.axis("off")

    fig.tight_layout()
    _save_figure(fig, output_path)
    plt.close(fig)
    return fig


def plot_gene_expression_map(
    spatial_data: Any,
    gene: str,
    output_path: str | Path,
    *,
    cmap: str = "Reds",
    figsize: tuple[float, float] = (8, 8),
) -> Any:
    """Plot spatial expression map for a single gene.

    Args:
        spatial_data: A SpatialDataset (or compatible object with .expression,
            .coordinates, .gene_names attributes).
        gene: Gene name to plot.
        output_path: Path to save the plot.
        cmap: Colormap for expression.
        figsize: Figure dimensions.

    Returns:
        Matplotlib Figure object.

    Raises:
        ValueError: If gene is not found in the dataset.
    """
    _ensure_plotting_deps()

    gene_names = spatial_data.gene_names
    if gene not in gene_names:
        raise ValueError(f"Gene '{gene}' not found. Available: {gene_names[:10]}...")

    gene_idx = gene_names.index(gene)

    expression = spatial_data.expression
    if sp_sparse is not None and sp_sparse.issparse(expression):
        expr_values = np.asarray(expression[:, gene_idx].toarray()).flatten()
    elif hasattr(expression, "toarray"):
        expr_values = np.asarray(expression[:, gene_idx].toarray()).flatten()
    else:
        expr_values = np.asarray(expression[:, gene_idx]).flatten()

    return plot_spatial_scatter(
        spatial_data.coordinates,
        expr_values,
        output_path,
        cmap=cmap,
        title=f"Expression: {gene}",
        colorbar_label="Expression",
        figsize=figsize,
    )


def plot_cell_type_map(
    spatial_data: Any,
    cell_types: Any,
    output_path: str | Path,
    *,
    figsize: tuple[float, float] = (10, 10),
    point_size: float = 15.0,
) -> Any:
    """Plot spatial distribution of cell types.

    Args:
        spatial_data: A SpatialDataset or compatible object with .coordinates attribute.
        cell_types: Cell type labels (length n_spots), string or integer.
        output_path: Path to save the plot.
        figsize: Figure dimensions.
        point_size: Scatter point size.

    Returns:
        Matplotlib Figure object.
    """
    _ensure_plotting_deps()

    coords = np.asarray(spatial_data.coordinates if hasattr(spatial_data, "coordinates") else spatial_data)
    types = np.asarray(cell_types)

    return plot_spatial_scatter(
        coords,
        types,
        output_path,
        title="Cell Type Distribution",
        point_size=point_size,
        figsize=figsize,
    )


def plot_neighborhood_graph(
    coordinates: Any,
    spatial_graph: Any,
    output_path: str | Path,
    *,
    node_colors: Any | None = None,
    title: str = "Spatial Neighborhood Graph",
    figsize: tuple[float, float] = (10, 10),
    node_size: float = 5.0,
    edge_alpha: float = 0.2,
    edge_width: float = 0.3,
) -> Any:
    """Plot the spatial neighborhood graph on tissue coordinates.

    Draws edges between connected spots and colors nodes by cluster or expression.

    Args:
        coordinates: Spatial coordinates (n x 2).
        spatial_graph: Sparse adjacency matrix (n x n).
        output_path: Path to save the plot.
        node_colors: Values to color nodes by (length n). If None, uniform color.
        title: Plot title.
        figsize: Figure dimensions.
        node_size: Node marker size.
        edge_alpha: Edge line transparency.
        edge_width: Edge line width.

    Returns:
        Matplotlib Figure object.
    """
    _ensure_plotting_deps()

    coords = np.asarray(coordinates, dtype=np.float64)
    n = coords.shape[0]

    if sp_sparse is not None and sp_sparse.issparse(spatial_graph):
        adj = spatial_graph
    else:
        adj = sp_sparse.csr_matrix(np.asarray(spatial_graph)) if sp_sparse else spatial_graph

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Draw edges
    if sp_sparse is not None:
        sources, targets = adj.nonzero()
        for s, t in zip(sources, targets):
            if s < t:  # draw each edge once
                ax.plot(
                    [coords[s, 1], coords[t, 1]],
                    [coords[s, 0], coords[t, 0]],
                    color="gray",
                    alpha=edge_alpha,
                    linewidth=edge_width,
                    zorder=1,
                )

    # Draw nodes
    if node_colors is not None:
        nc = np.asarray(node_colors)
        is_cat = nc.dtype.kind in ("U", "S", "O")
        if is_cat:
            unique_vals = sorted(set(nc.tolist()))
            cmap_obj = plt.cm.get_cmap("tab20", len(unique_vals))
            color_map = {v: cmap_obj(i) for i, v in enumerate(unique_vals)}
            for v in unique_vals:
                mask = nc == v
                ax.scatter(
                    coords[mask, 1], coords[mask, 0],
                    c=[color_map[v]], s=node_size,
                    label=str(v), edgecolors="none", zorder=2,
                )
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False, fontsize=8)
        else:
            nc_float = nc.astype(np.float64)
            sc = ax.scatter(
                coords[:, 1], coords[:, 0],
                c=nc_float, cmap="viridis", s=node_size,
                edgecolors="none", zorder=2,
            )
            fig.colorbar(sc, ax=ax, shrink=0.5)
    else:
        ax.scatter(
            coords[:, 1], coords[:, 0],
            c="steelblue", s=node_size,
            edgecolors="none", zorder=2,
        )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xlabel("X coordinate", fontsize=10)
    ax.set_ylabel("Y coordinate", fontsize=10)
    ax.tick_params(labelsize=8)

    fig.tight_layout()
    _save_figure(fig, output_path)
    plt.close(fig)
    return fig


def plot_interaction_heatmap(
    interaction_matrix: Any,
    output_path: str | Path,
    *,
    cell_type_names: list[str] | None = None,
    title: str = "Cell Type Interaction",
    cmap: str = "RdBu_r",
    figsize: tuple[float, float] = (8, 7),
    annot: bool = True,
) -> Any:
    """Plot a cell type interaction heatmap.

    Args:
        interaction_matrix: Interaction scores (n_types x n_types).
        output_path: Path to save the plot.
        cell_type_names: Labels for rows/columns.
        title: Plot title.
        cmap: Colormap.
        figsize: Figure dimensions.
        annot: If True, annotate cells with values.

    Returns:
        Matplotlib Figure object.
    """
    _ensure_plotting_deps()

    mat = np.asarray(interaction_matrix, dtype=np.float64)
    n = mat.shape[0]

    if cell_type_names is None:
        cell_type_names = [str(i) for i in range(n)]

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    if sns is not None:
        sns.heatmap(
            mat,
            xticklabels=cell_type_names,
            yticklabels=cell_type_names,
            cmap=cmap,
            center=0,
            annot=annot,
            fmt=".2f" if annot else "",
            square=True,
            ax=ax,
            cbar_kws={"shrink": 0.8},
        )
    else:
        im = ax.imshow(mat, cmap=cmap, aspect="equal")
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(cell_type_names, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(cell_type_names, fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.8)

        if annot:
            for i in range(n):
                for j in range(n):
                    ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=7)

    ax.set_title(title, fontsize=14, fontweight="bold")

    fig.tight_layout()
    _save_figure(fig, output_path)
    plt.close(fig)
    return fig


def plot_deconvolution_pie(
    coordinates: Any,
    fractions: Any,
    output_path: str | Path,
    *,
    cell_type_names: list[str] | None = None,
    title: str = "Cell Type Composition",
    figsize: tuple[float, float] = (10, 10),
    pie_radius: float | None = None,
    min_fraction: float = 0.05,
) -> Any:
    """Plot pie charts per spatial spot showing cell type fractions.

    Each spot is represented as a small pie chart showing the estimated
    cell type proportions from deconvolution.

    Args:
        coordinates: Spatial coordinates (n_spots x 2).
        fractions: Cell type fractions (n_spots x n_types), rows sum to 1.
        output_path: Path to save the plot.
        cell_type_names: Names of cell types.
        title: Plot title.
        figsize: Figure dimensions.
        pie_radius: Radius of each pie chart. If None, auto-calculated.
        min_fraction: Minimum fraction to display (smaller merged into "other").

    Returns:
        Matplotlib Figure object.
    """
    _ensure_plotting_deps()

    if Wedge is None:
        raise ImportError("matplotlib is required: uv pip install matplotlib")

    coords = np.asarray(coordinates, dtype=np.float64)
    frac = np.asarray(fractions, dtype=np.float64)
    n_spots = coords.shape[0]
    n_types = frac.shape[1]

    if cell_type_names is None:
        cell_type_names = [f"Type {i}" for i in range(n_types)]

    # Auto-calculate pie radius based on average spacing
    if pie_radius is None:
        if n_spots > 1:
            from scipy.spatial import KDTree

            tree = KDTree(coords)
            dists, _ = tree.query(coords, k=2)
            pie_radius = float(np.median(dists[:, 1])) * 0.4
        else:
            pie_radius = 10.0

    # Color palette
    cmap_obj = plt.cm.get_cmap("tab20", n_types)
    colors = [cmap_obj(i) for i in range(n_types)]

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    for si in range(n_spots):
        cx, cy = coords[si, 1], coords[si, 0]  # x=col, y=row
        spot_frac = frac[si, :]

        # Filter small fractions
        significant = spot_frac >= min_fraction
        if not significant.any():
            continue

        # Draw pie wedges
        start_angle = 0.0
        for ti in range(n_types):
            if spot_frac[ti] < min_fraction:
                continue
            sweep = spot_frac[ti] * 360.0
            wedge = Wedge(
                (cx, cy), pie_radius,
                start_angle, start_angle + sweep,
                facecolor=colors[ti],
                edgecolor="white",
                linewidth=0.3,
            )
            ax.add_patch(wedge)
            start_angle += sweep

    # Legend
    legend_handles = []
    for ti in range(n_types):
        patch = plt.Rectangle((0, 0), 1, 1, facecolor=colors[ti])
        legend_handles.append(patch)
    ax.legend(
        legend_handles, cell_type_names,
        bbox_to_anchor=(1.05, 1), loc="upper left",
        frameon=False, fontsize=8,
    )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlim(coords[:, 1].min() - pie_radius * 3, coords[:, 1].max() + pie_radius * 3)
    ax.set_ylim(coords[:, 0].min() - pie_radius * 3, coords[:, 0].max() + pie_radius * 3)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xlabel("X coordinate", fontsize=10)
    ax.set_ylabel("Y coordinate", fontsize=10)
    ax.tick_params(labelsize=8)

    fig.tight_layout()
    _save_figure(fig, output_path)
    plt.close(fig)
    return fig


def plot_spatial_autocorrelation(
    coordinates: Any,
    local_scores: Any,
    output_path: str | Path,
    *,
    cluster_labels: list[str] | None = None,
    title: str = "Local Spatial Autocorrelation (LISA)",
    figsize: tuple[float, float] = (10, 10),
    point_size: float = 15.0,
) -> Any:
    """Plot LISA cluster map showing spatial autocorrelation patterns.

    Colors spots by their LISA cluster classification:
    HH (red), LL (blue), HL (pink), LH (lightblue), NS (gray).

    Args:
        coordinates: Spatial coordinates (n x 2).
        local_scores: Local Moran's I values (length n) or similar local statistic.
        output_path: Path to save the plot.
        cluster_labels: LISA cluster labels ("HH", "LL", "HL", "LH", "NS").
            If None, uses local_scores as continuous values.
        title: Plot title.
        figsize: Figure dimensions.
        point_size: Scatter point size.

    Returns:
        Matplotlib Figure object.
    """
    _ensure_plotting_deps()

    coords = np.asarray(coordinates, dtype=np.float64)

    if cluster_labels is not None:
        # Categorical LISA map
        labels = np.asarray(cluster_labels)

        lisa_colors = {
            "HH": "#d73027",    # Hot-spot: red
            "LL": "#4575b4",    # Cold-spot: blue
            "HL": "#f4a582",    # High-Low outlier: pink
            "LH": "#abd9e9",    # Low-High outlier: light blue
            "NS": "#e0e0e0",    # Not significant: gray
        }

        fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Plot NS first (background), then significant on top
        for label_type in ["NS", "LH", "HL", "LL", "HH"]:
            mask = labels == label_type
            if not mask.any():
                continue
            color = lisa_colors.get(label_type, "#e0e0e0")
            ax.scatter(
                coords[mask, 1], coords[mask, 0],
                c=color, s=point_size,
                label=f"{label_type} (n={mask.sum()})",
                edgecolors="none", alpha=0.8,
            )

        ax.legend(
            bbox_to_anchor=(1.05, 1), loc="upper left",
            frameon=False, fontsize=9, markerscale=2,
        )

        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set_xlabel("X coordinate", fontsize=10)
        ax.set_ylabel("Y coordinate", fontsize=10)
        ax.tick_params(labelsize=8)

        fig.tight_layout()
        _save_figure(fig, output_path)
        plt.close(fig)
        return fig

    else:
        # Continuous local statistics map
        scores = np.asarray(local_scores, dtype=np.float64)
        return plot_spatial_scatter(
            coords,
            scores,
            output_path,
            cmap="RdBu_r",
            title=title,
            point_size=point_size,
            colorbar_label="Local Moran's I",
            figsize=figsize,
        )
