"""PCA-related population structure visualization for GWAS.

This module provides PCA plots, scree plots, multi-panel PCA, and 3D PCA
for visualizing population stratification and structure.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def _load_pca_from_file(pca_file: str | Path) -> Dict[str, Any]:
    """Load PCA data from a TSV file.

    Expects a TSV with columns: sample_id, PC1, PC2, ... PCn.

    Args:
        pca_file: Path to TSV file containing PCA data.

    Returns:
        Dictionary with 'pcs', 'explained_variance', 'sample_names'.
    """
    pca_file = Path(pca_file)
    lines = pca_file.read_text().strip().split("\n")
    header = lines[0].split("\t")
    pc_cols = [c for c in header if c.upper().startswith("PC")]
    pc_indices = [header.index(c) for c in pc_cols]

    sample_names = []
    pcs_data = []
    for line in lines[1:]:
        parts = line.split("\t")
        sample_names.append(parts[0])
        pcs_data.append([float(parts[i]) for i in pc_indices])

    pcs_array = np.array(pcs_data)
    # Compute approximate explained variance from variance of each PC column
    variances = np.var(pcs_array, axis=0)
    total_var = np.sum(variances)
    explained_variance = variances / total_var if total_var > 0 else np.zeros(len(pc_cols))

    return {
        "pcs": pcs_array,
        "explained_variance": explained_variance,
        "sample_names": sample_names,
    }


def pca_plot(
    pca_data: Dict[str, Any] | str | Path,
    output_file: Optional[str | Path] = None,
    title: str = "PCA Plot",
    pc1: int = 0,
    pc2: int = 1,
) -> Dict[str, Any]:
    """Create a PCA plot for population structure visualization.

    Args:
        pca_data: Dictionary containing PCA data with keys:
                 'pcs': 2D array of principal components
                 'explained_variance': array of explained variance ratios
                 'sample_names': list of sample names (optional)
                 'populations': list of population labels (optional)
                 OR a file path (str/Path) to a TSV file with PCA data.
        output_file: Optional output file path
        title: Plot title
        pc1: Index of first principal component to plot (0-based)
        pc2: Index of second principal component to plot (0-based)

    Returns:
        Dictionary with status and metadata.

    Example:
        >>> data = {'pcs': pcs_array, 'explained_variance': var_ratios}
        >>> result = pca_plot(data)
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib or seaborn not available for PCA plot")
        return {"status": "failed", "error": "matplotlib/seaborn not available"}

    try:
        # If pca_data is a file path, load it
        if isinstance(pca_data, (str, Path)):
            pca_data = _load_pca_from_file(pca_data)

        # Validate input data
        required_keys = ["pcs", "explained_variance"]
        missing_keys = [key for key in required_keys if key not in pca_data]
        if missing_keys:
            logger.error(f"Missing required keys in pca_data: {missing_keys}")
            return {"status": "failed", "error": f"Missing required keys: {missing_keys}"}

        pcs = pca_data["pcs"]
        explained_var = pca_data["explained_variance"]

        if len(pcs.shape) != 2:
            logger.error("PCA data 'pcs' must be a 2D array")
            return {"status": "failed", "error": "PCA data 'pcs' must be a 2D array"}

        if pc1 >= pcs.shape[1] or pc2 >= pcs.shape[1]:
            logger.error(f"PC indices {pc1}, {pc2} out of range for {pcs.shape[1]} components")
            return {"status": "failed", "error": f"PC indices out of range"}

        fig, ax = plt.subplots(1, 1, figsize=(10, 8))

        # Extract PC data
        x_data = pcs[:, pc1]
        y_data = pcs[:, pc2]

        # Check for population labels
        populations = pca_data.get("populations")
        sample_names = pca_data.get("sample_names")

        if populations is not None and len(populations) == len(x_data):
            # Color by population
            unique_pops = list(set(populations))
            colors = plt.cm.tab10(np.linspace(0, 1, len(unique_pops)))

            for i, pop in enumerate(unique_pops):
                mask = [p == pop for p in populations]
                ax.scatter(x_data[mask], y_data[mask], c=[colors[i]], label=pop, alpha=0.7, s=50)
            ax.legend(title="Population", bbox_to_anchor=(1.05, 1), loc="upper left")
        else:
            # Simple scatter plot
            ax.scatter(x_data, y_data, alpha=0.7, s=50, c="blue")

        # Add sample labels if provided and not too many samples
        if sample_names is not None and len(sample_names) <= 50:
            for i, name in enumerate(sample_names):
                ax.annotate(
                    name, (x_data[i], y_data[i]), xytext=(5, 5), textcoords="offset points", fontsize=8, alpha=0.8
                )

        # Labels and title
        var1 = explained_var[pc1] * 100
        var2 = explained_var[pc2] * 100

        ax.set_xlabel(f"PC{pc1+1} ({var1:.1f}% variance)", fontsize=12)
        ax.set_ylabel(f"PC{pc2+1} ({var2:.1f}% variance)", fontsize=12)
        ax.set_title(title, fontsize=14, pad=20)
        ax.grid(True, alpha=0.3)

        # Add center lines
        ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
        ax.axvline(x=0, color="gray", linestyle="--", alpha=0.5)

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            logger.info(f"Saved PCA plot to {output_file}")

        plt.close(fig)
        return {
            "status": "success",
            "n_samples": pcs.shape[0],
            "n_components": pcs.shape[1],
            "output_path": str(output_file) if output_file else None,
        }

    except Exception as e:
        logger.error(f"Error creating PCA plot: {e}")
        return {"status": "failed", "error": str(e)}


def pca_scree_plot(
    explained_variance: List[float] | str | Path,
    output_file: Optional[str | Path] = None,
    title: str = "PCA Scree Plot",
) -> Dict[str, Any]:
    """Create a scree plot showing explained variance by principal component.

    Args:
        explained_variance: List of explained variance ratios (0-1),
            OR a file path to a text file with one eigenvalue per line.
        output_file: Optional output file path
        title: Plot title

    Returns:
        Dictionary with status and metadata.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for scree plot")
        return {"status": "failed", "error": "matplotlib not available"}

    try:
        # If explained_variance is a file path, load eigenvalues from it
        if isinstance(explained_variance, (str, Path)):
            eigenval_path = Path(explained_variance)
            lines = eigenval_path.read_text().strip().split("\n")
            explained_variance = [float(line.strip()) for line in lines if line.strip()]

        if not explained_variance:
            logger.error("No explained variance data provided")
            return {"status": "failed", "error": "No explained variance data provided"}

        # Convert to percentages if not already
        if max(explained_variance) <= 1.0:
            explained_variance_pct = [v * 100 for v in explained_variance]
        else:
            explained_variance_pct = list(explained_variance)

        fig, ax = plt.subplots(figsize=(10, 6))

        # Create scree plot
        components = list(range(1, len(explained_variance) + 1))
        bars = ax.bar(components, explained_variance_pct, color="skyblue", edgecolor="navy", alpha=0.7)

        # Add cumulative variance line
        cumulative = np.cumsum(explained_variance_pct)
        ax.plot(components, cumulative, "r-o", linewidth=2, markersize=4, label="Cumulative")

        # Add value labels on bars
        for bar, pct in zip(bars, explained_variance_pct):
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2.0, height + 0.5, f"{pct:.1f}%", ha="center", va="bottom", fontsize=8
            )

        ax.set_xlabel("Principal Component", fontsize=12)
        ax.set_ylabel("Explained Variance (%)", fontsize=12)
        ax.set_title(title, fontsize=14, pad=20)
        ax.set_xticks(components)
        ax.grid(True, alpha=0.3, axis="y")
        ax.legend()

        # Add elbow detection heuristic (find where additional variance drops below 5%)
        if len(explained_variance_pct) > 1:
            for i in range(1, len(explained_variance_pct)):
                if explained_variance_pct[i] < 5.0:
                    ax.axvline(x=i + 0.5, color="red", linestyle="--", alpha=0.7, label=f"Elbow at PC{i+1}")
                    break

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            logger.info(f"Saved PCA scree plot to {output_file}")

        plt.close(fig)
        return {
            "status": "success",
            "n_components": len(explained_variance),
            "output_path": str(output_file) if output_file else None,
        }

    except Exception as e:
        logger.error(f"Error creating scree plot: {e}")
        return {"status": "failed", "error": str(e)}


def pca_multi_panel(
    pca_data: Dict[str, Any],
    metadata: Optional[Dict[str, Dict]] = None,
    output_file: Optional[str | Path] = None,
    pairs: Optional[List[tuple[int, int]]] = None,
    color_by: str = "population",
    title: str = "PCA Multi-Panel",
) -> Dict[str, Any]:
    """Create a grid of PCA scatter plots for different PC pairs.

    Generates a multi-panel figure with one scatter plot per PC pair.  When
    ``metadata`` is provided each sample is coloured by the ``color_by`` field.

    Args:
        pca_data: Dictionary with keys:
            - ``"pcs"``: list of lists (samples x components) or 2-D ndarray.
            - ``"explained_variance_ratio"``: list of floats (one per component).
            - ``"sample_ids"`` (optional): list of str identifying samples.
        metadata: Optional mapping from sample id to a dict of annotations.
            The ``color_by`` key inside each annotation dict determines colour.
        output_file: Optional file path for saving the figure.
        pairs: List of (pc_i, pc_j) tuples to plot.  Defaults to
            ``[(0, 1), (0, 2), (1, 2)]``, pruned to available components.
        color_by: Metadata field name used for colouring points.
        title: Super-title for the figure.

    Returns:
        Dictionary with ``"status"``, ``"output_path"``, ``"n_samples"``, and
        ``"n_panels"`` keys.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for PCA multi-panel plot")
        return {"status": "skipped", "output_path": None, "n_samples": 0, "n_panels": 0}

    try:
        pcs = np.asarray(pca_data.get("pcs", []), dtype=np.float64)
        explained = pca_data.get("explained_variance_ratio", [])
        sample_ids = pca_data.get("sample_ids", [])

        if pcs.ndim != 2 or pcs.shape[0] == 0:
            logger.error("pca_data['pcs'] must be a non-empty 2D structure")
            return {"status": "failed", "output_path": None, "n_samples": 0, "n_panels": 0}

        n_samples, n_pcs = pcs.shape

        # Build default pairs and prune to available components
        if pairs is None:
            pairs = [(0, 1), (0, 2), (1, 2)]
        valid_pairs = [(a, b) for a, b in pairs if a < n_pcs and b < n_pcs]
        if not valid_pairs:
            logger.error("No valid PC pairs for the number of available components")
            return {"status": "failed", "output_path": None, "n_samples": n_samples, "n_panels": 0}

        n_panels = len(valid_pairs)
        n_cols = min(n_panels, 3)
        n_rows = (n_panels + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows), squeeze=False)

        # Resolve per-sample colour labels from metadata
        color_labels: Optional[List[str]] = None
        if metadata is not None and sample_ids:
            color_labels = []
            for sid in sample_ids:
                entry = metadata.get(sid, {})
                color_labels.append(str(entry.get(color_by, "unknown")))

        unique_labels: Optional[List[str]] = None
        label_to_color: Optional[Dict[str, Any]] = None
        if color_labels is not None:
            unique_labels = sorted(set(color_labels))
            cmap = plt.cm.tab10(np.linspace(0, 1, max(len(unique_labels), 1)))
            label_to_color = {lab: cmap[i] for i, lab in enumerate(unique_labels)}

        for idx, (pc_a, pc_b) in enumerate(valid_pairs):
            row, col = divmod(idx, n_cols)
            ax = axes[row][col]

            if color_labels is not None and label_to_color is not None and unique_labels is not None:
                for lab in unique_labels:
                    mask = np.array([cl == lab for cl in color_labels])
                    ax.scatter(pcs[mask, pc_a], pcs[mask, pc_b], label=lab, alpha=0.7, s=40)
            else:
                ax.scatter(pcs[:, pc_a], pcs[:, pc_b], alpha=0.7, s=40, c="steelblue")

            xlabel = f"PC{pc_a + 1}"
            ylabel = f"PC{pc_b + 1}"
            if pc_a < len(explained):
                xlabel += f" ({explained[pc_a] * 100:.1f}%)"
            if pc_b < len(explained):
                ylabel += f" ({explained[pc_b] * 100:.1f}%)"
            ax.set_xlabel(xlabel, fontsize=10)
            ax.set_ylabel(ylabel, fontsize=10)
            ax.grid(True, alpha=0.3)

        # Hide unused axes
        for idx in range(n_panels, n_rows * n_cols):
            row, col = divmod(idx, n_cols)
            axes[row][col].set_visible(False)

        # Shared legend from the first axes
        if unique_labels is not None and len(unique_labels) > 0:
            handles, labels = axes[0][0].get_legend_handles_labels()
            if handles:
                fig.legend(handles, labels, loc="upper right", title=color_by, fontsize=9)

        fig.suptitle(title, fontsize=14, y=1.02)
        fig.tight_layout()

        output_path_str: Optional[str] = None
        if output_file:
            out = Path(output_file)
            out.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out, dpi=300, bbox_inches="tight")
            output_path_str = str(out)
            logger.info(f"Saved PCA multi-panel plot to {out}")

        plt.close(fig)
        return {
            "status": "success",
            "output_path": output_path_str,
            "n_samples": n_samples,
            "n_panels": n_panels,
        }

    except Exception as e:
        logger.error(f"Error creating PCA multi-panel plot: {e}")
        return {"status": "failed", "output_path": None, "n_samples": 0, "n_panels": 0}


def pca_3d(
    pca_data: Dict[str, Any],
    metadata: Optional[Dict[str, Dict]] = None,
    output_file: Optional[str | Path] = None,
    components: tuple[int, int, int] = (0, 1, 2),
    color_by: str = "population",
    title: str = "3D PCA",
) -> Dict[str, Any]:
    """Create a static 3D PCA scatter plot.

    Uses matplotlib's ``Axes3D`` projection to render three principal components
    in a single 3-D scatter.  Points are coloured by the ``color_by`` metadata
    field when ``metadata`` is supplied.

    Args:
        pca_data: Dictionary with keys:
            - ``"pcs"``: list of lists (samples x components) or 2-D ndarray.
            - ``"explained_variance_ratio"``: list of floats (one per component).
            - ``"sample_ids"`` (optional): list of str identifying samples.
        metadata: Optional mapping from sample id to a dict of annotations.
        output_file: Optional file path for saving the figure.
        components: Tuple of three 0-based PC indices to plot.
        color_by: Metadata field name used for colouring points.
        title: Title for the figure.

    Returns:
        Dictionary with ``"status"``, ``"output_path"``, ``"n_samples"``, and
        ``"n_components"`` keys.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    except ImportError:
        logger.warning("matplotlib/mpl_toolkits not available for 3D PCA plot")
        return {"status": "skipped", "output_path": None, "n_samples": 0, "n_components": 0}

    try:
        pcs = np.asarray(pca_data.get("pcs", []), dtype=np.float64)
        explained = pca_data.get("explained_variance_ratio", [])
        sample_ids = pca_data.get("sample_ids", [])

        if pcs.ndim != 2 or pcs.shape[0] == 0:
            logger.error("pca_data['pcs'] must be a non-empty 2D structure")
            return {"status": "failed", "output_path": None, "n_samples": 0, "n_components": 0}

        n_samples, n_pcs = pcs.shape

        # Validate that all requested components exist
        max_comp = max(components)
        if max_comp >= n_pcs:
            # Graceful degradation: fall back to available components
            available = tuple(c for c in components if c < n_pcs)
            if len(available) < 3:
                logger.warning(f"Only {n_pcs} PCs available, need 3 for 3D plot. " f"Returning status 'failed'.")
                return {
                    "status": "failed",
                    "output_path": None,
                    "n_samples": n_samples,
                    "n_components": n_pcs,
                    "error": f"Need at least 3 PCs but only {n_pcs} available",
                }
            components = (available[0], available[1], available[2])

        pc_a, pc_b, pc_c = components

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection="3d")

        # Resolve per-sample colour labels from metadata
        color_labels: Optional[List[str]] = None
        if metadata is not None and sample_ids:
            color_labels = [str(metadata.get(sid, {}).get(color_by, "unknown")) for sid in sample_ids]

        if color_labels is not None:
            unique_labels = sorted(set(color_labels))
            cmap = plt.cm.tab10(np.linspace(0, 1, max(len(unique_labels), 1)))
            label_to_color = {lab: cmap[i] for i, lab in enumerate(unique_labels)}
            for lab in unique_labels:
                mask = np.array([cl == lab for cl in color_labels])
                ax.scatter(
                    pcs[mask, pc_a],
                    pcs[mask, pc_b],
                    pcs[mask, pc_c],
                    label=lab,
                    alpha=0.7,
                    s=40,
                    c=[label_to_color[lab]],
                )
            ax.legend(title=color_by, fontsize=9)
        else:
            ax.scatter(pcs[:, pc_a], pcs[:, pc_b], pcs[:, pc_c], alpha=0.7, s=40, c="steelblue")

        # Axis labels with explained variance
        def _label(idx: int) -> str:
            base = f"PC{idx + 1}"
            if idx < len(explained):
                base += f" ({explained[idx] * 100:.1f}%)"
            return base

        ax.set_xlabel(_label(pc_a), fontsize=10)
        ax.set_ylabel(_label(pc_b), fontsize=10)
        ax.set_zlabel(_label(pc_c), fontsize=10)
        ax.set_title(title, fontsize=14, pad=20)

        output_path_str: Optional[str] = None
        if output_file:
            out = Path(output_file)
            out.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out, dpi=300, bbox_inches="tight")
            output_path_str = str(out)
            logger.info(f"Saved 3D PCA plot to {out}")

        plt.close(fig)
        return {
            "status": "success",
            "output_path": output_path_str,
            "n_samples": n_samples,
            "n_components": n_pcs,
        }

    except Exception as e:
        logger.error(f"Error creating 3D PCA plot: {e}")
        return {"status": "failed", "output_path": None, "n_samples": 0, "n_components": 0}
